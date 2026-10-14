"""Tests for the prefix-bucketed index (``malva.malva_index``) and spatial barcodes.

The index was refactored from the legacy HDF5 layout to the prefix-bucketed
``pi.bin`` / ``suffixes.bin`` / ``data.bin`` layout, with ``MalvaIndex`` moved to
``malva.malva_index`` and the barcode index to ``malva.barcodes``. These tests
exercise the current API end to end on small, deterministic datasets.
"""

import gzip

import numpy as np

from malva.barcodes import (
    SpatialIndex,
    create_singlecell_index,
    create_spatial_index,
)
from malva.malva_index import MalvaIndex

# 24-mer query that every crafted R2 read carries.
QUERY = "ACGTACGTACGTACGTACGTACGT"
# A 24-mer absent from every crafted dataset (exercises the empty-result path).
ABSENT = "TTTTTTTTTTTTTTTTTTTTTTTT"

OPENST_BARCODES = [
    "ACGTACGTACGTACGTACGTACGTA",
    "TGCATGCATGCATGCATGCATGCAT",
    "GATTACAGATTACAGATTACAGATT",
    "CCCCGGGGAAAATTTTCCCCGGGGA",
    "TTTTAAAACCCCGGGGTTTTAAAAC",
]
SC_BARCODES = [
    "ACGTACGTACGTACGT",
    "TGCATGCATGCATGCA",
    "GATTACAGATTACAGA",
    "CCCCGGGGAAAATTTT",
]


def _write_fastq_pair(path, records):
    """Write one paired FASTQ pair from ``(name, r1_seq, r2_seq)`` records."""
    r1_path = path / "R1.fastq.gz"
    r2_path = path / "R2.fastq.gz"
    with gzip.open(r1_path, "wt") as f1, gzip.open(r2_path, "wt") as f2:
        for name, r1_seq, r2_seq in records:
            f1.write(f"{name}\n{r1_seq}\n+\n{'I' * len(r1_seq)}\n")
            f2.write(f"{name}\n{r2_seq}\n+\n{'I' * len(r2_seq)}\n")
    return r1_path, r2_path


def _build_index(path, barcodes, barcode_span, r2_sequence, spatial=True):
    """Build a malva index over a crafted dataset and return its directory."""
    reads = []
    r1_len = barcode_span[1] + 10
    for i, barcode in enumerate(barcodes):
        pad = "A" * (r1_len - barcode_span[1] - barcode_span[0])
        r1_seq = "N" * barcode_span[0] + barcode + pad
        reads.append((f"@read{i}", r1_seq, r2_sequence))
    r1_path, r2_path = _write_fastq_pair(path, reads)

    if spatial:
        spatial_file = path / "spatial.csv"
        with open(spatial_file, "w") as f:
            f.write("barcode,x,y\n")
            for i, barcode in enumerate(barcodes):
                f.write(f"{barcode},{i * 10.0},{i * 5.0}\n")
        sindex = create_spatial_index(str(spatial_file))
    else:
        whitelist = path / "whitelist.txt"
        whitelist.write_text("\n".join(barcodes) + "\n")
        sindex = create_singlecell_index(str(whitelist))

    index_dir = path / "index"
    index_dir.mkdir()
    index = MalvaIndex(str(index_dir), kmer_size_initialize=24)
    if spatial:
        index.set_spatial_index(sindex)
    else:
        index.set_barcode_index(sindex)
    index.add_reads(
        [str(r1_path), str(r2_path)],
        cell=f"r1[{barcode_span[0]}:{barcode_span[1]}]",
        chunksize=10,
    )
    index.close()
    return index_dir


class TestSpatialIndex:
    def test_spatial_index_from_csv(self, tmp_path):
        spatial_file = tmp_path / "spatial.csv"
        with open(spatial_file, "w") as f:
            f.write("barcode,x,y\n")
            for i, barcode in enumerate(OPENST_BARCODES):
                f.write(f"{barcode},{i * 10.0},{i * 5.0}\n")

        sindex = create_spatial_index(str(spatial_file))
        assert isinstance(sindex, SpatialIndex)
        assert sindex.num_items() == len(OPENST_BARCODES)
        coords = sindex.get_coords()
        assert coords.shape == (len(OPENST_BARCODES), 2)
        assert coords.dtype == np.float32

    def test_singlecell_index_from_whitelist(self, tmp_path):
        whitelist = tmp_path / "whitelist.txt"
        whitelist.write_text("\n".join(SC_BARCODES) + "\n")
        sindex = create_singlecell_index(str(whitelist))
        assert sindex.num_items() == len(SC_BARCODES)


class TestIndexBuildAndQuery:
    def test_build_reopen_and_query_spatial(self, tmp_path):
        index_dir = _build_index(tmp_path, OPENST_BARCODES, (2, 27), QUERY, spatial=True)

        index = MalvaIndex(str(index_dir))
        index.open()
        try:
            assert index.n_spatial == len(OPENST_BARCODES)
            assert "spatial_coord" in index.index

            results = index.where(QUERY, sliding_size=24, pct_threshold=0.65, count_at_least=1)
            # One sequence -> one group; a group is (locations, counts, matches).
            assert len(results) == 1
            locations, counts, matches = results[0]
            assert isinstance(locations, np.ndarray) and locations.dtype == np.uint32
            assert isinstance(counts, np.ndarray) and counts.dtype == np.uint32
            assert isinstance(matches, list)

            # Every crafted read carries QUERY. Cell id 0 is the "not found"
            # sentinel in the barcode map, so the barcodes that survive are 1..n-1.
            assert set(locations.tolist()) == set(range(1, len(OPENST_BARCODES)))
        finally:
            index.close()

    def test_absent_sequence_returns_empty(self, tmp_path):
        index_dir = _build_index(tmp_path, OPENST_BARCODES, (2, 27), QUERY, spatial=True)
        index = MalvaIndex(str(index_dir))
        index.open()
        try:
            locations, counts, _ = index.where(
                ABSENT, sliding_size=24, pct_threshold=0.65, count_at_least=1
            )[0]
            assert locations.size == 0
            assert counts.size == 0
        finally:
            index.close()

    def test_absent_kmers_after_present_kmers_do_not_crash(self, tmp_path):
        """Regression: a batch with absent kmers after a present batch must not
        dereference the stale CSR buffers left by the previous query."""
        index_dir = _build_index(tmp_path, OPENST_BARCODES, (2, 27), QUERY, spatial=True)
        index = MalvaIndex(str(index_dir))
        index.open()
        try:
            present = index.where(QUERY, sliding_size=24, pct_threshold=0.65, count_at_least=1)
            assert present[0][0].size > 0

            absent = index.where(ABSENT, sliding_size=24, pct_threshold=0.65, count_at_least=1)
            assert absent[0][0].size == 0

            # A mixed batch keeps the present hits and drops the absent ones.
            mixed = index.where(
                [QUERY, ABSENT], sliding_size=24, pct_threshold=0.65, count_at_least=1
            )
            assert len(mixed) == 1
        finally:
            index.close()

    def test_where_groups_list_of_sequences_as_one_group(self, tmp_path):
        index_dir = _build_index(tmp_path, OPENST_BARCODES, (2, 27), QUERY, spatial=True)
        index = MalvaIndex(str(index_dir))
        index.open()
        try:
            # A flat list is a single group; a list of lists is one group each.
            assert len(index.where([QUERY, ABSENT], sliding_size=24)) == 1
            assert len(index.where([[QUERY], [ABSENT]], sliding_size=24)) == 2
        finally:
            index.close()

    def test_build_reopen_and_query_singlecell(self, tmp_path):
        index_dir = _build_index(tmp_path, SC_BARCODES, (0, 16), QUERY, spatial=False)
        index = MalvaIndex(str(index_dir))
        index.open()
        try:
            assert index.n_spatial == len(SC_BARCODES)
            locations, _, _ = index.where(
                QUERY, sliding_size=24, pct_threshold=0.65, count_at_least=1
            )[0]
            assert locations.size >= 1
            # Single-cell ids are 1-based (0 is the not-found sentinel).
            assert np.all(locations >= 1)
            assert np.all(locations <= len(SC_BARCODES))
        finally:
            index.close()
