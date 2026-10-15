import gzip
import json

import pytest

from malva.cli import get_index_parser
from malva.index import _run_index
from malva.malva_index import MalvaIndex
from malva.utils import FormatError
from tests.index_helpers import OPENST_BARCODES, QUERY

# cDNA long enough that the openst [2:27] slice is always in range.
CDNA = QUERY + "A" * 12


def _write_r2_barcode_dataset(path, barcodes, cdna=CDNA, span=(2, 27)):
    """R1 carries the cDNA, R2 carries the barcode at ``span``."""
    r1 = path / "R1.fastq.gz"
    r2 = path / "R2.fastq.gz"
    with gzip.open(r1, "wt") as f1, gzip.open(r2, "wt") as f2:
        for i, barcode in enumerate(barcodes):
            r2_seq = "N" * span[0] + barcode + "A" * 10
            f1.write(f"@read{i}\n{cdna}\n+\n{'I' * len(cdna)}\n")
            f2.write(f"@read{i}\n{r2_seq}\n+\n{'I' * len(r2_seq)}\n")
    spatial = path / "spatial.csv"
    spatial.write_text(
        "barcode,x,y\n" + "".join(f"{b},{i * 10.0},{i * 5.0}\n" for i, b in enumerate(barcodes))
    )
    return r1, r2, spatial


def _run(reads_in, spatial, index_out, *extra):
    args = get_index_parser().parse_args(
        [
            "--reads-in",
            *[str(p) for p in reads_in],
            "--spatial-bc-in",
            str(spatial),
            "--index-out",
            str(index_out),
            # small k -> small l_prefix -> a tiny pi.bin (2^(2*l_prefix) entries)
            "--kmer-length",
            "12",
            *extra,
        ]
    )
    _run_index(args)


class TestCellOverride:
    def test_cell_r2_indexes_barcodes_from_the_second_read(self, tmp_path):
        r1, r2, spatial = _write_r2_barcode_dataset(tmp_path, OPENST_BARCODES)
        index_out = tmp_path / "index"

        _run((r1, r2), spatial, index_out, "--cell", "r2[2:27]")

        meta = json.loads((index_out / "meta.json").read_text())
        assert meta["flavor"]["cell"] == "r2[2:27]"
        assert meta["metrics"]["reads_indexed"] > 0
        assert meta["metrics"]["kmer_cell_pairs"] > 0

        index = MalvaIndex(str(index_out))
        index.open()
        try:
            assert index.n_spatial == len(OPENST_BARCODES)
            hits, _, _ = index.where(CDNA, sliding_size=12, pct_threshold=0.65, count_at_least=1)[0]
            assert hits.size >= 1
        finally:
            index.close()

    def test_default_flavor_does_not_find_barcodes_in_read_two(self, tmp_path):
        r1, r2, spatial = _write_r2_barcode_dataset(tmp_path, OPENST_BARCODES)
        index_out = tmp_path / "index"

        _run((r1, r2), spatial, index_out)  # openst reads the barcode from r1[2:27]

        meta = json.loads((index_out / "meta.json").read_text())
        assert meta["metrics"]["reads_indexed"] == 0

    def test_invalid_cell_spec_is_rejected(self, tmp_path):
        r1, r2, spatial = _write_r2_barcode_dataset(tmp_path, OPENST_BARCODES)
        with pytest.raises(FormatError):
            _run((r1, r2), spatial, tmp_path / "index", "--cell", "not-a-spec")


class TestMetricsPersistence:
    def test_flavour_and_counters_are_stored_in_meta(self, tmp_path):
        r1, r2, spatial = _write_r2_barcode_dataset(tmp_path, OPENST_BARCODES)
        index_out = tmp_path / "index"

        _run((r1, r2), spatial, index_out, "--cell", "r2[2:27]")

        meta = json.loads((index_out / "meta.json").read_text())
        assert meta["flavor"]["cell"] == "r2[2:27]"
        assert meta["flavor"]["bam_tags"] == "CB:{cell}"
        metrics = meta["metrics"]
        for key in (
            "reads_total",
            "reads_indexed",
            "percent_reads_indexed",
            "nucleotides_read",
            "nucleotides_indexed",
            "kmer_cell_pairs",
            "distinct_barcodes",
            "barcodes_total",
            "reads_per_second_p95",
            "total_seconds",
            "n_kmers",
        ):
            assert key in metrics
        assert metrics["reads_total"] == len(OPENST_BARCODES)
        assert metrics["nucleotides_read"] == len(OPENST_BARCODES) * len(CDNA)
