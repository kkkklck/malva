"""Shared fixtures for building small, deterministic malva indices in tests."""

import gzip

from malva.barcodes import create_singlecell_index, create_spatial_index
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


def write_fastq_pair(path, records):
    """Write one paired FASTQ pair from ``(name, r1_seq, r2_seq)`` records."""
    r1_path = path / "R1.fastq.gz"
    r2_path = path / "R2.fastq.gz"
    with gzip.open(r1_path, "wt") as f1, gzip.open(r2_path, "wt") as f2:
        for name, r1_seq, r2_seq in records:
            f1.write(f"{name}\n{r1_seq}\n+\n{'I' * len(r1_seq)}\n")
            f2.write(f"{name}\n{r2_seq}\n+\n{'I' * len(r2_seq)}\n")
    return r1_path, r2_path


def build_index(path, barcodes, barcode_span, r2_sequence, spatial=True):
    """Build a malva index over a crafted dataset and return its directory."""
    path.mkdir(parents=True, exist_ok=True)
    reads = []
    r1_len = barcode_span[1] + 10
    for i, barcode in enumerate(barcodes):
        pad = "A" * (r1_len - barcode_span[1] - barcode_span[0])
        r1_seq = "N" * barcode_span[0] + barcode + pad
        reads.append((f"@read{i}", r1_seq, r2_sequence))
    r1_path, r2_path = write_fastq_pair(path, reads)

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
    # A small l_prefix keeps pi.bin tiny; the test only needs a working index.
    index = MalvaIndex(str(index_dir), kmer_size_initialize=24, l_prefix=8)
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
