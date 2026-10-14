import json

import h5py
import numpy as np

from malva.convert import convert_malva_to_prefix, verify_indices
from malva.indexes import PrefixIndex
from malva.kmer_processing import encode_kmer

KMER_A = encode_kmer("ACGTACGTACGTACGTACGTACGT")
KMER_B = encode_kmer("TGCATGCATGCATGCATGCATGCA")


def _write_legacy_index(index_dir, kmers, indptr, cells, kmer_size=24, n_spatial=4):
    index_dir.mkdir(parents=True, exist_ok=True)
    with h5py.File(index_dir / "malva_index.h5", "w", driver="split") as f:
        f.attrs["kmer_size"] = kmer_size
        f.attrs["n_spatial"] = n_spatial
        f.attrs["n_chunks"] = 1
        f.create_dataset("index_0_indices", data=np.asarray(sorted(kmers), dtype=np.uint64))
        f.create_dataset("index_0_indptr", data=np.asarray(indptr, dtype=np.uint64))
        f.create_dataset("index_0_data", data=np.asarray(cells, dtype=np.uint32))
        f.create_dataset("spatial_coord", data=np.zeros((n_spatial, 2), dtype=np.float32))


class TestConvertMalvaToPrefix:
    def test_convert_writes_a_queryable_prefix_index(self, tmp_path):
        legacy = tmp_path / "legacy"
        _write_legacy_index(legacy, [KMER_A, KMER_B], [0, 2, 4], [1, 2, 3, 4])
        out = tmp_path / "prefix"

        convert_malva_to_prefix(str(legacy), str(out), l_prefix=8, verbose=False)

        assert (out / "pi.bin").exists()
        assert (out / "spatial_coord.npy").exists()
        meta = json.loads((out / "meta.json").read_text())
        assert meta["n_cells"] == 4

        index = PrefixIndex()
        index.open(str(out))
        try:
            result = index.query(np.array([KMER_A, KMER_B], dtype=np.uint64), count_at_least=1)
            assert isinstance(result, dict)
        finally:
            index.close()


class TestVerifyIndices:
    def test_verify_reports_a_match(self, tmp_path):
        legacy = tmp_path / "legacy"
        _write_legacy_index(legacy, [KMER_A, KMER_B], [0, 2, 4], [1, 2, 3, 4])
        out = tmp_path / "prefix"
        convert_malva_to_prefix(str(legacy), str(out), l_prefix=8, verbose=False)

        assert verify_indices(str(legacy), str(out), verbose=False) is True
