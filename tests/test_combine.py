import shutil

import pytest

from malva.combine import combine_indices
from malva.malva_index import MalvaIndex
from tests.index_helpers import OPENST_BARCODES, build_index

QUERY_B = "TGCATGCATGCATGCATGCATGCA"


class TestCombineIndices:
    def test_merges_two_indices_into_one(self, tmp_path):
        index_a = build_index(tmp_path / "A", OPENST_BARCODES, (2, 27), "ACGTACGTACGTACGTACGTACGT")
        index_b = build_index(tmp_path / "B", OPENST_BARCODES, (2, 27), QUERY_B)

        combine_dir = tmp_path / "combine"
        combine_dir.mkdir()
        shutil.copytree(index_a, combine_dir / "idxA")
        shutil.copytree(index_b, combine_dir / "idxB")

        project_mapping, n_indices = combine_indices(str(combine_dir))
        assert n_indices == 2
        assert set(project_mapping) == {0, 1}

        merged = MalvaIndex(str(combine_dir))
        merged.open()
        try:
            assert merged.n_chunks == 1
            assert (combine_dir / "pi.bin").exists()
            hits, _, _ = merged.where(
                "ACGTACGTACGTACGTACGTACGT",
                sliding_size=24,
                pct_threshold=0.65,
                count_at_least=1,
            )[0]
            assert hits.size >= 1
        finally:
            merged.close()

    def test_empty_dir_raises(self, tmp_path):
        empty = tmp_path / "empty"
        empty.mkdir()
        with pytest.raises(ValueError):
            combine_indices(str(empty))
