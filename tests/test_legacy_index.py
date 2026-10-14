import json

import h5py
import numpy as np
import pytest

from malva._legacy_index import MalvaIndex


@pytest.fixture
def legacy_index_dir(tmp_path):
    index_dir = tmp_path / "legacy"
    index_dir.mkdir()
    index_file = index_dir / "malva_index.h5"

    with h5py.File(index_file, "w", driver="split") as f:
        f.attrs["kmer_size"] = 24
        f.attrs["coord_lims"] = (0.0, 1.0, 0.0, 1.0)
        f.attrs["n_spatial"] = 3
        f.attrs["n_chunks"] = 1
        f.attrs["project_id_shift"] = 22
        f.attrs["cell_id_mask"] = 4194303
        f.attrs["has_merged_projects"] = False
        f.attrs["project_mapping"] = json.dumps({"0": "projA"})
        f.create_dataset("index_0_data", data=np.arange(5))
        f.create_dataset("spatial_coord", data=np.zeros((3, 2), dtype=np.float32))

    return index_dir


class TestLegacyMalvaIndex:
    def test_open_reads_attrs_and_arrays(self, legacy_index_dir):
        index = MalvaIndex(str(legacy_index_dir))
        index.open()
        try:
            assert index.kmer_size == 24
            assert index.n_spatial == 3
            assert index.n_chunks == 1
            assert index.data_lengths == [5]
            assert index.coord_lims == (0.0, 1.0, 0.0, 1.0)
            assert index.project_mapping == {0: "projA"}
            assert index.PROJECT_ID_SHIFT == 22
            assert index.CELL_ID_MASK == 4194303
            assert index.spatial_coord is not None
        finally:
            index.close()

    def test_close_is_idempotent(self, legacy_index_dir):
        index = MalvaIndex(str(legacy_index_dir))
        index.open()
        index.close()
        index.close()
        assert index.index is None

    def test_defaults_before_open(self, legacy_index_dir):
        index = MalvaIndex(str(legacy_index_dir))
        assert index.kmer_size == 24
        assert index.n_spatial == 0
        assert index.n_chunks == 0
        assert index.index is None
