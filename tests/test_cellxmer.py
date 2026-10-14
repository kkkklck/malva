import numpy as np

from malva.cellxmer import _get_csr_arrays, _get_csr_arrays_chunk, malva_to_cellxmer
from malva.malva_index import MalvaIndex
from tests.index_helpers import OPENST_BARCODES, QUERY, build_index


class TestCsrArrays:
    def test_full_arrays(self, tmp_path):
        index_dir = build_index(tmp_path, OPENST_BARCODES, (2, 27), QUERY, spatial=True)
        index = MalvaIndex(str(index_dir))
        indices, indptr, data, num_cells = _get_csr_arrays(index, verbose=False)

        assert num_cells == len(OPENST_BARCODES)
        assert indptr.dtype == np.int64
        assert len(indptr) >= 1

    def test_chunk_is_a_prefix_of_the_whole(self, tmp_path):
        index_dir = build_index(tmp_path, OPENST_BARCODES, (2, 27), QUERY, spatial=True)
        index = MalvaIndex(str(index_dir))
        indices, indptr, data, num_cells, total = _get_csr_arrays_chunk(index, 0, 2, verbose=False)

        assert num_cells == len(OPENST_BARCODES)
        assert total >= len(indices)
        assert indptr.dtype == np.int64

    def test_malva_to_cellxmer_builds_anndata(self, tmp_path):
        import anndata

        index_dir = build_index(tmp_path, OPENST_BARCODES, (2, 27), QUERY, spatial=True)
        index = MalvaIndex(str(index_dir))
        adata = malva_to_cellxmer(index, verbose=False)

        assert isinstance(adata, anndata.AnnData)
        assert adata.n_obs == len(OPENST_BARCODES)
