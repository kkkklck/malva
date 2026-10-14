import io

from malva.quant import write_mtx_header


class TestWriteMtxHeader:
    def test_writes_matrix_market_header_and_shape(self):
        buf = io.BytesIO()
        write_mtx_header(buf, (3, 5, 7))

        content = buf.getvalue().decode()
        lines = content.splitlines()
        assert lines[0] == "%%MatrixMarket matrix coordinate real general"
        assert lines[1] == "%"
        assert lines[2].split() == ["3", "5", "7"]


class TestProcessBatch:
    def test_writes_feature_names_and_counts(self, tmp_path):
        from malva.malva_index import MalvaIndex
        from malva.quant import process_batch
        from tests.index_helpers import OPENST_BARCODES, QUERY, build_index

        index_dir = build_index(tmp_path, OPENST_BARCODES, (2, 27), QUERY, spatial=True)
        index = MalvaIndex(str(index_dir))

        mtx = io.BytesIO()
        features = io.BytesIO()
        total_nnz = process_batch(
            index,
            [QUERY],
            ["mygene"],
            mtx,
            features,
            start_col=1,
            count_at_least=1,
            use_background_model=False,
        )

        assert isinstance(total_nnz, int)
        assert b"mygene" in features.getvalue()
