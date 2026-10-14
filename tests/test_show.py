import numpy as np
import pytest

from malva.malva_index import MalvaIndex
from malva.show import MalvaPlot
from tests.index_helpers import OPENST_BARCODES, QUERY, build_index


class TestMalvaPlot:
    def test_rejects_non_index(self):
        with pytest.raises(ValueError):
            MalvaPlot(object())

    @pytest.fixture
    def plotter(self, tmp_path):
        index_dir = build_index(tmp_path, OPENST_BARCODES, (2, 27), QUERY, spatial=True)
        return MalvaPlot(MalvaIndex(str(index_dir)))

    def test_loads_index_geometry(self, plotter):
        assert plotter.n_spatial == len(OPENST_BARCODES)
        assert plotter.xmax > plotter.xmin
        assert plotter.ymax > plotter.ymin

    def test_image_without_smoothing(self, plotter):
        image = plotter.image(np.array([1, 2]), np.array([1, 1]), render_smoothing=0)
        assert image.ndim == 2

    def test_image_with_smoothing_and_normalize(self, plotter):
        image = plotter.image(np.array([1, 2]), np.array([2, 2]), normalize=True)
        assert image.dtype == np.uint8

    def test_image_normalize_all_zero(self, plotter):
        image = plotter.image(np.array([], dtype=int), np.array([], dtype=int), normalize=True)
        assert image.dtype == np.uint8
