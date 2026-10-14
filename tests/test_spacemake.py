import anndata
import numpy as np
import pytest

from malva.spacemake import (
    aggregate_adata_by_indices,
    binning_hexagon,
    create_mesh,
    create_meshed_adata,
    nonsingular,
)


class TestNonsingular:
    def test_passthrough(self):
        assert nonsingular(1.0, 2.0) == (1.0, 2.0)

    def test_swaps_when_increasing(self):
        assert nonsingular(2.0, 1.0) == (1.0, 2.0)

    def test_no_swap_when_not_increasing(self):
        assert nonsingular(2.0, 1.0, increasing=False) == (2.0, 1.0)

    def test_non_finite_returns_expanded_interval(self):
        assert nonsingular(float("nan"), 1.0) == (-0.001, 0.001)
        assert nonsingular(0.0, float("inf")) == (-0.001, 0.001)

    def test_zeros_are_expanded(self):
        assert nonsingular(0.0, 0.0) == (-0.001, 0.001)


class TestCreateMesh:
    def test_shape_and_origin(self):
        xy = create_mesh(width=10, height=10, diameter=1, distance=2)
        assert xy.ndim == 2
        assert xy.shape[1] == 2
        assert tuple(xy[0]) == (0.0, 0.0)

    def test_push_offsets_the_origin(self):
        xy = create_mesh(width=10, height=10, diameter=1, distance=2, push_x=1, push_y=2)
        assert tuple(xy[0]) == (1.0, 2.0)


class TestBinningHexagon:
    def test_returns_coordinates_and_bins(self):
        coords, accumulated = binning_hexagon(
            [0.5, 1.0, 1.5], [0.5, 1.0, 1.5], gridsize=2, extent=(0, 2, 0, 2)
        )
        assert coords.ndim == 2 and coords.shape[1] == 2
        assert isinstance(accumulated, list)
        binned = [i for bucket in accumulated for i in bucket]
        assert set(binned).issubset({0, 1, 2})
        assert len(binned) >= 1

    def test_accepts_tuple_gridsize(self):
        coords, accumulated = binning_hexagon([0.5], [0.5], gridsize=(2, 1), extent=(0, 2, 0, 2))
        assert coords.ndim == 2


class TestAggregateAdataByIndices:
    def test_sums_observations_into_bins(self):
        adata = anndata.AnnData(X=np.arange(8, dtype=float).reshape(4, 2))
        adata.obsm["spatial"] = np.array([[0, 0], [1, 1], [2, 2], [3, 3]], dtype=float)

        aggregated = aggregate_adata_by_indices(
            adata,
            idx_to_aggregate=np.array([0, 1, 2, 3]),
            idx_aggregated=np.array([0, 0, 1, 1]),
            coordinates_aggregated=np.array([[0.0, 0.0], [1.0, 1.0]]),
        )

        assert aggregated.shape == (2, 2)
        np.testing.assert_array_equal(aggregated.X, [[2, 4], [10, 12]])
        assert list(aggregated.obs["n_joined"]) == [2, 2]


class TestCreateMeshedAdata:
    @pytest.fixture
    def spatial_adata(self):
        coords = np.array([[x * 10.0, y * 10.0] for x in range(10) for y in range(10)])
        adata = anndata.AnnData(X=np.ones((len(coords), 3), dtype=float))
        adata.obsm["spatial"] = coords
        return adata

    def test_hexagon_mesh(self, spatial_adata):
        meshed = create_meshed_adata(spatial_adata, px_by_um=1.0, mesh_type="hexagon")
        assert isinstance(meshed, anndata.AnnData)

    def test_rejects_unknown_mesh_type(self, spatial_adata):
        with pytest.raises(ValueError):
            create_meshed_adata(spatial_adata, px_by_um=1.0, mesh_type="triangle")
