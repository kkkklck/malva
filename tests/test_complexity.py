from malva.complexity import (
    complexity,
    compute_rep_vector,
    mask_low_complexity,
    overlapping_windows,
)


class TestOverlappingWindows:
    def test_returns_all_windows(self):
        assert overlapping_windows("ACGTACG", 4) == ["ACGT", "CGTA", "GTAC", "TACG"]

    def test_window_equal_to_sequence(self):
        assert overlapping_windows("ACGT", 4) == ["ACGT"]

    def test_window_longer_than_sequence(self):
        assert overlapping_windows("ACGT", 5) == []


class TestComputeRepVector:
    def test_pads_and_sorts_descending(self):
        # Only A and G are seen, so the vector is padded to length N and sorted.
        assert compute_rep_vector("AAGG", 4) == [2, 2, 0, 0]

    def test_stops_once_n_residues_are_seen(self):
        # A (3) and G (2) are seen; with N=2 the later C run is never counted.
        assert compute_rep_vector("AAAGGCCCC", 2) == [3, 2]

    def test_homopolymer(self):
        assert compute_rep_vector("AAAA", 4) == [4, 0, 0, 0]


class TestComplexity:
    def test_homopolymer_has_zero_entropy(self):
        assert complexity("AAAAAAAA", 4) == 0

    def test_diverse_sequence_has_higher_entropy(self):
        assert complexity("AAAAAAAA", 4) < complexity("ACGTACGT", 4)


class TestMaskLowComplexity:
    def test_preserves_length_and_type(self):
        seq = "ACGT" * 10
        masked = mask_low_complexity(seq, N=4, L=12)
        assert isinstance(masked, str)
        assert len(masked) == len(seq)

    def test_masks_a_long_homopolymer_run(self):
        seq = "ACGTACGTACGT" + "A" * 30 + "ACGTACGTACGT"
        masked = mask_low_complexity(seq, maskchar="N", N=4, L=12)
        assert "N" in masked
        assert len(masked) == len(seq)
