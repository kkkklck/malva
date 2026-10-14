from malva.kmer_processing import (
    FastKmerExtractor,
    decode_kmer,
    encode_kmer,
    get_kmers_numeric,
    get_sliding_kmers_numeric,
)


class TestKmerProcessing:
    def test_encode_decode_kmer(self):
        test_kmers = ["ACGT", "AAAA", "CCCC", "GGGG", "TTTT"]
        for kmer in test_kmers:
            encoded = encode_kmer(kmer)
            decoded = decode_kmer(encoded, len(kmer))
            assert decoded == kmer

    def test_get_kmers_numeric(self):
        sequence = "ACGTACGT"
        k = 4
        kmers = get_kmers_numeric(sequence, k)
        assert len(kmers) == 2

        # Test with non-complex sequences
        sequence_with_n = "ACGTNACGT"
        kmers = get_kmers_numeric(sequence_with_n, k, remove_noncomplex=True)
        assert 0 in kmers

    def test_get_sliding_kmers_numeric(self):
        sequence = "ACGTACGT"
        k = 4
        kmers = get_sliding_kmers_numeric(sequence, k)
        assert len(kmers) == 5


class TestKmerProcessingExtra:
    def test_encode_special_bases(self):
        # U encodes like T, N encodes to 3.
        assert encode_kmer("U") == encode_kmer("T")
        assert encode_kmer("N") == 3

    def test_decode_zero_is_poly_a(self):
        assert decode_kmer(0, 4) == "AAAA"

    def test_get_kmers_numeric_pads_partial_trailing_kmer(self):
        # 6 bases with k=4 -> one full k-mer and one built from overlap.
        assert len(get_kmers_numeric("ACGTAC", 4)) == 2

    def test_get_kmers_numeric_without_removing_n(self):
        kmers = get_kmers_numeric("ACGTNACGT", 4, remove_noncomplex=False)
        assert 0 not in kmers

    def test_get_kmers_numeric_remove_noncomplex_marks_n(self):
        kmers = get_kmers_numeric("ACGTNACGT", 4, remove_noncomplex=True)
        assert 0 in kmers

    def test_sliding_kmers_shorter_than_k(self):
        assert get_sliding_kmers_numeric("ACGT", 8) == []

    def test_sliding_kmers_remove_noncomplex(self):
        kmers = get_sliding_kmers_numeric("ACGTN", 4, remove_noncomplex=True)
        assert kmers == [encode_kmer("ACGT"), 0]

    def test_fast_kmer_extractor_constructs(self):
        assert FastKmerExtractor(4, True) is not None
