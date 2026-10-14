import pytest

from malva.serve.modeling import parse_gene_lists


class TestParseGeneLists:
    def test_single_list(self):
        assert parse_gene_lists('["TP53", "KRAS"]') == (["TP53", "KRAS"], None)

    def test_two_lists(self):
        assert parse_gene_lists('["A", "B"]\n["C"]') == (["A", "B"], ["C"])

    def test_list_embedded_in_prose(self):
        assert parse_gene_lists('Here you go: ["EGFR"] done') == (["EGFR"], None)

    def test_numbers_are_parsed(self):
        assert parse_gene_lists("[1, 2, 3]") == ([1, 2, 3], None)

    def test_missing_list_raises(self):
        with pytest.raises(ValueError):
            parse_gene_lists("there is no list in here")
