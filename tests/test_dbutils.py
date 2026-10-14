import pytest

from malva.dbutils import (
    parse_multifasta,
    process_dna_string,
    process_ensembl_string,
    process_gene_string,
    validate_and_infer_query,
)


class TestProcessDnaString:
    def test_plain_sequence_is_uppercased(self):
        assert process_dna_string("acgt") == {"sequence": "ACGT"}

    def test_uracil_becomes_thymine(self):
        assert process_dna_string("ACGU") == {"sequence": "ACGT"}

    def test_invalid_characters_are_dropped(self):
        assert process_dna_string("AC-GT!") == {"sequence": "ACGT"}

    def test_whitespace_is_removed(self):
        assert process_dna_string("AC GT\nAC") == {"sequence": "ACGTAC"}

    def test_fasta_header_routes_to_multifasta(self):
        assert process_dna_string(">hdr\nACGT\n") == {"sequence": ["ACGT"]}


class TestParseMultifasta:
    def test_multiple_entries(self):
        assert parse_multifasta(">a\nACGT\n>b\nUUUU\n") == ["ACGT", "TTTT"]

    def test_no_entries(self):
        assert parse_multifasta("not a fasta") == []


class TestValidateAndInferQuery:
    def test_dna_sequence(self):
        assert validate_and_infer_query("ACGTN") == "ACGTN"

    def test_gene_id(self):
        assert validate_and_infer_query("TP53") == "gene:TP53"

    def test_invalid_input_raises(self):
        with pytest.raises(ValueError):
            validate_and_infer_query("ACGT?!")


class TestProcessGeneString:
    def test_defaults(self):
        assert process_gene_string("gene:TP53") == {
            "gene_id": "TP53",
            "species": "homo_sapiens",
            "split": [0, None],
            "seqtype": "genomic",
        }

    def test_all_fields(self):
        assert process_gene_string("gene:TP53;species:mus_musculus;type:cdna;split:2,10") == {
            "gene_id": "TP53",
            "species": "mus_musculus",
            "split": [2, 10],
            "seqtype": "cdna",
        }

    def test_requires_gene_prefix(self):
        with pytest.raises(ValueError):
            process_gene_string("TP53")

    def test_multiple_gene_ids_raise(self):
        with pytest.raises(ValueError):
            process_gene_string("gene:TP53;gene:KRAS")

    def test_bad_split_arity_raises(self):
        with pytest.raises(ValueError):
            process_gene_string("gene:TP53;split:1")

    def test_non_integer_split_raises(self):
        with pytest.raises(ValueError):
            process_gene_string("gene:TP53;split:a,b")


class TestProcessEnsemblString:
    def test_defaults(self):
        assert process_ensembl_string("ensembl:ENSG00000141510") == {
            "ensembl_id": "ENSG00000141510",
            "split": [0, None],
            "seqtype": "genomic",
        }

    def test_all_fields(self):
        assert process_ensembl_string("ensembl:ENSG00000141510;type:cdna;split:1,5") == {
            "ensembl_id": "ENSG00000141510",
            "split": [1, 5],
            "seqtype": "cdna",
        }

    def test_requires_ensembl_prefix(self):
        with pytest.raises(ValueError):
            process_ensembl_string("ENSG00000141510")

    def test_multiple_ids_raise(self):
        with pytest.raises(ValueError):
            process_ensembl_string("ensembl:ENSG1;ensembl:ENSG2")
