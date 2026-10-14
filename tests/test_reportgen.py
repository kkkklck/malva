import io
import zipfile

from malva.serve.reportgen import HTMLReportGenerator


class TestHTMLReportGenerator:
    def test_report_zip_contains_all_components(self):
        generator = HTMLReportGenerator(
            {
                "query_term": "gene:TP53",
                "sliding_size": 128,
                "pct_threshold": 0.65,
                "low_complexity_filter": True,
                "query_seq": "ACGTACGT",
                "where_abundant": True,
            }
        )
        data = generator.create_report_zip()

        with zipfile.ZipFile(io.BytesIO(data)) as archive:
            names = set(archive.namelist())
            assert {"report.html", "queries.fa", "malva.log", "summary.txt"} <= names
            html = archive.read("report.html").decode()
            assert "Malva Search Report" in html
            assert "gene:TP53" in html
            fasta = archive.read("queries.fa").decode()
            assert "ACGTACGT" in fasta

    def test_sequences_as_list_of_dicts(self):
        generator = HTMLReportGenerator(
            {
                "query_seq": [
                    {"header": ">a", "sequence": "ACGT"},
                    {"header": ">b", "sequence": "TTTT"},
                ]
            }
        )
        fasta = generator._create_fasta()
        assert ">a\nACGT" in fasta
        assert ">b\nTTTT" in fasta

    def test_long_sequence_is_wrapped_in_html(self):
        generator = HTMLReportGenerator({"query_seq": "A" * 130, "query_term": "probe"})
        html = generator._create_html()
        assert ">query_sequence_1" in html
        assert len(generator._create_fasta().splitlines()) == 2

    def test_summary_counts_sequences(self):
        generator = HTMLReportGenerator({"query_term": "x", "query_seq": ["ACGT"]})
        summary = generator._create_summary()
        assert "MALVA SEARCH SUMMARY" in summary
        assert "Number of Sequences: 1" in summary
        assert "Query Term: x" in summary
