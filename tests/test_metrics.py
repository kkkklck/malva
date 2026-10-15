from malva.metrics import METRIC_PREFIX, IndexingMetrics, percentiles


class TestPercentiles:
    def test_empty(self):
        assert percentiles([]) == (0.0, 0.0, 0.0)

    def test_known_values(self):
        p5, p50, p95 = percentiles(list(range(1, 101)))
        assert round(p5) == 6
        assert round(p50) == 50
        assert round(p95) == 95


class TestIndexingMetrics:
    def test_derived_fields(self):
        metrics = IndexingMetrics(
            input_type="fastq",
            flavor="openst",
            kmer_size=24,
            l_prefix=12,
            jump_amount=24,
            overlap=False,
            version="1.1.2",
        )
        metrics.reads_total = 100
        metrics.reads_indexed = 80
        metrics.nucleotides_read = 1000
        metrics.nucleotides_indexed = 800
        metrics.read_seconds = 2.0
        metrics.throughput_samples = [100.0, 200.0, 300.0]

        data = metrics.as_dict()
        assert data["reads_unindexed"] == 20
        assert data["percent_reads_indexed"] == 80.0
        assert data["reads_per_second"] == 50.0
        assert data["flavor"] == "openst"
        assert data["malva_version"] == "1.1.2"

    def test_summary_lines_are_prefixed_key_value(self):
        metrics = IndexingMetrics(kmer_size=24, l_prefix=12)
        metrics.reads_total = 10
        metrics.reads_indexed = 10

        lines = metrics.summary_lines()
        assert all(line.startswith(METRIC_PREFIX + " ") for line in lines)
        assert f"{METRIC_PREFIX} reads_total=10" in lines
        assert f"{METRIC_PREFIX} reads_indexed=10" in lines

    def test_percent_indexed_zero_when_no_reads(self):
        assert IndexingMetrics().percent_indexed == 0.0
