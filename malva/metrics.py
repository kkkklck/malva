# Copyright (c) 2025 Daniel León-Periñán and Nikolaos Karaiskos
#                    Rajewsky Lab, Max Delbrück Center for Molecular Medicine (MDC), Berlin
#
# Non-commercial and academic use only. See LICENSE for full terms.

"""Counters and a key=value summary for an indexing run.

The readers (FASTQ/SRA/BAM) fill an :class:`IndexingMetrics` as they go; the
indexing command emits :meth:`IndexingMetrics.summary_lines` just before the
success message and stores the counters (and the resolved flavour) in the
index's ``meta.json``.
"""

import numpy as np

# Every summary line starts with this marker so a log can be parsed with a
# single ``grep '^MALVA_INDEX_METRIC '`` without matching anything else.
METRIC_PREFIX = "MALVA_INDEX_METRIC"


def percentiles(values, points=(5, 50, 95)):
    """Return the requested percentiles of ``values`` (zeros if empty)."""
    if not values:
        return tuple(0.0 for _ in points)
    arr = np.asarray(values, dtype=float)
    return tuple(float(np.percentile(arr, p)) for p in points)


class IndexingMetrics:
    """Accumulates per-run counters and renders the final summary.

    Readers set the ``reads_*``, ``nucleotides_*``, ``kmer_cell_pairs``,
    ``chunks``, ``distinct_barcodes``, ``throughput_samples`` and ``*_seconds``
    fields; the command fills ``flavor``/config and the post-build index stats.
    """

    def __init__(
        self,
        input_type="unknown",
        flavor=None,
        kmer_size=None,
        l_prefix=None,
        jump_amount=None,
        overlap=False,
        version=None,
    ):
        self.input_type = input_type
        self.flavor = flavor
        self.kmer_size = kmer_size
        self.l_prefix = l_prefix
        self.jump_amount = jump_amount
        self.overlap = overlap
        self.version = version

        self.reads_total = 0
        self.reads_indexed = 0
        self.nucleotides_read = 0
        self.nucleotides_indexed = 0
        self.kmer_cell_pairs = 0
        self.chunks = 0
        self.distinct_barcodes = 0
        self.barcodes_total = 0
        # reads/s sampled once per reporting interval (for the percentiles)
        self.throughput_samples = []
        self.read_seconds = 0.0
        self.write_seconds = 0.0
        self.build_seconds = 0.0
        self.total_seconds = 0.0

        # filled after the index is written
        self.n_kmers = 0
        self.n_buckets = 0
        self.n_prefixes = 0
        self.n_cells = 0
        self.index_bytes = 0

        # the resolved flavour configuration (kept in meta.json)
        self.flavor_config = None

    @property
    def reads_unindexed(self):
        return max(self.reads_total - self.reads_indexed, 0)

    @property
    def percent_indexed(self):
        if self.reads_total == 0:
            return 0.0
        return 100.0 * self.reads_indexed / self.reads_total

    def as_dict(self):
        p5, p50, p95 = percentiles(self.throughput_samples)
        reads_per_s = self.reads_total / self.read_seconds if self.read_seconds else 0.0
        nt_per_s = self.nucleotides_read / self.read_seconds if self.read_seconds else 0.0
        return {
            "input_type": self.input_type,
            "flavor": self.flavor,
            "malva_version": self.version,
            "kmer_size": self.kmer_size,
            "l_prefix": self.l_prefix,
            "jump_amount": self.jump_amount,
            "overlap": int(bool(self.overlap)),
            "reads_total": self.reads_total,
            "reads_indexed": self.reads_indexed,
            "reads_unindexed": self.reads_unindexed,
            "percent_reads_indexed": round(self.percent_indexed, 2),
            "nucleotides_read": self.nucleotides_read,
            "nucleotides_indexed": self.nucleotides_indexed,
            "kmer_cell_pairs": self.kmer_cell_pairs,
            "chunks": self.chunks,
            "distinct_barcodes": self.distinct_barcodes,
            "barcodes_total": self.barcodes_total,
            "reads_per_second": round(reads_per_s, 2),
            "nucleotides_per_second": round(nt_per_s, 2),
            "reads_per_second_p5": round(p5, 2),
            "reads_per_second_p50": round(p50, 2),
            "reads_per_second_p95": round(p95, 2),
            "read_seconds": round(self.read_seconds, 3),
            "write_seconds": round(self.write_seconds, 3),
            "build_seconds": round(self.build_seconds, 3),
            "total_seconds": round(self.total_seconds, 3),
            "n_kmers": self.n_kmers,
            "n_buckets": self.n_buckets,
            "n_prefixes": self.n_prefixes,
            "n_cells": self.n_cells,
            "index_bytes": self.index_bytes,
        }

    def summary_lines(self):
        """One ``MALVA_INDEX_METRIC key=value`` line per metric."""
        return [f"{METRIC_PREFIX} {key}={value}" for key, value in self.as_dict().items()]
