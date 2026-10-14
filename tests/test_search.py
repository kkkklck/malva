from types import SimpleNamespace
from unittest.mock import MagicMock, patch

import requests

from malva.search import (
    _run_job_status,
    _run_list_datasets,
    _run_search,
    _run_search_data,
    check_server_connection,
    format_results,
    get_datasets,
    get_job_status,
    read_sequences_file,
    submit_search_job,
    wait_for_completion,
)


class TestReadSequencesFile:
    def test_reads_fasta_sequences(self, tmp_path):
        fasta = tmp_path / "seqs.fa"
        fasta.write_text(">a\nACGT\n>b\nTTTT\n")
        assert read_sequences_file(str(fasta)) == ["ACGT", "TTTT"]


class TestFormatResults:
    BASE = {"job_id": "j1", "dataset_id": "d1", "status": "completed"}

    def test_json(self):
        assert '"job_id": "j1"' in format_results(self.BASE, "json")

    def test_text_completed_with_metadata(self):
        results = dict(
            self.BASE,
            created_at=1.0,
            completed_at=2.5,
            query={"query_type": "gene"},
            results={
                "ACGT": {
                    "cells": [0, 1, 2, 3, 4, 5, 6],
                    "counts": [10, 9, 8, 7, 6, 5, 4],
                    "metadata": {0: {"cell_type": "T"}},
                }
            },
        )
        text = format_results(results, "text")
        assert "Job ID: j1" in text
        assert "Execution time: 1.50 seconds" in text
        assert "cell_type: T" in text
        assert "and 2 more cells" in text

    def test_text_error(self):
        text = format_results(dict(self.BASE, status="error", error="boom"), "text")
        assert "Error: boom" in text


class TestServerCalls:
    def test_check_server_connection_ok(self):
        with patch("malva.search.requests.get") as get:
            get.return_value = SimpleNamespace(status_code=200)
            assert check_server_connection("http://x") is True

    def test_check_server_connection_failure(self):
        with patch("malva.search.requests.get", side_effect=requests.RequestException):
            assert check_server_connection("http://x") is False

    def test_submit_search_job(self):
        with patch("malva.search.requests.post") as post:
            post.return_value = MagicMock(json=lambda: {"job_id": "j1"})
            result = submit_search_job("http://x", "ACGT", dataset_id="d1", include_metadata=True)
        assert result == {"job_id": "j1"}

    def test_get_datasets_and_status(self):
        with patch("malva.search.requests.get") as get:
            get.return_value = MagicMock(json=lambda: [{"dataset_id": "d1"}])
            assert get_datasets("http://x") == [{"dataset_id": "d1"}]
            get.return_value = MagicMock(json=lambda: {"status": "completed"})
            assert get_job_status("http://x", "j1") == {"status": "completed"}


class TestWaitForCompletion:
    def test_returns_when_completed(self):
        with patch("malva.search.get_job_status", return_value={"status": "completed"}):
            assert wait_for_completion("http://x", "j1", poll_interval=0) == {"status": "completed"}


class TestRunCommands:
    def test_list_datasets(self, capsys):
        datasets = [
            {
                "dataset_id": "d1",
                "name": "n",
                "description": "d",
                "organism": "human",
                "cell_count": 10,
                "kmer_size": 24,
            }
        ]
        args = SimpleNamespace(server="http://x")
        with (
            patch("malva.search.check_server_connection", return_value=True),
            patch("malva.search.get_datasets", return_value=datasets),
        ):
            assert _run_list_datasets(args) == 0
        assert "Available datasets (1)" in capsys.readouterr().out

    def test_job_status_to_stdout(self, capsys):
        args = SimpleNamespace(server="http://x", job_id="j1", output=None, format="text")
        with (
            patch("malva.search.check_server_connection", return_value=True),
            patch(
                "malva.search.get_job_status",
                return_value={"job_id": "j1", "dataset_id": "d1", "status": "completed"},
            ),
        ):
            assert _run_job_status(args) == 0
        assert "Job ID: j1" in capsys.readouterr().out

    def test_search_without_query_returns_error(self, tmp_path):
        args = SimpleNamespace(
            server="http://x",
            sequence=None,
            file=None,
            dataset=None,
            metadata=False,
            wait=False,
            output=None,
            format="text",
        )
        with patch("malva.search.check_server_connection", return_value=True):
            assert _run_search_data(args) == 1

    def test_search_submits_a_job(self):
        args = SimpleNamespace(
            server="http://x",
            sequence="ACGT",
            file=None,
            dataset="d1",
            metadata=False,
            wait=False,
            output=None,
            format="text",
        )
        with (
            patch("malva.search.check_server_connection", return_value=True),
            patch("malva.search.submit_search_job", return_value={"job_id": "j1"}),
        ):
            assert _run_search_data(args) == 0

    def test_search_unreachable_server(self):
        args = SimpleNamespace(server="http://x", sequence="ACGT", file=None)
        with patch("malva.search.check_server_connection", return_value=False):
            assert _run_search_data(args) == 1

    def test_dispatch_unknown_command(self):
        assert _run_search(SimpleNamespace(command="nope")) == 1

    def test_dispatch_status(self):
        args = SimpleNamespace(
            command="status", server="http://x", job_id="j1", output=None, format="text"
        )
        with (
            patch("malva.search.check_server_connection", return_value=True),
            patch(
                "malva.search.get_job_status",
                return_value={"job_id": "j1", "dataset_id": "d", "status": "completed"},
            ),
        ):
            assert _run_search(args) == 0
