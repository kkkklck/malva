from unittest.mock import patch

from malva.__main__ import run_malva


def test_run_malva_delegates_to_the_cli():
    with patch("malva.__main__.cmdline_main") as cmdline_main:
        run_malva()
    cmdline_main.assert_called_once_with()
