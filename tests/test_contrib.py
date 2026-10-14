import malva.contrib as contrib


def test_module_metadata():
    assert contrib.__version__ == 1.0
    assert contrib.__author__ == ["Daniel León-Periñán"]
    assert contrib.__email__ == ["daniel.leonperinan@mdc-berlin.de"]
    assert "MDC Berlin" in contrib.__license__
