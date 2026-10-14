import pytest

import malva.serve.serve as serve


@pytest.fixture
def app_client(monkeypatch, tmp_path):
    # flask_session writes a ./flask_session directory; keep it in the temp dir.
    monkeypatch.chdir(tmp_path)
    had_state = hasattr(serve, "global_state")
    original = getattr(serve, "global_state", None)
    serve.global_state = serve.GlobalState()
    try:
        yield serve.create_app(init_state=False).test_client()
    finally:
        if had_state:
            setattr(serve, "global_state", original)
        else:
            delattr(serve, "global_state")


def test_health_reports_initializing(app_client):
    response = app_client.get("/health")
    assert response.status_code == 200
    assert response.get_json()["status"] == "initializing"


def test_parse_queried_without_session_data(app_client):
    response = app_client.post("/parse_queried")
    assert response.status_code == 200
    assert response.get_json()["success"] is False
