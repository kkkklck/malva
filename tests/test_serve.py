import pytest


@pytest.fixture
def app_client(monkeypatch, tmp_path):
    # flask_session resolves its storage dir at import time from the cwd, so
    # import the app *after* chdir to keep ./flask_session inside the temp dir
    # (and out of the repository).
    monkeypatch.chdir(tmp_path)
    import malva.serve.serve as serve

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
