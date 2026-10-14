# Publishing to PyPI

Malva ships a source distribution plus Linux (manylinux) and macOS wheels. The
sdist builds the Cython extensions from source; the wheels are built with
[cibuildwheel](https://cibuildwheel.readthedocs.io) — the same tool used by
`build_scripts/build_with_docker.sh` — so the Linux wheels carry the manylinux
tags PyPI accepts.

## Releasing from CI (recommended)

`.github/workflows/release-pypi.yml` builds the sdist and the wheels and uploads
everything to PyPI. It publishes whatever `version` is in `pyproject.toml`, so
bump the version (and regenerate `setup.py`, see `build_scripts/README.md`)
before releasing.

Trigger it either way:

- GitHub → **Actions** → **release-pypi** → **Run workflow** (leave *Publish to
  PyPI* checked), or
- push a tag: `git tag v1.1.0 && git push origin v1.1.0`.

Run it from a green `main`. PyPI rejects a version that already exists, so a
release is one-shot per version bump.

### One-time authentication setup

Pick one of the two options below.

**Trusted publishing (recommended — nothing to rotate).** On
<https://pypi.org/manage/account/publishing/> add a *pending publisher*:

| Field | Value |
|---|---|
| PyPI project name | `malva` |
| Owner | `malva-bio` |
| Repository name | `malva` |
| Workflow name | `release-pypi.yml` |
| Environment name | `pypi` |

The workflow already requests `id-token: write` and targets the `pypi`
environment, so that is all. (Add the same pending publisher on TestPyPI if you
want a dry run.)

**API token.** Create a project token at
<https://pypi.org/manage/account/token/>, add it as the repository secret
`PYPI_API_TOKEN` (Settings → Secrets and variables → Actions → New repository
secret), then in the `publish` job swap the trusted-publishing step for:

```yaml
    steps:
      - uses: actions/download-artifact@v4
        with:
          path: dist
          merge-multiple: true
      - uses: pypa/gh-action-pypi-publish@release/v1
        with:
          password: ${{ secrets.PYPI_API_TOKEN }}
```

## Publishing locally

Keep this for a dry run or an offline release:

```bash
pip install build twine

cd /tmp
python -m build /local/malva/repos/malva --outdir /local/malva/repos/malva/dist
python -m twine upload /local/malva/repos/malva/dist/*
```

Or build the platform wheels with the scripts in `build_scripts/` (see its
README) and upload them with `twine upload wheelhouse/*`.
