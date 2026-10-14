#!/usr/bin/env python3
# Copyright (c) 2025 Daniel León-Periñán and Nikolaos Karaiskos
#                    Rajewsky Lab, Max Delbrück Center for Molecular Medicine (MDC), Berlin
#
# Non-commercial and academic use only. See LICENSE for full terms.
"""Regenerate the committed ``setup.py`` from ``pyproject.toml``.

``pyproject.toml`` is the single source of truth for the package metadata
(name, version, dependencies, Python range, entry points). Poetry's
``[tool.poetry.build] generate-setup-file = true`` emits a ``setup.py`` when it
builds an sdist, but that copy lives inside the tarball, and Poetry's
``setup_py()`` deliberately reuses an existing root ``setup.py`` rather than
rewriting it.

The repository keeps ``setup.py`` at its root so the setuptools paths keep
working (``setup.py build_ext --inplace`` and ``setup.py bdist_wheel`` in
``build_scripts/build_with_apptainer.sh``). This script asks Poetry for the
generated content, re-adds the licence header Poetry's template omits, and
writes it to ``setup.py`` so the committed file cannot drift from
``pyproject.toml``.

Usage::

    uvx --from poetry python build_scripts/generate_setup.py
    uvx --from poetry python build_scripts/generate_setup.py --check

``--check`` regenerates in memory and fails if ``setup.py`` is out of date, so
CI can enforce the invariant once the pipeline lands.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

HEADER = (
    "# Copyright (c) 2025 Daniel León-Periñán and Nikolaos Karaiskos\n"
    "#                    Rajewsky Lab, Max Delbrück Center for Molecular Medicine (MDC), Berlin\n"
    "#\n"
    "# Non-commercial and academic use only. See LICENSE for full terms.\n"
)
SETUP_PATH = Path(__file__).resolve().parent.parent / "setup.py"
PROJECT_ROOT = SETUP_PATH.parent


def render() -> str:
    """Return the setup.py Poetry would generate, with the licence header."""
    from poetry.core.factory import Factory
    from poetry.core.masonry.builders.sdist import SdistBuilder

    poetry = Factory().create_poetry(str(PROJECT_ROOT))
    generated = SdistBuilder(poetry).build_setup().decode()
    first_line, _, rest = generated.partition("\n")
    return f"{first_line}\n{HEADER}\n{rest}"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--check",
        action="store_true",
        help="exit non-zero if setup.py is not up to date instead of writing it",
    )
    args = parser.parse_args()

    rendered = render()
    current = SETUP_PATH.read_text(encoding="utf-8") if SETUP_PATH.exists() else None

    if args.check:
        if current == rendered:
            print("setup.py is up to date with pyproject.toml")
            return 0
        print(
            "setup.py is out of date; run build_scripts/generate_setup.py",
            file=sys.stderr,
        )
        return 1

    if current == rendered:
        print("setup.py already up to date")
        return 0

    SETUP_PATH.write_text(rendered, encoding="utf-8")
    print(f"wrote {SETUP_PATH.relative_to(PROJECT_ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
