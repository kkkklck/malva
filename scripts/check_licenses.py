"""Fail a release when an installed dependency is not open source.

Reads a ``pip-licenses --format=json`` report (from a path or stdin) and
classifies every package's license. By default only a verifiably closed-source
license fails; a license that cannot be classified is reported but allowed, so an
unlabelled package does not break the build. Pass ``--strict`` to fail on those
too.
"""

import argparse
import json
import sys

OPEN_SOURCE_TOKENS = (
    "mit",
    "expat",
    "bsd",
    "apache",
    "gpl",
    "lgpl",
    "agpl",
    "mpl",
    "mozilla",
    "isc",
    "psf",
    "python software foundation",
    "python-2.0",
    "zlib",
    "unlicense",
    "cc0",
    "cc-by",
    "cc by",
    "creative commons",
    "wtfpl",
    "boost",
    "artistic",
    "public domain",
    "0bsd",
    "hpnd",
    "blueoak",
    "mulanpsl",
    "eupl",
    "cddl",
    "epl",
    "zpl",
    "zope",
    "open software",
    "osl",
    "cdla",
    "odc",
)

NON_OPEN_SOURCE_TOKENS = (
    "proprietary",
    "commercial",
    "all rights reserved",
    "confidential",
    "not open source",
)

UNKNOWN_ALIASES = {"", "unknown", "none", "unlicensed", "unclear", "n/a", "not specified"}


def classify_license(license_text):
    text = (license_text or "").strip().lower()
    if text in UNKNOWN_ALIASES:
        return "unknown"
    if any(token in text for token in OPEN_SOURCE_TOKENS):
        return "open"
    if any(token in text for token in NON_OPEN_SOURCE_TOKENS):
        return "closed"
    return "unknown"


def read_report(path):
    if path:
        with open(path, encoding="utf-8") as fh:
            return json.load(fh)
    return json.load(sys.stdin)


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Fail when an installed dependency is not recognised as open source."
    )
    parser.add_argument(
        "report",
        nargs="?",
        help="pip-licenses --format=json output file (defaults to stdin)",
    )
    parser.add_argument(
        "--strict",
        action="store_true",
        help="also fail on licenses that cannot be classified, not just proprietary ones",
    )
    parser.add_argument(
        "--allow",
        action="append",
        default=[],
        help="package name to allow despite its license (repeatable)",
    )
    args = parser.parse_args(argv)

    records = read_report(args.report)
    allowed = {name.lower() for name in args.allow}

    not_open = [
        r
        for r in records
        if classify_license(r.get("License", "")) != "open"
        and str(r.get("Name", "")).lower() not in allowed
    ]
    unknown = [r for r in not_open if classify_license(r.get("License", "")) == "unknown"]
    closed = [r for r in not_open if classify_license(r.get("License", "")) == "closed"]

    print(
        f"checked {len(records)} packages: "
        f"{len(records) - len(not_open)} open source, "
        f"{len(unknown)} unknown, {len(closed)} non-open-source"
    )
    for record in closed:
        print(
            f"  NON-OPEN-SOURCE: {record.get('Name')} {record.get('Version', '')} "
            f"-> {record.get('License', '')!r}"
        )
    for record in unknown:
        print(
            f"  UNKNOWN LICENSE: {record.get('Name')} {record.get('Version', '')} "
            f"-> {record.get('License', '')!r}"
        )

    offenders = not_open if args.strict else closed
    if offenders:
        print("license audit FAILED")
        return 1
    print("license audit OK")
    return 0


if __name__ == "__main__":
    sys.exit(main())
