#!/usr/bin/env python3
"""Pre-publish guard for the PyPI release workflow.

Refuse to publish unless the pushed tag matches the version recorded in
``pyproject.toml``. Tag format: ``v<semver>-py`` (e.g. ``v0.2.1-py``).
The ``-py-client`` tag is handled by a different workflow and rejected
here.
"""

from __future__ import annotations

import os
import re
import sys
from pathlib import Path

import tomllib

ROOT = Path(__file__).resolve().parent.parent
PYPROJECT = ROOT / "pyproject.toml"

TAG_RE = re.compile(r"^v(.+)-py$")


def main() -> int:
    tag = os.environ.get("GITHUB_REF_NAME")
    if not tag:
        print(
            "verify-pypi-tag: GITHUB_REF_NAME is not set; refusing to publish.",
            file=sys.stderr,
        )
        return 1

    if tag.endswith("-py-client"):
        print(
            f"verify-pypi-tag: tag {tag!r} targets the client package, not iscn-authenticator.",
            file=sys.stderr,
        )
        return 1

    match = TAG_RE.match(tag)
    if not match:
        print(
            f"verify-pypi-tag: tag {tag!r} does not match expected v<semver>-py pattern.",
            file=sys.stderr,
        )
        return 1

    tag_version = match.group(1)

    with PYPROJECT.open("rb") as fh:
        data = tomllib.load(fh)
    pkg_version = data["project"]["version"]

    if tag_version != pkg_version:
        print(
            f"verify-pypi-tag: tag version ({tag_version}) does not match " f"pyproject.toml ({pkg_version}).",
            file=sys.stderr,
        )
        return 1

    print(f"verify-pypi-tag: ok ({tag} -> {pkg_version})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
