from __future__ import annotations

import json
from pathlib import Path


def test_release_versions_are_consistent() -> None:
    repository = Path(__file__).resolve().parents[1]
    version = (repository / "VERSION").read_text(encoding="utf-8").strip()
    build = json.loads((repository / "BUILD_INFO.json").read_text(encoding="utf-8"))
    zenodo = json.loads((repository / ".zenodo.json").read_text(encoding="utf-8"))
    citation = (repository / "CITATION.cff").read_text(encoding="utf-8")
    assert version == "2.0.1"
    assert build["software_release"] == version
    assert build["release_tag"] == f"v{version}"
    assert zenodo["version"] == version
    assert f"version: {version}" in citation


def test_distribution_metadata_is_complete() -> None:
    repository = Path(__file__).resolve().parents[1]
    zenodo = json.loads((repository / ".zenodo.json").read_text(encoding="utf-8"))
    assert zenodo["upload_type"] == "software"
    assert zenodo["access_right"] == "open"
    assert zenodo["license"] == "mit"
    assert zenodo["creators"][0]["orcid"] == "0000-0003-4447-4104"
    assert len(zenodo["description"]) > 100
