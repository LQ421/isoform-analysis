"""Fast, CPU-only tests for 3D exon HTML selection."""
from isoform_dashboard.viewer_selection import select_viewer_file, lookup_viewer_variant


def test_select_existing_exon(tmp_path):
    full = tmp_path / "viewer_model0.html"
    exon = tmp_path / "viewer_model0_exon5.html"
    full.write_text("full", encoding="utf-8")
    exon.write_text("exon", encoding="utf-8")
    geo = {"html": str(full), "exon_htmls": [None] * 5 + [str(exon)]}
    path, idx = select_viewer_file(geo, "TX.1", {"transcript_id": "TX.1", "exon_idx": 5})
    assert path == str(exon) and idx == 5
    assert lookup_viewer_variant(geo, "exon5") == str(exon)


def test_full_and_wrong_transcript(tmp_path):
    full = tmp_path / "viewer_model0.html"
    full.write_text("full", encoding="utf-8")
    geo = {"html": str(full), "exon_htmls": []}
    assert select_viewer_file(geo, "TX.1", None) == (str(full), None)
    assert select_viewer_file(geo, "TX.1", {"transcript_id": "TX.2", "exon_idx": 5}) == (str(full), None)
    assert lookup_viewer_variant(geo, "full") == str(full)


def test_missing_exon_falls_back(tmp_path):
    full = tmp_path / "viewer_model0.html"
    full.write_text("full", encoding="utf-8")
    geo = {"html": str(full), "exon_htmls": [None, str(tmp_path / "missing.html")]}
    assert select_viewer_file(geo, "TX", {"transcript_id": "TX", "exon_idx": 1}) == (str(full), None)
    assert lookup_viewer_variant(geo, "exon1") is None


def test_invalid_indices_and_paths(tmp_path):
    full = tmp_path / "viewer_model0.html"
    full.write_text("full", encoding="utf-8")
    geo = {"html": str(full), "exon_htmls": [None]}
    assert select_viewer_file(geo, "TX", {"transcript_id": "TX", "exon_idx": -1}) == (str(full), None)
    assert lookup_viewer_variant(geo, "../../file") is None
    assert lookup_viewer_variant(geo, "exon999") is None
    assert lookup_viewer_variant(geo, "exon-1") is None
