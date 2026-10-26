"""CPU-only tests: run `python -m pytest tests/test_esmfold_adapter.py -q`."""
import csv
from pathlib import Path

import pytest

from structure_pipeline.esmfold_adapter import (
    GEOMETRY_COLUMNS,
    convert,
    encode_transcript_id,
    exon_residue_indices,
    read_gtf_exons_and_cds,
    read_pdb_calpha,
)
from isoform_dashboard.alphafold_geometry import (
    build_alphafold_geometry_mapping,
    resolve_alphafold_geometry,
)


def _make_pdb(path: Path, count=6):
    with path.open("w") as fh:
        for i in range(1, count + 1):
            fh.write(
                f"ATOM   {i:4d}  CA  ALA A{i:4d}    {i*1.1:8.3f}{i*2.2:8.3f}{i*3.3:8.3f}"
                f"  1.00{i*10.0:6.2f}           C\n"
            )
        fh.write("TER\nEND\n")


def _make_gtf(path: Path, strand="+"):
    # 5 coding bases in exon 0 and 13 in exon 1 -> 6 residues.
    # Codon 2 is encoded by both exons; both should highlight residue 2.
    with path.open("w") as fh:
        for feature, s, e in [("exon", 100, 104), ("CDS", 100, 104),
                              ("exon", 200, 212), ("CDS", 200, 212)]:
            fh.write(f'chr1\tunit\t{feature}\t{s}\t{e}\t.\t{strand}\t0\tgene_id "g"; transcript_id "T.1";\n')


def test_encoding():
    assert encode_transcript_id("G60740.176.nnc") == "g60740176nnc"
    assert encode_transcript_id("G60740.186_ENST00000361752.8") == "g60740186_enst000003617528"
    with pytest.raises(ValueError):
        encode_transcript_id("../escape")


def test_pdb_csv_loader(tmp_path):
    pdb = tmp_path / "six.pdb"
    _make_pdb(pdb)
    out = convert(pdb, "T.1", tmp_path / "results", csv_only=True)
    csv_path = out / "geometry_model0.csv"
    with csv_path.open(newline="") as fh:
        reader = csv.DictReader(fh)
        rows = list(reader)
    assert reader.fieldnames == list(GEOMETRY_COLUMNS)
    assert len(rows) == 6
    assert rows[0]["res_seq"] == "1"
    assert rows[0]["mean_plddt"] == "10.0000"
    assert all(row["pae_self"] == "" for row in rows)
    geometry_mapping = build_alphafold_geometry_mapping(str(tmp_path / "results"))
    resolved = resolve_alphafold_geometry("T.1", geometry_mapping)
    assert resolved is not None
    assert Path(resolved["csv"]) == csv_path
    assert resolved["html"] is None  # CSV-only does not generate a full viewer


def test_gtf_cds_plus_strand_codon_boundary(tmp_path):
    gtf = tmp_path / "mini.gtf"
    _make_gtf(gtf, "+")
    exons, cds = read_gtf_exons_and_cds(gtf, "T.1")
    mapping = exon_residue_indices(exons, cds, protein_length=6)
    assert mapping[0] == [1, 2]
    assert mapping[1] == [2, 3, 4, 5, 6]


def test_gtf_cds_minus_strand_codon_boundary(tmp_path):
    gtf = tmp_path / "mini.gtf"
    _make_gtf(gtf, "-")
    exons, cds = read_gtf_exons_and_cds(gtf, "T.1")
    mapping = exon_residue_indices(exons, cds, protein_length=6)
    assert mapping[1] == [1, 2, 3, 4, 5]
    assert mapping[0] == [5, 6]


def test_reject_incorrect_transcript_or_length(tmp_path):
    pdb = tmp_path / "six.pdb"
    gtf = tmp_path / "mini.gtf"
    _make_pdb(pdb, count=6)
    _make_gtf(gtf)
    with pytest.raises(ValueError, match="no exon"):
        convert(pdb, "Other", tmp_path / "out", gtf=gtf, csv_only=True)
    with pytest.raises(ValueError, match="CDS length"):
        exon_residue_indices(*read_gtf_exons_and_cds(gtf, "T.1"), protein_length=10)


def test_protect_existing_outputs(tmp_path):
    pdb = tmp_path / "six.pdb"
    _make_pdb(pdb)
    convert(pdb, "T.1", tmp_path / "results", csv_only=True)
    with pytest.raises(FileExistsError):
        convert(pdb, "T.1", tmp_path / "results", csv_only=True)


def test_no_bfactor_confidence(tmp_path):
    pdb = tmp_path / "none.pdb"
    _make_pdb(pdb)
    pdb.write_text(pdb.read_text().replace(" 10.00", "  0.00").replace(" 20.00", "  0.00")
                   .replace(" 30.00", "  0.00").replace(" 40.00", "  0.00")
                   .replace(" 50.00", "  0.00").replace(" 60.00", "  0.00"))
    with pytest.raises(ValueError, match="zero"):
        read_pdb_calpha(pdb)


def test_html_viewer_names_and_exon_discovery(tmp_path, monkeypatch):
    """Contract test for HTML file discovery; uses a CPU-only py3Dmol stub."""
    import sys
    from types import SimpleNamespace
    from isoform_dashboard.alphafold_geometry import discover_exon_viewers

    class FakeView:
        def addModel(self, *_):
            return self
        def setStyle(self, *_):
            return self
        def zoomTo(self):
            return self
        def spin(self, *_):
            return self
        def _make_html(self):
            return "<div>synthetic viewer placeholder</div>"

    monkeypatch.setitem(sys.modules, "py3Dmol", SimpleNamespace(view=lambda **_: FakeView()))
    pdb = tmp_path / "six.pdb"
    gtf = tmp_path / "mini.gtf"
    _make_pdb(pdb)
    _make_gtf(gtf)
    out = convert(pdb, "T.1", tmp_path / "results", gtf=gtf)
    assert (out / "viewer_model0.html").is_file()
    assert (out / "viewer_model0_exon0.html").is_file()
    assert (out / "viewer_model0_exon1.html").is_file()
    m = build_alphafold_geometry_mapping(str(tmp_path / "results"))
    assert discover_exon_viewers(m) == 2
    assert len(m["t1"]["exon_htmls"]) == 2
    assert m["t1"]["html"] == str(out / "viewer_model0.html")
