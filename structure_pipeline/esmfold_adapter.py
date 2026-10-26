"""Convert a single-chain ESMFold PDB into ExSquared geometry and viewers.

Run from the repository root::

    python -m structure_pipeline.esmfold_adapter --pdb model.pdb \
        --transcript-id G60740.176.nnc --output-dir output_geometry \
        --gtf expressed_isoforms.gtf

The output layout and field names match the existing ExSquared geometry loader.
ESMFold PDB C-alpha B factors are interpreted as pLDDT. ESMFold does not
provide PAE: pae_self is left empty rather than invented.

Only PDB-to-geometry conversion is included; no GPU inference is run.
"""

from __future__ import annotations

import argparse
import csv
import html
import logging
import math
from dataclasses import dataclass
from pathlib import Path
import re
import sys
from typing import Sequence

LOG = logging.getLogger(__name__)

GEOMETRY_COLUMNS = (
    "transcript_id", "model_index", "chain", "res_seq", "res_name",
    "x", "y", "z", "b_factor", "mean_plddt", "pae_self",
)

_GTF_TRANSCRIPT = re.compile(r'(?:^|;)\s*transcript_id\s+"([^"]+)"')
_COLORS = ((90, "royalblue"), (70, "cornflowerblue"),
           (50, "yellow"), (0, "orange"))


@dataclass(frozen=True)
class Residue:
    chain: str
    res_seq: int
    res_name: str
    x: float
    y: float
    z: float
    plddt: float


@dataclass(frozen=True)
class Exon:
    start: int
    end: int
    strand: str


def encode_transcript_id(transcript_id: str) -> str:
    """Match ExSquared's existing lossy, no-dots/lowercase folder convention."""
    if not transcript_id or not transcript_id.strip():
        raise ValueError("Transcript ID cannot be empty")
    if transcript_id in (".", "..") or any(c in transcript_id for c in "/\\\x00"):
        raise ValueError("Invalid transcript ID for folder naming")
    encoded = transcript_id.replace(".", "").lower()
    if not encoded:
        raise ValueError("Transcript ID becomes empty after encoding")
    return encoded


def read_pdb_calpha(pdb_path: Path, chain: str | None = None) -> list[Residue]:
    """Read C-alpha atoms in the first PDB MODEL; require 1-based contiguous IDs.

    ESMFold generally predicts a single chain with residue indices 1..L.
    Reject gaps, multiple chains and insertion codes rather than silently
    misaligning geometry residues and transcript/CDS coordinates.
    """
    seen: dict[tuple[str, int], Residue] = {}
    chains = set()
    first_model_started = False
    for raw in pdb_path.open("r", encoding="utf-8", errors="replace"):
        if raw.startswith("MODEL "):
            if first_model_started:
                break
            first_model_started = True
            continue
        if raw.startswith("ENDMDL"):
            break
        if not raw.startswith("ATOM  "):
            continue
        if raw[12:16].strip() != "CA" or raw[16:17] not in (" ", "A"):
            continue
        if len(raw) < 66:
            raise ValueError("Malformed PDB C-alpha record (line too short)")
        pdb_chain = raw[21:22].strip() or "A"
        chains.add(pdb_chain)
        if chain is not None and pdb_chain != chain:
            continue
        if raw[26:27].strip():
            raise ValueError("PDB insertion codes are not supported for isoform mapping")
        try:
            num = int(raw[22:26])
            x, y, z = (float(raw[a:b]) for a, b in ((30, 38), (38, 46), (46, 54)))
            conf = float(raw[60:66])
        except ValueError as exc:
            raise ValueError(f"Invalid PDB numeric field: {raw.rstrip()!r}") from exc
        if not (math.isfinite(conf) and 0 <= conf <= 100):
            raise ValueError(f"Invalid pLDDT/B-factor at residue {num}: {conf}")
        if not all(math.isfinite(v) for v in (x, y, z)):
            raise ValueError(f"Nonfinite C-alpha coordinate at residue {num}")
        key = (pdb_chain, num)
        if key in seen:
            raise ValueError(f"Duplicate residue/chain {key} in PDB")
        seen[key] = Residue(pdb_chain, num, raw[17:20].strip(), x, y, z, conf)
    if not seen:
        raise ValueError(f"No C-alpha residues found for chain {chain or '<auto>'} in {pdb_path}")
    selected_chains = {r.chain for r in seen.values()}
    if chain is None and len(selected_chains) != 1:
        raise ValueError(f"Multiple protein chains {sorted(selected_chains)}: use --chain")
    if chain is not None and chain not in selected_chains:
        raise ValueError(f"Chain {chain!r} not found in {pdb_path}")
    residues = sorted(seen.values(), key=lambda r: r.res_seq)
    expected = list(range(1, len(residues) + 1))
    if [r.res_seq for r in residues] != expected:
        raise ValueError("PDB residues must be 1..L without gaps for accurate exon mapping")
    if all(r.plddt == 0 for r in residues):
        raise ValueError("All PDB C-alpha B factors are zero; cannot interpret as ESMFold pLDDT")
    return residues


def read_gtf_exons_and_cds(gtf_path: Path, transcript_id: str) -> tuple[list[Exon], list[Exon]]:
    """Read exact transcript's exon and CDS intervals, in genomic order."""
    exons: set[Exon] = set()
    cds: set[Exon] = set()
    for line in gtf_path.open("r", encoding="utf-8"):
        if not line.strip() or line.startswith("#"):
            continue
        cols = line.rstrip("\n").split("\t")
        if len(cols) < 9 or cols[2] not in ("exon", "CDS"):
            continue
        match = _GTF_TRANSCRIPT.search(cols[8])
        if not match or match.group(1) != transcript_id:
            continue
        if cols[6] not in ("+", "-"):
            raise ValueError(f"Invalid strand for {transcript_id}")
        pos = Exon(int(cols[3]), int(cols[4]), cols[6])
        if pos.start > pos.end:
            raise ValueError("Invalid GTF coordinates")
        (exons if cols[2] == "exon" else cds).add(pos)
    if not exons:
        raise ValueError(f"Transcript {transcript_id!r} has no exon rows in {gtf_path}")
    if not cds:
        raise ValueError(f"Transcript {transcript_id!r} has no CDS rows in {gtf_path}")
    if len({x.strand for x in exons | cds}) != 1:
        raise ValueError("Inconsistent strands in GTF exon/CDS rows")
    exon_list = sorted(exons, key=lambda x: (x.start, x.end))
    cds_list = sorted(cds, key=lambda x: (x.start, x.end))
    for x, y in zip(exon_list, exon_list[1:]):
        if x.end >= y.start:
            raise ValueError("Overlapping exon intervals; cannot assign exon indices unambiguously")
    for x, y in zip(cds_list, cds_list[1:]):
        if x.end >= y.start:
            raise ValueError("Overlapping CDS intervals; cannot assign codons unambiguously")
    for c in cds_list:
        if not any(e.start <= c.start and c.end <= e.end for e in exon_list):
            raise ValueError("CDS interval not wholly contained in an exon")
    return exon_list, cds_list


def exon_residue_indices(exons: Sequence[Exon], cds: Sequence[Exon],
                         protein_length: int) -> dict[int, list[int]]:
    """Map genomic-order 0-based exon indices to residue indices.

    Includes a residue in both exons when its codon crosses an exon junction.
    GTF base positions are inclusive; negative-strand CDS order is reversed.
    One terminal stop codon (not in the predicted protein) is permitted.
    """
    if not exons or not cds:
        return {}
    strand = exons[0].strand
    order = range(len(exons)) if strand == "+" else range(len(exons) - 1, -1, -1)
    cumulative = 0
    mapping: dict[int, list[int]] = {}
    for idx in order:
        exon = exons[idx]
        subsegments = [c for c in cds if exon.start <= c.start and c.end <= exon.end]
        subsegments.sort(key=lambda c: c.start, reverse=(strand == "-"))
        matched: set[int] = set()
        for region in subsegments:
            bases = region.end - region.start + 1
            first = cumulative // 3 + 1
            last = (cumulative + bases - 1) // 3 + 1
            matched.update(range(first, last + 1))
            cumulative += bases
        if matched:
            mapping[idx] = sorted(x for x in matched if x <= protein_length)
    if abs(cumulative - protein_length * 3) > 3:
        raise ValueError(
            f"CDS length ({cumulative} nt) does not match PDB ({protein_length} aa). "
            "Check transcript identity/translation; no exon viewers generated."
        )
    return mapping


def _color(plddt: float) -> str:
    for lower, color in _COLORS:
        if plddt >= lower:
            return color
    return "orange"


def _render_viewer(pdb_text: str, residues: Sequence[Residue], title: str,
                   output: Path, highlight: set[int] | None = None) -> None:
    try:
        import py3Dmol
    except ImportError as exc:
        raise RuntimeError(
            "py3Dmol is needed to produce HTML viewers; install it in the active "
            "Conda environment, or use --csv-only"
        ) from exc
    view = py3Dmol.view(width=700, height=600)
    view.addModel(pdb_text, "pdb")
    if highlight is not None:
        view.setStyle({"cartoon": {"color": "lightgrey", "opacity": 0.25}})
    else:
        view.setStyle({"cartoon": {"color": "white"}})
    for residue in residues:
        if highlight is not None and residue.res_seq not in highlight:
            continue
        view.setStyle(
            {"chain": residue.chain, "resi": residue.res_seq},
            {"cartoon": {"color": _color(residue.plddt), "opacity": 1.0}},
        )
    view.zoomTo()
    view.spin(False)
    body = view._make_html()
    output.write_text(
        "<!DOCTYPE html>\n<html><head><meta charset='utf-8'><title>"
        + html.escape(title) + "</title></head>"
        + "<body style='margin:0;padding:10px;background:#fff;'>"
        + "<h3>" + html.escape(title) + "</h3>" + body
        + "</body></html>", encoding="utf-8"
    )


def convert(pdb: Path, transcript_id: str, output_dir: Path, *,
            gtf: Path | None = None, model_index: int = 0,
            chain: str | None = None, csv_only: bool = False,
            overwrite: bool = False) -> Path:
    """Convert one PDB, returning the transcript output directory."""
    pdb = Path(pdb)
    if not pdb.is_file():
        raise FileNotFoundError(pdb)
    if gtf is not None and not Path(gtf).is_file():
        raise FileNotFoundError(gtf)
    if model_index < 0:
        raise ValueError("model_index must be >= 0")
    folder_name = encode_transcript_id(transcript_id)
    residues = read_pdb_calpha(pdb, chain=chain)
    mapping: dict[int, list[int]] = {}
    if gtf is not None:
        exons, cds = read_gtf_exons_and_cds(Path(gtf), transcript_id)
        mapping = exon_residue_indices(exons, cds, len(residues))
    if not csv_only:
        try:
            import py3Dmol  # noqa: F401
        except ImportError as exc:
            raise RuntimeError("Install py3Dmol or supply --csv-only") from exc
    out = Path(output_dir) / folder_name
    csv_file = out / f"geometry_model{model_index}.csv"
    html_file = out / f"viewer_model{model_index}.html"
    viewer_files = [out / f"viewer_model{model_index}_exon{i}.html" for i in mapping if mapping[i]]
    targets = [csv_file] + ([] if csv_only else [html_file, *viewer_files])
    if not overwrite and any(p.exists() for p in targets):
        raise FileExistsError(
            f"Existing output in {out}. Use a new --output-dir or --overwrite."
        )
    out.mkdir(parents=True, exist_ok=True)
    with csv_file.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=GEOMETRY_COLUMNS)
        writer.writeheader()
        for r in residues:
            writer.writerow(dict(transcript_id=transcript_id, model_index=model_index,
                                 chain=r.chain, res_seq=r.res_seq, res_name=r.res_name,
                                 x=f"{r.x:.4f}", y=f"{r.y:.4f}", z=f"{r.z:.4f}",
                                 b_factor=f"{r.plddt:.4f}", mean_plddt=f"{r.plddt:.4f}",
                                 pae_self=""))
    LOG.info("Wrote %d C-alpha rows: %s", len(residues), csv_file)
    if not csv_only:
        pdb_text = pdb.read_text(encoding="utf-8", errors="replace")
        _render_viewer(pdb_text, residues, transcript_id, html_file)
        for exon_idx, matched in mapping.items():
            if matched:
                _render_viewer(pdb_text, residues, f"{transcript_id} – exon {exon_idx + 1}",
                               out / f"viewer_model{model_index}_exon{exon_idx}.html",
                               set(matched))
        LOG.info("Wrote %d exon HTML viewers", len(viewer_files))
    return out


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pdb", type=Path, required=True)
    parser.add_argument("--transcript-id", required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--gtf", type=Path, help="Exact transcript GTF for exon HTML viewers")
    parser.add_argument("--model-index", type=int, default=0)
    parser.add_argument("--chain", help="Select a PDB chain (normally A for ESMFold)")
    parser.add_argument("--csv-only", action="store_true", help="Skip HTML viewers")
    parser.add_argument("--overwrite", action="store_true", help="Replace existing outputs")
    args = parser.parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(levelname)s | %(message)s")
    try:
        out = convert(args.pdb, args.transcript_id, args.output_dir,
                      gtf=args.gtf, model_index=args.model_index,
                      chain=args.chain, csv_only=args.csv_only,
                      overwrite=args.overwrite)
    except (OSError, ValueError, RuntimeError) as exc:
        LOG.error("Conversion failed: %s", exc)
        return 1
    print(f"ExSquared geometry directory: {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
