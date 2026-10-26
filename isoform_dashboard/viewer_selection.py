"""Select pre-rendered ExSquared 3D viewer files without Dash dependencies.

The geometry directory is indexed when the dashboard starts. The returned paths
are always drawn from those indexed entries, never built from user paths.
"""

from __future__ import annotations

from pathlib import Path


def select_viewer_file(geo: dict, transcript_id: str, selected_exon: dict | None):
    """Return (path, exon_index_or_None), falling back to the full viewer.

    Exon indices are zero-based *genomic* indices and must match the sparse
    exon_htmls list created by discover_exon_viewers().
    """
    full = geo.get("html")
    if selected_exon and selected_exon.get("transcript_id") == transcript_id:
        idx = selected_exon.get("exon_idx")
        exons = geo.get("exon_htmls") or []
        if isinstance(idx, int) and not isinstance(idx, bool) and 0 <= idx < len(exons):
            path = exons[idx]
            if path and Path(path).is_file():
                return str(path), idx
    return (str(full) if full else None), None


def lookup_viewer_variant(geo: dict, variant: str):
    """Resolve a requested viewer variant from an already-indexed geometry entry.

    User-supplied names never become filesystem paths. Missing files return None.
    """
    if variant == "full":
        path = geo.get("html")
    elif variant.startswith("exon") and variant[4:].isdecimal():
        idx = int(variant[4:])
        exons = geo.get("exon_htmls") or []
        path = exons[idx] if idx < len(exons) else None
    else:
        return None
    return str(path) if path and Path(path).is_file() else None
