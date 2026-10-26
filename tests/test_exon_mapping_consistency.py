"""Regression tests: ExSquared legacy mapping and ESMFold adapter agree."""
import pytest
from utilities.extract_3d_geometry import compute_exon_residue_range
from structure_pipeline.esmfold_adapter import Exon, exon_residue_indices


def fixtures(lengths, strand="+"):
    exons, cds, legacy = [], [], []
    for i, nt in enumerate(lengths):
        start = 1000 + 1000 * i
        end = start + max(nt, 1) - 1
        exons.append(Exon(start, start + 900, strand))
        if nt:
            cds.append(Exon(start, end, strand))
        legacy.append(dict(exon_start=start, exon_end=start + 900,
                           cds_start=start if nt else None,
                           cds_end=end if nt else None, strand=strand))
    return exons, cds, legacy


@pytest.mark.parametrize("strand", ["+", "-"])
def test_coding_exons_match_adapter(strand):
    exons, cds, legacy = fixtures([0,120,117,144,89,300,75,16], strand)
    mapping = exon_residue_indices(exons, cds, 286)
    for i in range(len(exons)):
        expected = mapping.get(i)
        span = compute_exon_residue_range(i, legacy, protein_length=286)
        assert span == ((min(expected),max(expected)) if expected else None)


def test_one_base_coding_exon_and_shared_junction():
    exons, cds, legacy = fixtures([2,1,3])
    mapping = exon_residue_indices(exons, cds, 2)
    assert mapping == {0:[1], 1:[1], 2:[2]}
    assert compute_exon_residue_range(1, legacy, protein_length=2) == (1,1)


def test_stop_only_exon_is_not_highlighted():
    _, _, legacy = fixtures([6,3])
    assert compute_exon_residue_range(0, legacy, protein_length=2) == (1,2)
    assert compute_exon_residue_range(1, legacy, protein_length=2) is None


def test_noncoding_exon():
    _, _, legacy = fixtures([0,6,0,3])
    assert compute_exon_residue_range(0, legacy, protein_length=2) is None
    assert compute_exon_residue_range(2, legacy, protein_length=2) is None
