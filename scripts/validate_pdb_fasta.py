r"""Check exact amino-acid agreement between a transcript FASTA and ESMFold PDBs.

Example (Windows CMD, from the ExSquared repository):
    python scripts\validate_pdb_fasta.py --fasta "path\to\proteins.fa" --pdb-dir "path\to\pdb_files" --transcript-id "TRANSCRIPT_ID_1" --transcript-id "TRANSCRIPT_ID_2"

Only reads user data. It never modifies the FASTA or PDB files.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

from Bio import SeqIO
from Bio.PDB import PDBParser
from Bio.SeqUtils import seq1


def check_one(transcript_id: str, fasta_records: dict, pdb_dir: Path) -> bool:
    pdb_path = pdb_dir / f"{transcript_id}.pdb"
    print(f"\nTRANSCRIPT: {transcript_id}")

    if transcript_id not in fasta_records:
        print("FAIL: transcript ID not found in FASTA (exact ID match required)")
        return False
    if not pdb_path.is_file():
        print(f"FAIL: PDB file not found: {pdb_path}")
        return False

    ref_raw = str(fasta_records[transcript_id].seq).upper()
    ref = ref_raw.removesuffix("*")
    if "*" in ref:
        print("FAIL: internal stop codon in FASTA")
        return False

    structure = PDBParser(QUIET=True).get_structure(transcript_id, str(pdb_path))
    model = next(structure.get_models(), None)
    if model is None:
        print("FAIL: PDB contains no model")
        return False

    chains = [chain for chain in model if any("CA" in residue for residue in chain)]
    if not chains:
        print("FAIL: PDB has no C-alpha atoms")
        return False
    if len(chains) != 1:
        print(f"FAIL: expected one protein chain; found {[c.id for c in chains]}")
        return False

    chain = chains[0]
    residues = [residue for residue in chain if "CA" in residue]
    pdb_seq = "".join(seq1(res.get_resname(), custom_map={"MSE": "M"}, undef_code="X") for res in residues)
    numbers = [res.id[1] for res in residues]
    expected_numbers = list(range(1, len(residues) + 1))

    print(f"FASTA: {len(ref)} aa | PDB: {len(pdb_seq)} aa | chain: {chain.id}")
    print(f"PDB residue numbering: {numbers[0]}..{numbers[-1]}")
    seq_match = ref == pdb_seq
    numbered = numbers == expected_numbers
    if seq_match:
        print("Sequence identity: 100.00% (exact match)")
    else:
        mismatches = [(i, a, b) for i, (a, b) in enumerate(zip(ref, pdb_seq), start=1) if a != b]
        matches = sum(a == b for a, b in zip(ref, pdb_seq))
        print(f"Sequence identity over corresponding positions: {matches}/{max(len(ref),len(pdb_seq))}")
        print("First mismatches (position, FASTA, PDB):", mismatches[:12])
        if len(ref) != len(pdb_seq):
            print("WARNING: lengths differ; corresponding-position comparison is not an alignment")
    if not numbered:
        print("FAIL: PDB C-alpha residue numbers are not consecutive 1..N")

    passed = seq_match and numbered
    print("RESULT:", "PASS" if passed else "FAIL")
    return passed


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fasta", type=Path, required=True)
    parser.add_argument("--pdb-dir", type=Path, required=True)
    parser.add_argument("--transcript-id", action="append", required=True, help="Repeat for each transcript ID")
    args = parser.parse_args(argv)
    if not args.fasta.is_file():
        parser.error(f"FASTA not found: {args.fasta}")
    if not args.pdb_dir.is_dir():
        parser.error(f"PDB directory not found: {args.pdb_dir}")
    records = SeqIO.to_dict(SeqIO.parse(str(args.fasta), "fasta"))
    results = [check_one(t, records, args.pdb_dir) for t in args.transcript_id]
    print(f"\nSUMMARY: {sum(results)}/{len(results)} passed")
    return 0 if all(results) else 1


if __name__ == "__main__":
    sys.exit(main())
