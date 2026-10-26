# ESMFold PDB to ExSquared: structure geometry and exon viewers

This optional module converts **precomputed ESMFold PDB files** into the geometry CSV and interactive HTML layout consumed by ExSquared. It does **not** perform protein structure inference. High-throughput Snakemake/SLURM orchestration is outside the scope of this integration.

## Requirements

- Python 3.11+ and the packages in the repository's existing environment.
- `py3Dmol` to generate interactive HTML viewers. Without it, use `--csv-only` to export CSV geometry.
- Optional, for exon-specific viewers: a GTF file containing **exon** and **CDS** records for the exact transcript ID.
- PDB with a contiguous 1-based single-chain C-alpha trace and confidence values in the C-alpha B-factor field, as exported by ESMFold. Other PDBs with ordinary experimental B factors must not be interpreted as pLDDT.

## Convert one predicted structure

Run from the project root (substitute your own approved input files and destinations):

```bash
python -m structure_pipeline.esmfold_adapter \
  --pdb path/to/predicted.pdb \
  --transcript-id TRANSCRIPT_ID \
  --gtf path/to/transcripts.gtf \
  --output-dir path/to/output_geometry
```

Omit `--gtf` for a full-protein viewer only. Add `--csv-only` to skip HTML generation. Use `--chain A` if selection of a PDB chain is necessary. Existing outputs are protected by default; use a separate output directory or `--overwrite` deliberately.

Expected layout:

```text
output_geometry/
  <transcript-id-lowercase-with-dots-removed>/
    geometry_model0.csv
    viewer_model0.html
    viewer_model0_exonN.html  # only for coding exons when --gtf is supplied
```

The geometry CSV has these columns:

```text
transcript_id,model_index,chain,res_seq,res_name,x,y,z,b_factor,mean_plddt,pae_self
```

The output keeps `pae_self` empty because ESMFold PDB outputs do not provide PAE. `mean_plddt` is taken from the PDB C-alpha B-factor values **only for compatible ESMFold outputs**.

## Launch ExSquared with the converted results

```bash
python -m isoform_dashboard.dashboard_app \
  --exons path/to/transcripts.gtf \
  --proteins path/to/proteins.fasta \
  --geometry-dir path/to/output_geometry \
  --port 8052
```

Within ExSquared select a gene and a transcript that has geometry, then click a coding exon. The full protein should appear pale grey except for residues highlighted using the pLDDT colour scheme. Clicking **Full structure** restores the full-colour model.

The viewer endpoint serves **only already-indexed HTML files**, using a model variant such as `full` or `exon5`; it must not use a URL parameter as an arbitrary filesystem path. Viewer files are rendered as HTML and JavaScript, so only load trusted prediction outputs.

## Quality control and known constraints

Before interpreting isoform changes:

1. Compare the PDB-derived amino-acid sequence with the exact transcript's translated protein FASTA. Require matching sequence, not merely equal length.
2. Verify PDB residue numbering is contiguous and begins at 1.
3. Validate total CDS length against translated protein length, allowing for a terminal stop codon where appropriate.
4. Check exon-junction codons, GTF CDS phase, strand and transcript identity. A residue may overlap two adjacent exons when its codon straddles an exon junction. Genomic-coordinate and residue-boundary agreement alone is **not** complete CDS-phase validation.
5. Confirm the HTML exon viewer highlights the same residues as the protein sequence panel.

The inherited ExSquared folder-key convention removes dots and lowercases transcript IDs. It is **not collision-free**: do not process distinct IDs that collapse to the same key in one output directory. Outputs are limited to one selected PDB chain with no missing C-alpha residues or insertion codes.

No model weights, unpublished input sequences, raw CELLGEN project files, private logs, HPC account identifiers, or real PDBs are included in this feature. Provide synthetic fixtures for tests. Check project permissions before sharing derived output files.

## Tests

Use the appropriate environment in which these tests' dependencies are installed:

```bash
python -m pytest tests/test_esmfold_adapter.py -q
python -m pytest tests/test_viewer_selection.py tests/test_exon_mapping_consistency.py -q
python -m py_compile isoform_dashboard/app_layout.py utilities/extract_3d_geometry.py structure_pipeline/esmfold_adapter.py
```

A passing unit-test suite does not substitute for a browser-based smoke test of full/exon switching or validation against the reference GTF and FASTA.
