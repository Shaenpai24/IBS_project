# PRNP Mutation Pipeline (with AlphaFold Structural Features)

This workspace builds PRNP single-point mutation data and engineered features in staged steps.

## Current Pipeline

1. Deduplicate cleaned CDS sequences.
2. Select canonical PRNP reference and generate all single amino-acid substitutions.
3. Build Tier-2 biochemical and conservation features.
4. Optionally enrich Tier-2 output with AlphaFold-derived structural features.

## Files and Purpose

- `remdups.py`: deduplicate input CDS FASTA.
- `generate_prnp_mutations.py`: choose canonical PRNP reference and generate mutation table.
- `generate_prnp_mutation_features_tier2.py`: generate Tier-2 mutation-impact features.
- `generate_prnp_alphafold_features.py`: merge AlphaFold structural features into Tier-2 table.

## Reproducible Run Order

1. Deduplicate CDS:

```bash
python remdups.py
```

2. Generate reference protein and all single substitutions:

```bash
python generate_prnp_mutations.py \
	--input_fasta prnp_clean_nr.fasta \
	--output_csv prnp_single_point_mutations.csv \
	--reference_fasta prnp_reference_protein.fasta
```

3. Generate Tier-2 features:

```bash
python generate_prnp_mutation_features_tier2.py \
	--mutations_csv prnp_single_point_mutations.csv \
	--reference_fasta prnp_reference_protein.fasta \
	--source_cds_fasta prnp_clean_nr.fasta \
	--output_csv prnp_mutation_features_tier2.csv
```

## AlphaFold Integration (Beginner-Friendly)

### 1) Get AlphaFold structure

Download a PRNP AlphaFold model file (`.pdb` or `.cif`) and place it in this folder, for example:

- `AF-P04156-F1-model_v4.pdb`

### 2) Merge structural features into mutation table

```bash
python generate_prnp_alphafold_features.py \
	--af_structure AF-P04156-F1-model_v4.pdb \
	--reference_fasta prnp_reference_protein.fasta \
	--tier2_csv prnp_mutation_features_tier2.csv \
	--output_csv prnp_mutation_features_tier2_structural.csv \
	--position_csv prnp_alphafold_position_features.csv
```

### 3) What structural columns are added

The merge script adds per-position AlphaFold columns, then broadcasts them across all mutations at that position:

- `af2_plddt`: residue confidence score (from AlphaFold B-factor field).
- `af2_plddt_window_mean_5`, `af2_plddt_window_min_5`: local confidence context.
- `af2_neighbor_count_8A`, `af2_neighbor_count_12A`: C-alpha contact density.
- `af2_mean_neighbor_distance_8A`, `af2_min_neighbor_distance`: local packing geometry.
- `af2_is_buried_proxy`, `af2_is_exposed_proxy`: simple burial/exposure proxies from contact density.
- `af2_is_low_confidence`, `af2_is_high_confidence`: confidence flags.
- `af2_mapped`, `af2_struct_aa`, `af2_struct_seq_id`: sequence-to-structure mapping diagnostics.

## Notes

- `Bio.pairwise2` can emit a deprecation warning in newer Biopython versions; the current scripts still run.
- If multiple chains exist in the structure file, use `--chain` in the AlphaFold merge script.
