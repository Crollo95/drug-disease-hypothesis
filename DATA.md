## Data setup

This project **does not ship large datasets**. To run on real data, you must download a few
public datasets and place them under `data/raw/` in a specific layout.

You are responsible for complying with the licenses/terms of the source databases.

### Directory layout (expected by scripts)

After downloading, your `data/raw` directory should look like:

```text
data/raw/
├── chembl_activities/
│   ├── chembl_activities_part1.csv
│   ├── chembl_activities_part2.csv
│   ├── chembl_activities_part3.csv
│   └── chembl_activities_part4.csv
├── association_overall_direct/
│   ├── part-00000-...snappy.parquet
│   ├── part-00001-...snappy.parquet
│   └── ...
├── targets/
│   ├── part-00000-...snappy.parquet
│   ├── part-00001-...snappy.parquet
│   └── ...
├── molecular_interactions/
│   ├── part-00000-...snappy.parquet
│   ├── ...
├── drug_moa/
│   └── part-00000-...snappy.parquet
└── (optional) drug_indications/
    └── part-00000-...snappy.parquet
````

### Optional datasets (names / labels / evaluation)

The core DDH pipeline only requires:

- `chembl_activities/` (ChEMBL activities)
- `association_overall_direct/` (Open Targets Associations – direct overall)
- `targets/` (Open Targets Target)
- `molecular_interactions/` (Open Targets Molecular interactions)

In this repo we also use several **optional datasets** that improve interpretability
and enable ML evaluation.

#### Disease / Phenotype metadata

- Path: `data/raw/diseases/disease.parquet`
- Source: Open Targets “Disease / Phenotype”
- Use:
  - Enrich `data/real/diseases.csv` with better disease names, synonyms, ontology info.
- Not required for scoring, but recommended for nicer human-readable outputs.

#### Drug / clinical candidates metadata

- Path: `data/raw/ot_drugs/*.parquet`
- Source: Open Targets “Drug / Clinical candidates”
- Use:
  - Enrich `data/real/drugs.csv` with better drug names beyond bare ChEMBL IDs.
- Optional but very helpful for interpreting ranked drug lists.

#### Known drug (labels for ML & evaluation)

- Path: `data/raw/known_drug/*.parquet`
- Source: Open Targets “Known drug”
- Use:
  - Extract `(drug_id, disease_id)` pairs corresponding to approved / clinical indications.
  - Used to label candidate pairs with `label = 1` (known indication), which enables:
    - evaluating enrichment of known indications in the ranked list
    - training the frozen logistic ML model used in `ddh.ml_scoring`.
- **Recommended** if you want to use the ML ranking and evaluation; not strictly required for generating raw hypotheses.

#### Drug – indications (experimental)

- Path: `data/raw/drug_indications/*.parquet`
- Source: Open Targets “Drug – indications”
- Use:
  - We experimented with using this to derive known (drug, disease) pairs.
  - Due to its nested schema and lack of a simple disease ID field, it is **not currently used**
    in the main pipeline; the canonical source for labels here is the “Known drug” dataset instead.
- Safe to omit if you only care about the main pipeline.


### 1. Drug–target activities (ChEMBL)

1. Go to the ChEMBL downloads page and obtain the **activities** export (CSV).

2. Split files (if any) should be placed under:

   ```text
   data/raw/chembl_activities/chembl_activities_part*.csv
   ```

3. Merge these into a single CSV for convenience:

   ```bash
   cat data/raw/chembl_activities/chembl_activities_part*.csv \
     > data/raw/chembl_activities_merged.csv
   ```

4. Then run:

   ```bash
   python scripts/prepare_drug_targets_from_csv.py \
     --input data/raw/chembl_activities_merged.csv \
     --output-drug-targets data/real/drug_targets.csv \
     --output-drugs data/real/drugs.csv \
     --drug-id-col "Molecule ChEMBL ID" \
     --drug-name-col "Molecule Name" \
     --gene-id-col "Target Name" \
     --score-col "pChEMBL Value"
   ```

### 2. Gene–disease associations (Open Targets – Associations overall direct)

1. Download the **“Associations – direct (overall score)”** dataset (parquet).

2. Place all parquet files under:

   ```text
   data/raw/association_overall_direct/*.parquet
   ```

3. Convert to a flat CSV:

   ```bash
   python scripts/opentargets_associations_to_csv.py \
     --input-dir data/raw/association_overall_direct \
     --output data/raw/gene_disease_raw.csv
   ```

4. Normalize into canonical tables:

   ```bash
   python scripts/prepare_gene_disease_from_csv.py \
     --input data/raw/gene_disease_raw.csv \
     --output-gene-disease data/real/gene_disease.csv \
     --output-diseases data/real/diseases.csv \
     --gene-id-col gene_id \
     --disease-id-col disease_id \
     --disease-name-col disease_name \
     --score-col score
   ```

### 3. Target metadata (Open Targets – Target)

1. Download the **“Target”** dataset (parquet).

2. Place all parts under:

   ```text
   data/raw/targets/*.parquet
   ```

3. Build a mapping from various target IDs to Ensembl:

   ```bash
   python scripts/opentargets_targets_to_mapping.py \
     --input-dir data/raw/targets \
     --output data/real/target_mapping.csv
   ```

4. Use this to map the ChEMBL target names to Ensembl IDs:

   ```bash
   python scripts/map_drug_targets_to_ensembl.py
   ```

   This will write `data/real/drug_targets_mapped.csv`.

### 4. Molecular interactions (Open Targets – Molecular interactions)

1. Download the **“Molecular interactions”** dataset (parquet).

2. Place all parts under:

   ```text
   data/raw/molecular_interactions/*.parquet
   ```

3. Convert to a simple PPI CSV:

   ```bash
   python scripts/opentargets_molecular_interactions_to_ppi.py \
     --input-dir data/raw/molecular_interactions \
     --output data/real/ppi.csv \
     --gene1-col targetA \
     --gene2-col targetB \
     --weight-col scoring \
     --min-weight 0.4
   ```

### 5. Mechanism-of-action (MoA) (Open Targets – Drug mechanism of action)

1. Download the **“Drug – mechanism of action”** dataset (parquet).

2. Place it under:

   ```text
   data/raw/drug_moa/*.parquet
   ```

3. Extract `(drug_id, gene_id)` MoA pairs:

   ```bash
   python scripts/opentargets_moa_to_pairs.py \
     --input-dir data/raw/drug_moa \
     --target-mapping data/real/target_mapping.csv \
     --output-raw data/real/moa_pairs_raw.csv \
     --output-ensembl data/real/moa_pairs_ensembl.csv
   ```

4. (Optional) Filter drug–target edges to MoA-backed edges:

   ```bash
   python scripts/filter_drug_targets_by_moa.py \
     --drug-targets data/real/drug_targets_mapped.csv \
     --moa-pairs data/real/moa_pairs_ensembl.csv \
     --output data/real/drug_targets_moa.csv
   ```

### 6. Filtering & gene universe

To apply basic strength filters and build a unified gene list:

```bash
python scripts/filter_real_data.py
python scripts/prepare_genes_from_real_data.py
```

This will create:

* `data/real/drug_targets_filtered.csv`
* `data/real/gene_disease_filtered.csv`
* `data/real/genes.csv`

### 7. PPI distance matrix (required for fast proximity scoring)

Finally, precompute the gene–gene distance matrix:

```bash
python scripts/precompute_gene_distance_matrix.py \
  --ppi-csv data/real/ppi.csv \
  --drug-targets-csv data/real/drug_targets_filtered.csv \
  --gene-disease-csv data/real/gene_disease_filtered.csv \
  --out-index data/real/gene_index.csv \
  --out-matrix data/real/gene_distances.uint16.dat \
  --cutoff 6
```

At this point, the real-data pipeline can be run via:

```bash
python -m ddh.cli csv \
  --drugs-csv data/real/drugs.csv \
  --genes-csv data/real/genes.csv \
  --diseases-csv data/real/diseases.csv \
  --drug-targets-csv data/real/drug_targets_filtered.csv \
  --gene-disease-csv data/real/gene_disease_filtered.csv \
  --ppi-csv data/real/ppi.csv \
  --top-k 1000000 \
  --output-csv outputs/real_results_with_ppi.csv
```

