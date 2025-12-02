# Drug–Disease–Gene Hypothesis Pipeline (DDH)

A **local, scalable pipeline** to generate and rank:

* **Drug → Gene (target) hypotheses**
* **Drug → Disease repurposing hypotheses**

Designed for **researchers and developers** in healthcare, biology, and computational sciences who want to run everything **locally** with minimal dependencies, maximum control, reproducibility, and cluster-scale scalability.

---

## 🚀 Overview

This project integrates:

| Resource                                                     | Role                                                                           |
| ------------------------------------------------------------ | ------------------------------------------------------------------------------ |
| **Drug–Target associations** (ChEMBL-derived `pChEMBL`)      | Maps drugs to human gene targets                                               |
| **Gene–Disease associations** (Open Targets genetics scores) | Maps diseases to human risk genes                                              |
| **Protein–Protein Interaction (PPI) network**                | Computes gene–gene proximity across the interactome                            |
| **Combined scoring + ML ranker**                             | Prioritizes hypotheses with overlap, network proximity, and MoA-aware features |

All computation can run on:

* A **local workstation**, or
* A **cluster** (tested on 24 cores, 126 GB RAM), supporting multi-day jobs.

---

## 📦 Installation

```bash
pip install -e .
```

This installs the `ddh` package for import and provides a CLI under `python -m ddh.cli`.

---

# ✅ Before You Start: Validate Your Raw Data

DDH relies on several large external datasets (Open Targets, ChEMBL, MoA datasets, etc.).
Before running any pipeline steps, **ensure your `data/raw/` directory is correctly populated and has the expected schema**.

We provide an automated validator:

```bash
python scripts/check_raw_data.py
```

This script:

* Ensures all required raw datasets exist
* Detects missing or misnamed parquet shards
* Identifies malformed CSVs (e.g., ChEMBL semicolon-delimited exports)
* Verifies presence of required columns
* Produces clear ❌ errors / ⚠️ warnings when something is off

When everything is correct you will see:

```
✅ All required raw dataset checks passed.
```

Fix any issues reported before proceeding with data preparation.

---

## 🧠 Design Choices

### 1. Local-first

* No hidden online lookups
* All inputs downloaded manually / via scripts
* Ensures reproducibility, privacy, cluster usability

### 2. Conservative filtering

* Drug–target: `pChEMBL ≥ 5.0`
* Gene–disease: `genetics score ≥ 0.01`
* Encourages mechanistically plausible hypotheses

### 3. PPI precomputation

* We compute all **gene–gene distances once** (multi-hour)
* Store in a **uint16 memmap** for fast lookups across millions of pairs

### 4. Overlap-first candidate space

* We avoid full O(N_drugs × N_diseases) cartesian (billions)
* Instead generate ~40M biologically plausible overlap-derived candidates
* Then refine via PPI + ML

---

## ⚗️ What the Pipeline Provides

1. **Data ingestion & normalization**

   * ChEMBL drug–targets → canonical formats
   * Open Targets associations → flat CSVs
   * Ontology-aware disease table
   * Drug MoA mapping (ENSG targets)

2. **Feature computation**

   * Overlap stats
   * Graph-based proximity
   * Target and disease degrees
   * MoA features
   * ML feature engineering

3. **Distance matrix**

   * Precomputes all-pairs shortest paths
   * ~100–400 MB memmap
   * Reusable across runs

4. **High-throughput scoring**

   * Score tens of millions of (drug, disease) pairs
   * PPI proximity + logistic model
   * Cluster parallelism supported

5. **Interpretability utilities**

   * Inspect top-ranked repurposing hypotheses for a disease
   * View overlap genes, distances, MoA support

---

## 🧬 Running the Pipeline

### 1. Toy pipeline

```bash
python -m ddh.cli --toy-data-dir data/toy --top-k 10
```

### 2. Prepare real data (through the scripts in `scripts/`)

Once `check_raw_data.py` passes, run each step of the data-normalization pipeline.

### 3. Overlap-only baseline

```bash
python -m ddh.cli csv \
  --drugs-csv data/real/drugs.csv \
  --genes-csv data/real/genes.csv \
  --diseases-csv data/real/diseases.csv \
  --drug-targets-csv data/real/drug_targets_filtered.csv \
  --gene-disease-csv data/real/gene_disease_filtered.csv \
  --top-k 1000000 \
  --output-csv outputs/real_results_overlap.csv
```

### 4. Add PPI-based proximity

```bash
python scripts/annotate_pairs_with_ppi_from_matrix.py \
  --pairs-csv outputs/real_results_overlap.csv \
  --drug-targets-csv data/real/drug_targets_filtered.csv \
  --gene-disease-csv data/real/gene_disease_filtered.csv \
  --gene-index data/real/gene_index.csv \
  --dist-matrix data/real/gene_distances.uint16.dat \
  --max-pairs 100000000 \
  --chunk-size 10000 \
  --output outputs/real_results_with_ppi.csv
```

### 5. ML scoring (MoA-aware logistic model)

Load model via:

```python
from ddh.ml_scoring import score_pairs_with_frozen_moa_model
```

Call the per-disease inspector:

```python
df_ranked = show_top_drugs_for_disease(df_pairs, "EFO_XXXXX")
```

---

## 🧪 Tests

Run:

```bash
pytest
```

Tests include:

* PPI graph logic
* Overlap scoring
* Proximity scoring
* CLIs
* ML feature sanity checks

---

## 🔮 Future Directions

* PPI-only and diffusion-based candidate generation
* MONDO/EFO disease unification
* Improved MoA weighting
* Graph neural networks for repurposing
* Ontology-aware disease clustering
* Explainability + graph visualizations

---

## 🤝 Contributions

PRs are welcome!
Please run:

```bash
python scripts/check_raw_data.py
pytest
```

before contributing.

