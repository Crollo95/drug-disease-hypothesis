# Drug–Disease–Gene Hypothesis Pipeline (DDH)

A **local, scalable pipeline** to generate and rank:
- **Drug → Gene (target) hypotheses**
- **Drug → Disease repurposing hypotheses**

Designed for **researchers and developers** in healthcare, biology, and computational sciences who want to run everything **locally** with minimal dependencies, maximum control, reproducibility, and cluster-scale scalability.

---

## 🚀 Overview

This project implements a computational pipeline that integrates:

| Resource | Role |
|---|---|
| **Drug–Target associations** (e.g., ChEMBL-derived `pChEMBL`) | Maps drugs to human gene targets |
| **Gene–Disease associations** (e.g., Open Targets genetics scores) | Maps diseases to human risk genes |
| **Protein/Gene interaction network (PPI)** | Computes proximity between drug and disease genes |
| **Combined scoring scheme** | Ranks biased overlap hits + network-based mechanistic links |

All computation can run on:
- A **local machine**, or
- A **computational cluster** (tested on 24 cores, 126 GB RAM), supporting **multi-day jobs**.

---

## ✅ What We Have Built

1. `ddh/` is a valid **installable Python package** (via `pyproject.toml`)
2. A **toy pipeline** for quick development iteration
3. A **production prototype pipeline** that:
   - Filters weak signals (`score ≥ 5.0` for drug targets, `score ≥ 0.01` for gene-disease)
   - Precomputes all **gene–gene network distances** once into a `uint16` matrix for reuse
   - Supports a CLI (`python -m ddh.cli csv ...`) that runs ingestion + overlap + proximity + scoring
4. A standalone **cluster-ready PPI annotator** to score large batches of drug–disease pairs using matrix lookups (`annotate_pairs_with_ppi_from_matrix.py`)

---

## 🧠 Design Choices

### 1. **Local-first philosophy**
- Avoid unnecessary online execution
- Download external data manually or via user scripts, normalize, and operate locally
- Ensures privacy, reproducibility, and auditability

### 2. **Conservative scoring by default**
- Strong **drug–gene binding** required (`pChEMBL ≥ 5.0`)
- Weak genetic evidence is excluded (`gene–disease score ≥ 0.01`)
- This favors:
  - Interpretability
  - Clinically plausible repurposings
  - Genetic support for disease biology

### 3. **PPI proximity from scratch, but cached**
- We precompute **single-source BFS distances per gene once** into a `uint16 memmap matrix`
  - This allows:
    - Fast lookup per pair (no repeated BFS)
    - Multi-process parallel reads
    - Feasible scalability to tens or hundreds of millions of pairs

### 4. **Initial candidate class = overlap space**
- The first full prototype annotates only **overlap-driven drug–disease pairs with network refinement**, which is:
  - Not scientifically incorrect, but **more conservative**
- You can later expand to:
  - PPI-only pairs
  - Diffusion-based or embedding-based network medicine metrics

---

## ⚠️ Caveats, Nuances & Practical Notes

- The full cartesian drug×disease space (billions of pairs) is **not directly processed**
  - Instead, overlap generation yields ~40M biologically plausible pairs
  - Annotating 40M from the distance matrix takes ~10h on 1 core, proving feasibility of cluster scale
- This is a **prototype research pipeline**, not a validated clinical decision system
- pChEMBL scores are numeric but not necessarily comparable across target classes
  - Future improvement may include per-assay normalization or confidence weighting
- Disease IDs are currently based on **EFO ontology**
  - Meaningful branch-restriction filters may help in the future (e.g., oncology-focused, immune-focused)
- The PPI graph edges are considered undirected and deduplicated
  - Future work could incorporate directionality (SIGNOR, pathway edges) or multi-edge weighting

---

## 🧬 Example Usage

### Install package locally

```bash
pip install -e .
```

### Run toy pipeline

```bash
python -m ddh.cli --toy-data-dir data/toy --top-k 10
```

### Run real overlap-based pipeline (no PPI)

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

### Add cluster-scale PPI proximity (subset of interest)

```bash
python scripts/annotate_pairs_with_ppi_from_matrix.py \
  --pairs-csv outputs/real_results_overlap.csv \
  --drug-targets-csv data/real/drug_targets_filtered.csv \
  --gene-disease-csv data/real/gene_disease_filtered.csv \
  --gene-index data/real/gene_index.csv \
  --dist-matrix data/real/gene_distances.uint16.dat \
  --max-pairs 100000000 \
  --chunk-size 10000 \
  --alpha 1.0 \
  --beta 1.0 \
  --output outputs/real_results_with_ppi.csv
```

### 🧪 Tests

Run:
```bash
pytest
```
All core components have passing unit tests:
- PPI graph building
- Overlap scoring
- Proximity scoring
- CLI execution chain


## 🔮 Future Improvements (Roadmap)

## 🤝 Contributions





