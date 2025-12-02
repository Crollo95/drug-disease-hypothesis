
# Methods / System Overview — drug‑disease‑hypothesis (DDH) pipeline



## Purpose



The DDH pipeline is designed to generate and prioritize **drug–disease hypotheses** by integrating public large-scale data on drug–target binding, gene–disease associations, and protein–protein interactions (PPI), with a downstream machine-learning ranking module. The goal is to identify promising (and under-explored) repurposing candidates in a reproducible, scalable, transparent, and local-first manner.



The workflow aims to balance biological plausibility (target overlap / network proximity / MoA evidence) with computational tractability (avoid combinatorial explosion) and reproducibility, enabling hypothesis generation for preclinical or translational follow-up.



---



## Overview of Pipeline Steps



Broadly, the pipeline consists of the following stages:



1. **Data ingestion & normalization**

2. **Filtering and mapping to unified identifiers**

3. **Graph / network construction (PPI)**

4. **Computing per-pair features: overlap, network proximity, coverage fractions**

5. **(Optional) MoA anchoring — high-confidence curated drug–target edges**

6. **Ranking via logistic regression (frozen model)**

7. **Output and hypothesis export (CSV / notebook / CLI)**



Below is a more detailed description of each stage.



---



## 1. Data ingestion & normalization



* **Drug–target data**



  * Source: public databases (e.g. ChEMBL activities) containing biochemical assay results (IC₅₀, Kᵢ, K_d, pChEMBL, etc.)

  * All portions of the raw data (if split across parts) are merged before normalization.

  * Only one activity metric (e.g. `pChEMBL`) is used for scoring; weak bindings below a threshold are filtered out.



* **Gene–disease data**



  * Source: large-scale target–disease association datasets (e.g. from publicly available platforms integrating genetics, GWAS, functional evidence)

  * Each record maps a `gene_id` (e.g. Ensembl) to a `disease_id` (e.g. EFO, MONDO, Orphanet), plus a numeric confidence/association score.

  * The script normalizes and filters by a minimal confidence threshold (e.g. ≥ 0.01), to discard weak or spurious links.



* **PPI (protein–protein interactions)**



  * Source: integrated datasets covering physical, functional, or pathway-based interactions (multiple databases).

  * Interaction edges are filtered by a minimal weight (e.g. confidence score ≥ 0.4).

  * The resulting PPI graph connects genes/proteins in a large network to be used for network-based proximity scoring.



* **Identifier mapping / canonicalization**



  * Drugs: originally from ChEMBL chemical IDs; kept as canonical drug IDs.

  * Genes / Targets: mapped to stable identifiers (e.g. Ensembl `gene_id`) to unify across datasets.

  * Diseases: disease IDs (EFO, MONDO, Orphanet, etc.) are retained as given; user should be aware of duplication or overlap across ontologies (disease clustering recommended downstream).



---



## 2. Filtering and cleaning



* **Drug–target binding filter**: only target–gene associations with sufficiently strong binding (e.g. pChEMBL ≥ threshold) are kept.

* **Gene–disease filter**: only gene–disease associations with confidence ≥ threshold (e.g. 0.01) are kept.

* **Optional “MoA anchoring”**: from a curated drug–mechanism-of-action dataset, extract high-confidence `(drug, gene)` pairs; treat these as “ground-truth” targets for drugs where available.



These filters help reduce noise and focus on biologically plausible and therapeutically relevant connections.



---



## 3. Graph construction & distance matrix precomputation



Using the filtered PPI data:



* Build a **graph** where nodes are genes/proteins and edges are PPIs (weighted by confidence).

* Precompute a **distance (shortest-path) matrix** between all genes in the unified gene list.



  * Stored in an efficient binary format (e.g. `uint16`) for memory and speed.

  * Enables fast lookup of pairwise network distances, avoiding repeated expensive graph searches.



This setup provides the backbone for subsequent network-proximity scoring between drug targets and disease genes.



---



## 4. Feature computation for drug–disease pairs



For each candidate pair `(drug, disease)`, the pipeline computes:



* **Overlap-based features**:



  * `n_overlap`: number of genes common to the drug’s target set and the disease gene set.

  * `log1p_n_overlap`: log-transformed overlap (log(1 + n)).

  * `drug_deg`, `disease_deg`: number of unique targets per drug; number of disease-associated genes per disease.

  * `frac_drug_covered`: fraction of drug’s targets involved in the disease.

  * `frac_disease_covered`: fraction of disease-associated genes targeted by the drug.



* **Network-proximity feature**:



  * `ppi_proximity`: derived from the mean shortest path distance between drug targets and disease genes (via the precomputed distance matrix), transformed into a monotonic proximity score.



* **MoA-based features** (if MoA data exists):



  * `n_moa_targets`: number of curated, high-confidence targets for the drug.

  * `drug_has_moa`: binary flag (1 if at least one MoA target, else 0).



These features capture both **direct mechanistic overlap** (binding), **network neighborhood**, and **curated pharmacological relevance**.



---



## 5. Ranking via a frozen ML model



* A **logistic regression model** was trained on a sampled subset of drug–disease pairs, mixing positives (known indications) and a controlled number of negatives, using the full feature set (overlap + proximity + MoA).

* The model is then **frozen**: the scaler parameters (mean & scale), coefficients, and intercept are stored in `ml_scoring.py`.

* A utility function `score_pairs_with_frozen_moa_model(df)` allows applying the trained model to any table of candidate pairs, computing a normalized probability-like score (`ml_score_moa`) and enabling ranking.



This design provides a reproducible, stable ranking mechanism: you do not need to retrain the model every time — the same ranking logic can be applied to new data or parameter configurations.



---



## 6. Output, usability & modular design



* The pipeline is implemented as a **Python package** (`ddh/`) plus **scripts** (in `scripts/`) for data preparation, mapping, filtering, PPI construction, matrix precomputation, annotation, and ranking.

* A **CLI interface** allows chaining steps via command-line commands, enabling non-interactive use.

* Modular structure, scripted transformations, and test coverage (basic for graph construction and scoring) ensure maintainability.

* The **frozen ML model + minimal dependencies** makes deployment and reuse easier for collaborators.



---



## Design choices and trade-offs (rationale)



* **Overlap-first selection, then proximity** — avoids combinatorial explosion (billions of drug–disease pairs) and focuses on biologically plausible hypotheses; saves computational resources.

* **Use of public data only** — increases reproducibility and accessibility, and avoids licensing issues.

* **Precomputed distance matrix** — speeds up repeated proximity calculations, at cost of memory/time in precompute but huge savings in downstream bulk annotation.

* **MoA anchoring as optional feature** — accommodates both high-confidence curated targets and broad binding data; avoids over-filtering, while giving strong signal to known pharmacological drug–gene links.

* **Frozen ML ranking** — makes results reproducible, sharable, and auditable; avoids “black-box” retraining with unknown variability.



These design principles reflect broader best practices in computational pipelines: modularity, reproducibility, automation, transparency. ([OUP Academic][1])



---



## Known limitations and caveats



* **Disease identifier duplication / ontology fragmentation** — many disease IDs may refer to overlapping or similar disease constructs (e.g. different ontologies, subtypes); without disease clustering, the ranked output may have redundant hypotheses.

* **Binding data limitations** — biochemical assays (e.g. from ChEMBL) may not translate into effective in-vivo activity; affinity thresholds are arbitrary, and functional relevance (cell type, tissue, expression, PK/PD) is not considered.

* **PPI network context-blindness** — PPIs are aggregated across contexts (tissues, conditions) and do not reflect tissue-specific expression, pathway activation, or cell-type specificity.

* **No tissue / expression / pharmacokinetics filtering** — the model does not enforce that drug targets are expressed in disease-relevant tissues, nor checks PK/ADMET feasibility.

* **Sparse MoA coverage** — many drugs lack curated MoA data; their ranking relies solely on binding + network proximity, which may be less reliable.

* **No downstream biological validation** — predictions remain computational hypotheses; experimental / clinical follow-up required.

* **Scalability limitations** — although optimized vs full cartesian search, very large-scale analyses (e.g., all ~40M pairs + repeated reranking) still require significant computational resources.



---



## Summary



The DDH pipeline implements a **transparent, reproducible, data-driven drug repurposing engine**. By combining:



* drug–target binding data,

* gene–disease genetic and functional associations,

* protein–protein interaction network topology, and

* a simple but effective ML ranking model with MoA anchoring,



it generates ranked hypotheses for drug repurposing that balance **mechanistic plausibility**, **network context**, and **therapeutic relevance** — while being easy to rerun, extend, and share.



Given its modular design, open-source implementation, and public data foundation, DDH can be used as a **baseline engine** for further extensions (e.g. tissue-specific filtering, population-level data, PPI-only hypotheses, etc.), or as a starting point for experimental prioritization or translational repurposing efforts.



---



## Recommended next steps



To increase the utility, reproducibility, and shareability of the project, we recommend:



* freezing the current version as a release (with metadata)

* adding data-download instructions / manifest + schema checks for raw inputs

* implementing disease clustering / ontology canonicalization

* extending documentation (usage examples, limitations, future directions)

* optionally integrating downstream filtering (e.g. tissue expression, PK, MoA context)



Once voiced, these would turn DDH into a **fully reproducible, sharable, and community-ready drug-disease hypothesis engine**.






No file chosenNo file chosen
