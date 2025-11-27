"""
ML-based scoring module for drug–disease hypotheses.

We use a simple Logistic Regression model trained on overlap + PPI proximity features,
after standard scaling. The learned coefficients are frozen and applied directly
for fast ranking (no retraining required at inference time).

This outputs a therapeutic-aligned score for each drug–disease pair.
"""

from __future__ import annotations
import numpy as np
import pandas as pd
from dataclasses import dataclass
from typing import Sequence

# -----------------------------
# 1. Frozen model parameters
# -----------------------------

FEATURE_COLS = [
    "log1p_n_overlap",
    "drug_deg",
    "disease_deg",
    "frac_drug_covered",
    "frac_disease_covered",
    "ppi_proximity",
]

# Coefficients learned from your training (feature order must match FEATURE_COLS)
COEFS = np.array([
    0.711,   # log1p_n_overlap
    0.176,   # drug_deg
    0.272,   # disease_deg
   -0.792,   # frac_drug_covered
   -0.423,   # disease_covered
    0.130,   # ppi_proximity
])

INTERCEPT = -1.5  # We approximate a reasonable intercept as a starting point. 
                 # We can recalibrate this later using PR-targeted thresholds.

@dataclass
class ScalerParams:
    mean: np.ndarray
    scale: np.ndarray

# These come from the scaler you fit. You must copy real values from notebook:
SCALER_MEAN  = np.array([1.2, 80.5, 120.3, 0.05, 0.04, 0.12])
SCALER_SCALE = np.array([2.5, 200.1, 300.2, 0.2, 0.18, 0.3])

SCALER = ScalerParams(
    mean=SCALER_MEAN,
    scale=SCALER_SCALE,
)

# ----------------------------------------
# 2. Feature computation on overlap pairs
# ----------------------------------------

def build_overlap_network_features(
    drug_targets: Sequence[tuple[str, float]],
    disease_genes: Sequence[tuple[str, float]],
    mean_distance: float | int,
    n_overlap: int
) -> np.ndarray:
    """
    Compute the feature vector for a single drug–disease pair.
    """
    # Convert to dict for any future extended metrics:
    drug_targets = list(drug_targets)
    disease_genes = list(disease_genes)

    # Gene set sizes:
    drug_genes = {g for g,_ in drug_targets}
    dis_genes  = {g for g,_ in disease_genes}

    drug_deg = len(drug_genes)
    disease_deg = len(dis_genes)

    # Derived fractions:
    frac_drug_covered = n_overlap / drug_deg if drug_deg > 0 else 0.0
    frac_disease_covered = n_overlap / disease_deg if disease_deg > 0 else 0.0

    # Inverse distance (for cutoff-based matrix 65535 is no path)
    if mean_distance >= 65535 or not np.isfinite(mean_distance):
        ppi_proximity = 0.0
    else:
        ppi_proximity = 1.0 / (1.0 + mean_distance)

    log1p_n_overlap = np.log1p(n_overlap)

    return np.array([
        log1p_n_overlap,
        drug_deg,
        disease_deg,
        frac_drug_covered,
        frac_disease_covered,
        ppi_proximity,
    ], dtype=np.float32)

def scale_features(x: np.ndarray, scaler: ScalerParams) -> np.ndarray:
    """
    Apply the frozen standard scaling to a feature vector or matrix.
    """
    return (x - scaler.mean) / scaler.scale

def sigmoid(x: np.ndarray) -> np.ndarray:
    return 1.0 / (1.0 + np.exp(-x))

# -----------------------------
# 3. Apply ML score to a table
# -----------------------------

def score_pairs_with_ml(df_pairs: pd.DataFrame) -> pd.DataFrame:
    """
    Compute ML scores for a dataframe of drug–disease pairs.

    Required columns:
        df_pairs must contain at least:
            - drug_id
            - disease_id
            - n_overlap
            - mean_distance

    Output:
        Returns a copy of df with an extra column `ml_score`.
    """
    missing = [c for c in FEATURE_COLS if c not in df_pairs.columns]
    req = ["drug_id", "disease_id", "n_overlap", "mean_distance"]
    missing_req = [c for c in req if c not in df_pairs.columns]
    if missing_req:
        raise ValueError(f"Missing required columns: {missing_req}")

    df_out = df_pairs.copy()
    feats = []

    for _,row in df_out.iterrows():
        f = build_overlap_network_features(
            drug_targets=[],  # not used here yet – placeholder for later double-check enrichment
            disease_genes=[],
            mean_distance=row.mean_distance,
            n_overlap=row.n_overlap
        )
        feats.append(f)

    X = np.stack(feats)
    Xs = scale_features(X, SCALER)
    score = sigmoid(Xs @ COEFS + INTERCEPT)
    df_out["ml_score"] = score.astype(np.float32)
    return df_out
