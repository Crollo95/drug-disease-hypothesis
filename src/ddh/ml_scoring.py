"""
ML-based scoring for drug–disease hypotheses (with MoA features).

This module provides a *frozen* logistic regression scorer:

    - trained on overlap + PPI + MoA features
    - standard scaling applied to features
    - outputs a probability-like score `ml_score_moa`

The idea is that the model is trained once in a notebook, and the
resulting parameters (feature order, scaler mean/scale, coefficients,
intercept) are stored here for fast, reproducible inference.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import List

import numpy as np
import pandas as pd


# --------------------------------------------------------------------
# 1. Frozen model definition
# --------------------------------------------------------------------

# IMPORTANT: these must match exactly the order used in the notebook
FEATURE_COLS_MOA: List[str] = [
    "log1p_n_overlap",
    "drug_deg",
    "disease_deg",
    "frac_drug_covered",
    "frac_disease_covered",
    "ppi_proximity",
    "n_moa_targets",
    "drug_has_moa",
]

# Paste your actual values from the notebook here
SCALER_MEAN = np.array([np.float64(0.76770570366969), np.float64(2.1602803137456603), np.float64(866.3386160044575), np.float64(0.8375434231212616), np.float64(0.032570881299623856), np.float64(0.34160482839567646), np.float64(0.39023188033089024), np.float64(0.10011358278685012)], dtype=np.float64)

SCALER_SCALE = np.array([np.float64(0.2411371264991879), np.float64(4.125444617342613), np.float64(1346.4407245197522), np.float64(0.27842939222196394), np.float64(0.10527375044361369), np.float64(0.0680506756230397), np.float64(2.1653908388852714), np.float64(0.30015138401884917)], dtype=np.float64)

COEFS = np.array([np.float64(0.5606596061189063), np.float64(-0.17859301222661653), np.float64(0.2614834850038309), np.float64(0.03209003183540381), np.float64(-0.3931886376035342), np.float64(0.44448249320154837), np.float64(0.030300249073566703), np.float64(4.429415288172351)], dtype=np.float64)

INTERCEPT = float(
    -8.717857947648007
)

@dataclass
class Scaler:
    mean: np.ndarray
    scale: np.ndarray

    def transform(self, X: np.ndarray) -> np.ndarray:
        return (X - self.mean) / self.scale


SCALER = Scaler(mean=SCALER_MEAN, scale=SCALER_SCALE)


def _sigmoid(x: np.ndarray) -> np.ndarray:
    return 1.0 / (1.0 + np.exp(-x))


# --------------------------------------------------------------------
# 2. Public API: score a DataFrame
# --------------------------------------------------------------------

def score_pairs_with_frozen_moa_model(
    df: pd.DataFrame,
    score_col: str = "ml_score_moa",
) -> pd.DataFrame:
    """
    Score drug–disease pairs using the frozen logistic+scaling model.

    Parameters
    ----------
    df : DataFrame
        Must contain the feature columns in FEATURE_COLS_MOA:
        - log1p_n_overlap
        - drug_deg
        - disease_deg
        - frac_drug_covered
        - frac_disease_covered
        - ppi_proximity
        - n_moa_targets
        - drug_has_moa

        These features are expected to be *already computed* at the
        per-pair level (e.g. by your pipeline / notebook code).

    score_col : str
        Name of the column to write the ML score into.

    Returns
    -------
    DataFrame
        Copy of the input df with one extra column containing the
        model score in [0, 1].
    """
    missing = [c for c in FEATURE_COLS_MOA if c not in df.columns]
    if missing:
        raise ValueError(f"Missing required feature columns: {missing}")

    X = df[FEATURE_COLS_MOA].to_numpy(dtype=float)
    Xs = SCALER.transform(X)

    logits = Xs @ COEFS + INTERCEPT
    probs = _sigmoid(logits)

    df_out = df.copy()
    df_out[score_col] = probs.astype(np.float32)
    return df_out

