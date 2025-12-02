import pandas as pd
import numpy as np
from ddh.ml_scoring import score_pairs_with_frozen_moa_model, FEATURE_COLS_MOA

def test_frozen_model_basic():
    df = pd.DataFrame({
        "log1p_n_overlap": [0.0, 1.0],
        "drug_deg": [1, 10],
        "disease_deg": [3, 50],
        "frac_drug_covered": [0.0, 0.1],
        "frac_disease_covered": [0.0, 0.05],
        "ppi_proximity": [0.1, 0.3],
        "n_moa_targets": [0, 2],
        "drug_has_moa": [0, 1],
    })
    df_scored = score_pairs_with_frozen_moa_model(df)
    assert "ml_score_moa" in df_scored.columns
    assert np.all(np.isfinite(df_scored["ml_score_moa"]))
    assert ((df_scored["ml_score_moa"] >= 0) & (df_scored["ml_score_moa"] <= 1)).all()

