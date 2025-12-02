from __future__ import annotations

import argparse
from pathlib import Path
from typing import Optional

import pandas as pd


def parse_args(argv: Optional[list[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Normalize a raw gene–disease CSV into DDH canonical format.",
    )

    parser.add_argument(
        "--input",
        type=str,
        required=True,
        help="Path to the raw gene–disease CSV file.",
    )
    parser.add_argument(
        "--output-gene-disease",
        type=str,
        required=True,
        help="Path to write the normalized gene_disease.csv file.",
    )
    parser.add_argument(
        "--output-diseases",
        type=str,
        required=True,
        help="Path to write the normalized diseases.csv file.",
    )

    parser.add_argument(
        "--gene-id-col",
        type=str,
        required=True,
        help="Name of the column containing gene IDs/symbols in the raw CSV.",
    )
    parser.add_argument(
        "--disease-id-col",
        type=str,
        required=True,
        help="Name of the column containing disease IDs in the raw CSV.",
    )
    parser.add_argument(
        "--disease-name-col",
        type=str,
        required=True,
        help="Name of the column containing disease names in the raw CSV.",
    )
    parser.add_argument(
        "--score-col",
        type=str,
        default=None,
        help="Optional: name of the column containing a numeric association score.",
    )

    return parser.parse_args(argv)


def main(argv: Optional[list[str]] = None) -> None:
    args = parse_args(argv)

    input_path = Path(args.input)
    if not input_path.exists():
        raise SystemExit(f"Input file does not exist: {input_path}")

    df_raw = pd.read_csv(
        input_path,
        sep=None,
        engine="python",
        on_bad_lines="skip",
    )

    required_cols = {args.gene_id_col, args.disease_id_col, args.disease_name_col}
    missing = required_cols - set(df_raw.columns)
    if missing:
        raise SystemExit(f"Missing required columns in input: {missing}")

    # Build normalized gene_disease table
    out_cols = {
        "gene_id": args.gene_id_col,
        "disease_id": args.disease_id_col,
    }
    df_gd = df_raw[list(out_cols.values())].rename(columns={v: k for k, v in out_cols.items()})

    if args.score_col is not None and args.score_col in df_raw.columns:
        df_gd["score"] = df_raw[args.score_col]

    df_gd = df_gd.dropna(subset=["gene_id", "disease_id"])
    df_gd["gene_id"] = df_gd["gene_id"].astype(str).str.strip()
    df_gd["disease_id"] = df_gd["disease_id"].astype(str).str.strip()

    if "score" in df_gd.columns:
        df_gd = (
            df_gd.groupby(["gene_id", "disease_id"], as_index=False)
            .agg({"score": "max"})
        )
    else:
        df_gd = df_gd.drop_duplicates(subset=["gene_id", "disease_id"])

    # Build normalized diseases table
    df_diseases = (
        df_raw[[args.disease_id_col, args.disease_name_col]]
        .rename(columns={
            args.disease_id_col: "disease_id",
            args.disease_name_col: "disease_name",
        })
        .dropna(subset=["disease_id"])
        .drop_duplicates(subset=["disease_id"])
        .sort_values(by=["disease_id"])
        .reset_index(drop=True)
    )

    df_diseases["disease_id"] = df_diseases["disease_id"].astype(str).str.strip()
    df_diseases["disease_name"] = df_diseases["disease_name"].astype(str).str.strip()

    out_gd = Path(args.output_gene_disease)
    out_dis = Path(args.output_diseases)
    out_gd.parent.mkdir(parents=True, exist_ok=True)
    out_dis.parent.mkdir(parents=True, exist_ok=True)

    df_gd.to_csv(out_gd, index=False)
    df_diseases.to_csv(out_dis, index=False)

    print(f"Normalized gene_disease written to: {out_gd.resolve()}")
    print(f"Normalized diseases written to: {out_dis.resolve()}")


if __name__ == "__main__":
    main()

