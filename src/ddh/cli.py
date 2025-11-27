# src/ddh/cli.py

from __future__ import annotations

import argparse
from pathlib import Path
import sys

import pandas as pd

from .pipeline import run_toy_pipeline, run_csv_pipeline
from .config import CsvFilesConfig


def add_toy_subcommand(subparsers: argparse._SubParsersAction) -> None:
    parser_toy = subparsers.add_parser(
        "toy",
        help="Run the pipeline on the built-in toy dataset.",
        description="Run the drug–disease hypothesis pipeline on toy CSV data.",
    )

    parser_toy.add_argument(
        "--toy-data-dir",
        type=str,
        default="data/toy",
        help="Directory containing toy CSV files (default: data/toy).",
    )
    parser_toy.add_argument(
        "--alpha",
        type=float,
        default=1.0,
        help="Weight for normalized overlap term in combined score (default: 1.0).",
    )
    parser_toy.add_argument(
        "--beta",
        type=float,
        default=1.0,
        help="Weight for proximity term in combined score (default: 1.0).",
    )
    parser_toy.add_argument(
        "--top-k",
        type=int,
        default=20,
        help="Number of top-ranked pairs to print (default: 20).",
    )
    parser_toy.add_argument(
        "--output-csv",
        type=str,
        default=None,
        help="If provided, save the full ranking to this CSV file.",
    )

    parser_toy.set_defaults(func=run_toy_cli)


def add_csv_subcommand(subparsers: argparse._SubParsersAction) -> None:
    parser_csv = subparsers.add_parser(
        "csv",
        help="Run the pipeline on user-provided CSV files.",
        description="Run the drug–disease hypothesis pipeline on custom CSV input.",
    )

    # Required file paths
    parser_csv.add_argument(
        "--drugs-csv",
        type=str,
        required=True,
        help="CSV file with drug metadata (id + name).",
    )
    parser_csv.add_argument(
        "--genes-csv",
        type=str,
        required=True,
        help="CSV file with gene metadata (id + symbol).",
    )
    parser_csv.add_argument(
        "--diseases-csv",
        type=str,
        required=True,
        help="CSV file with disease metadata (id + name).",
    )
    parser_csv.add_argument(
        "--drug-targets-csv",
        type=str,
        required=True,
        help="CSV file with drug–target associations.",
    )
    parser_csv.add_argument(
        "--gene-disease-csv",
        type=str,
        required=True,
        help="CSV file with gene–disease associations.",
    )
    parser_csv.add_argument(
        "--ppi-csv",
        type=str,
        default=None,
        help="Optional CSV file with gene–gene interactions (PPI).",
    )

    # Optional overrides for column names
    parser_csv.add_argument(
        "--drug-id-col",
        type=str,
        default="drug_id",
        help="Column name for drug IDs in drugs CSV (default: drug_id).",
    )
    parser_csv.add_argument(
        "--drug-name-col",
        type=str,
        default="drug_name",
        help="Column name for drug names in drugs CSV (default: drug_name).",
    )
    parser_csv.add_argument(
        "--gene-id-col",
        type=str,
        default="gene_id",
        help="Column name for gene IDs in genes CSV (default: gene_id).",
    )
    parser_csv.add_argument(
        "--gene-symbol-col",
        type=str,
        default="symbol",
        help="Column name for gene symbols in genes CSV (default: symbol).",
    )
    parser_csv.add_argument(
        "--disease-id-col",
        type=str,
        default="disease_id",
        help="Column name for disease IDs in diseases CSV (default: disease_id).",
    )
    parser_csv.add_argument(
        "--disease-name-col",
        type=str,
        default="disease_name",
        help="Column name for disease names in diseases CSV (default: disease_name).",
    )

    parser_csv.add_argument(
        "--dt-drug-id-col",
        type=str,
        default="drug_id",
        help="Column name for drug IDs in drug–target CSV (default: drug_id).",
    )
    parser_csv.add_argument(
        "--dt-gene-id-col",
        type=str,
        default="gene_id",
        help="Column name for gene IDs in drug–target CSV (default: gene_id).",
    )
    parser_csv.add_argument(
        "--dt-score-col",
        type=str,
        default="score",
        help="Optional column name for score in drug–target CSV (default: score).",
    )

    parser_csv.add_argument(
        "--gd-gene-id-col",
        type=str,
        default="gene_id",
        help="Column name for gene IDs in gene–disease CSV (default: gene_id).",
    )
    parser_csv.add_argument(
        "--gd-disease-id-col",
        type=str,
        default="disease_id",
        help="Column name for disease IDs in gene–disease CSV (default: disease_id).",
    )
    parser_csv.add_argument(
        "--gd-score-col",
        type=str,
        default="score",
        help="Optional column name for score in gene–disease CSV (default: score).",
    )

    parser_csv.add_argument(
        "--ppi-gene1-col",
        type=str,
        default="gene1_id",
        help="Column name for first gene ID in PPI CSV (default: gene1_id).",
    )
    parser_csv.add_argument(
        "--ppi-gene2-col",
        type=str,
        default="gene2_id",
        help="Column name for second gene ID in PPI CSV (default: gene2_id).",
    )
    parser_csv.add_argument(
        "--ppi-weight-col",
        type=str,
        default="weight",
        help="Optional column name for interaction weight in PPI CSV (default: weight).",
    )

    # Scoring & output
    parser_csv.add_argument(
        "--alpha",
        type=float,
        default=1.0,
        help="Weight for normalized overlap term in combined score (default: 1.0).",
    )
    parser_csv.add_argument(
        "--beta",
        type=float,
        default=1.0,
        help="Weight for proximity term in combined score (default: 1.0).",
    )
    parser_csv.add_argument(
        "--top-k",
        type=int,
        default=20,
        help="Number of top-ranked pairs to print (default: 20).",
    )
    parser_csv.add_argument(
        "--output-csv",
        type=str,
        default=None,
        help="If provided, save the full ranking to this CSV file.",
    )

    parser_csv.set_defaults(func=run_csv_cli)


def run_toy_cli(args: argparse.Namespace) -> None:
    data_dir = Path(args.toy_data_dir)
    if not data_dir.exists():
        raise SystemExit(f"Toy data directory does not exist: {data_dir}")

    df = run_toy_pipeline(
        data_dir=data_dir,
        alpha=args.alpha,
        beta=args.beta,
    )

    if df.empty:
        print("No drug–disease pairs found.")
        return

    _print_and_maybe_save(df, args.top_k, args.output_csv)


def run_csv_cli(args: argparse.Namespace) -> None:
    ppi_csv_path = Path(args.ppi_csv) if args.ppi_csv is not None else None

    cfg = CsvFilesConfig(
        drugs_csv=Path(args.drugs_csv),
        genes_csv=Path(args.genes_csv),
        diseases_csv=Path(args.diseases_csv),
        drug_targets_csv=Path(args.drug_targets_csv),
        gene_disease_csv=Path(args.gene_disease_csv),
        ppi_csv=ppi_csv_path,
        drug_id_col=args.drug_id_col,
        drug_name_col=args.drug_name_col,
        gene_id_col=args.gene_id_col,
        gene_symbol_col=args.gene_symbol_col,
        disease_id_col=args.disease_id_col,
        disease_name_col=args.disease_name_col,
        dt_drug_id_col=args.dt_drug_id_col,
        dt_gene_id_col=args.dt_gene_id_col,
        dt_score_col=args.dt_score_col,
        gd_gene_id_col=args.gd_gene_id_col,
        gd_disease_id_col=args.gd_disease_id_col,
        gd_score_col=args.gd_score_col,
        ppi_gene1_col=args.ppi_gene1_col,
        ppi_gene2_col=args.ppi_gene2_col,
        ppi_weight_col=args.ppi_weight_col,
    )

    df = run_csv_pipeline(
        cfg=cfg,
        alpha=args.alpha,
        beta=args.beta,
    )

    if df.empty:
        print("No drug–disease pairs found for the given CSV inputs.")
        return

    _print_and_maybe_save(df, args.top_k, args.output_csv)


def _print_and_maybe_save(
    df: pd.DataFrame,
    top_k: int,
    output_csv: str | None,
) -> None:
    top_k = min(top_k, len(df))
    print(f"\nTop {top_k} drug–disease hypotheses:\n")

    cols_to_show = [
        "drug_id",
        "drug_name",
        "disease_id",
        "disease_name",
        "n_overlap",
        "overlapping_genes",
        "mean_distance",
        "proximity_score",
        "combined_score",
    ]
    cols_to_show = [c for c in cols_to_show if c in df.columns]

    print(df[cols_to_show].head(top_k).to_string(index=False))

    if output_csv is not None:
        out_path = Path(output_csv)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        df.to_csv(out_path, index=False)
        print(f"\nFull ranking saved to: {out_path.resolve()}")


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Drug–Disease Hypothesis Generator (DDH)",
    )

    subparsers = parser.add_subparsers(
        title="subcommands",
        dest="subcommand",
        required=True,
    )

    add_toy_subcommand(subparsers)
    add_csv_subcommand(subparsers)

    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> None:
    args = parse_args(argv)
    if not hasattr(args, "func"):
        # Should not happen with required=True on subparsers,
        # but just in case.
        raise SystemExit("No subcommand specified. Use 'toy' or 'csv'.")
    args.func(args)


if __name__ == "__main__":
    main(sys.argv[1:])

