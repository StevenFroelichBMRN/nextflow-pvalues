#!/usr/bin/env python3
"""
Extract Wilcoxon rank-sum p-values and BH FDR for each perturbation in a Tahoe h5ad file.

Each h5ad contains one cell line with ~1,130 perturbations and DMSO controls.
For each perturbation, computes per-gene differential expression vs DMSO,
then filters by |log2FC| > threshold and runs Wilcoxon + BH correction.

Output: Parquet file with columns:
    cell_line, perturbation, gene, log2fc, pval, fdr, n_cells_treated, n_cells_control
"""

import argparse
import json
import sys
import gc
import numpy as np
import pandas as pd
from scipy import stats
from scipy.sparse import issparse
from statsmodels.stats.multitest import multipletests
import anndata as ad


def detect_columns(obs_df):
    """Detect perturbation and cell line columns from obs."""
    pert_col = None
    for col in ["perturbation", "drugname_drugconc", "pert_label"]:
        if col in obs_df.columns:
            pert_col = col
            break

    cl_col = None
    for col in ["cell_name", "cell_line", "celltype"]:
        if col in obs_df.columns:
            cl_col = col
            break

    return pert_col, cl_col


def process_file(filepath, log2fc_threshold, min_cells_treated, min_cells_control):
    """Process one h5ad file, return DataFrame of p-values."""

    adata = ad.read_h5ad(filepath, backed='r')
    n_cells, n_genes = adata.shape

    pert_col, cl_col = detect_columns(adata.obs)
    if pert_col is None:
        raise ValueError(f"No perturbation column found. Available: {list(adata.obs.columns)}")

    cell_line = str(adata.obs[cl_col].iloc[0]) if cl_col else "unknown"
    gene_names = list(adata.var_names)

    obs_df = adata.obs[[pert_col]].copy()
    dmso_mask = obs_df[pert_col].str.contains("DMSO", na=False).values
    n_dmso = int(dmso_mask.sum())

    if n_dmso < min_cells_control:
        adata.file.close()
        return None, cell_line, {"status": "skipped", "reason": f"too few DMSO cells ({n_dmso})"}

    dmso_indices = np.where(dmso_mask)[0]
    print(f"Loading {n_dmso} DMSO control cells...", flush=True)
    X_dmso = adata.X[dmso_indices]
    if issparse(X_dmso):
        X_dmso = X_dmso.toarray()
    X_dmso = np.asarray(X_dmso, dtype=np.float32)
    mean_dmso = X_dmso.mean(axis=0)

    non_dmso_perts = obs_df[~dmso_mask][pert_col].unique()
    print(f"Processing {len(non_dmso_perts)} perturbations...", flush=True)

    all_results = []

    for idx, pert in enumerate(non_dmso_perts):
        pert_mask = (obs_df[pert_col] == pert).values
        pert_indices = np.where(pert_mask)[0]
        n_treated = len(pert_indices)

        if n_treated < min_cells_treated:
            continue

        X_treated = adata.X[pert_indices]
        if issparse(X_treated):
            X_treated = X_treated.toarray()
        X_treated = np.asarray(X_treated, dtype=np.float32)

        mean_treated = X_treated.mean(axis=0)
        log2fc = (mean_treated - mean_dmso) / np.log(2)

        sig_indices = np.where(np.abs(log2fc) > log2fc_threshold)[0]

        if len(sig_indices) == 0:
            del X_treated
            continue

        pvals = np.ones(len(sig_indices), dtype=np.float64)
        for j, gene_idx in enumerate(sig_indices):
            treated_vals = X_treated[:, gene_idx]
            dmso_vals = X_dmso[:, gene_idx]

            if np.std(treated_vals) == 0 and np.std(dmso_vals) == 0:
                continue
            try:
                _, pval = stats.ranksums(treated_vals, dmso_vals)
                pvals[j] = pval
            except Exception:
                pass

        _, fdr_vals, _, _ = multipletests(pvals, method='fdr_bh')

        for j, gene_idx in enumerate(sig_indices):
            all_results.append({
                "cell_line": cell_line,
                "perturbation": str(pert),
                "gene": gene_names[gene_idx],
                "log2fc": float(log2fc[gene_idx]),
                "pval": float(pvals[j]),
                "fdr": float(fdr_vals[j]),
                "n_cells_treated": n_treated,
                "n_cells_control": n_dmso,
            })

        del X_treated

        if (idx + 1) % 50 == 0:
            print(f"  {idx+1}/{len(non_dmso_perts)} perturbations, {len(all_results):,} results so far", flush=True)

    adata.file.close()
    del X_dmso
    gc.collect()

    if all_results:
        df = pd.DataFrame(all_results)
        return df, cell_line, {
            "status": "success",
            "cell_line": cell_line,
            "n_perturbations": len(non_dmso_perts),
            "n_perturbations_with_results": df["perturbation"].nunique(),
            "total_results": len(df),
            "n_cells": int(n_cells),
            "n_genes": int(n_genes),
            "n_dmso": n_dmso,
        }

    return None, cell_line, {
        "status": "no_results",
        "cell_line": cell_line,
        "n_perturbations": len(non_dmso_perts),
    }


def main():
    parser = argparse.ArgumentParser(description="Extract p-values from Tahoe h5ad file")
    parser.add_argument("--input", required=True, help="Path to h5ad file")
    parser.add_argument("--output", required=True, help="Output parquet path")
    parser.add_argument("--summary", required=True, help="Output summary JSON path")
    parser.add_argument("--log2fc-threshold", type=float, default=0.25)
    parser.add_argument("--min-cells-treated", type=int, default=5)
    parser.add_argument("--min-cells-control", type=int, default=10)
    parser.add_argument("--threads", type=int, default=4)
    args = parser.parse_args()

    print(f"Processing: {args.input}", flush=True)
    print(f"Parameters: log2fc_threshold={args.log2fc_threshold}, "
          f"min_cells_treated={args.min_cells_treated}, "
          f"min_cells_control={args.min_cells_control}", flush=True)

    try:
        df, cell_line, summary = process_file(
            args.input,
            args.log2fc_threshold,
            args.min_cells_treated,
            args.min_cells_control,
        )

        if df is not None and len(df) > 0:
            df.to_parquet(args.output, index=False)
            print(f"Wrote {len(df):,} rows to {args.output}", flush=True)
        else:
            pd.DataFrame(columns=[
                "cell_line", "perturbation", "gene", "log2fc", "pval", "fdr",
                "n_cells_treated", "n_cells_control"
            ]).to_parquet(args.output, index=False)
            print(f"No results for {cell_line}", flush=True)

        summary["file"] = args.input
        with open(args.summary, "w") as f:
            json.dump(summary, f, indent=2)

    except Exception as e:
        print(f"ERROR: {type(e).__name__}: {e}", file=sys.stderr, flush=True)
        summary = {"status": "error", "file": args.input, "error": str(e)[:500]}
        with open(args.summary, "w") as f:
            json.dump(summary, f, indent=2)
        pd.DataFrame(columns=[
            "cell_line", "perturbation", "gene", "log2fc", "pval", "fdr",
            "n_cells_treated", "n_cells_control"
        ]).to_parquet(args.output, index=False)
        sys.exit(1)


if __name__ == "__main__":
    main()
