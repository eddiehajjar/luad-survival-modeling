from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
from typing import Dict, List, Set, Tuple

import pandas as pd


def read_gmt(path: str) -> Dict[str, Set[str]]:
    """
    Read a GMT file into {pathway_name -> set(gene_symbols)}.
    GMT format: name <tab> description/url <tab> gene1 <tab> gene2 ...
    """
    gene_sets: Dict[str, Set[str]] = {}
    with open(path, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < 3:
                continue
            name = parts[0]
            genes = set(g for g in parts[2:] if g)
            gene_sets[name] = genes
    return gene_sets


def bh_fdr(pvals: List[float]) -> List[float]:
    """Benjamini-Hochberg FDR correction."""
    n = len(pvals)
    idx = sorted(range(n), key=lambda i: pvals[i])
    q = [0.0] * n
    prev = 1.0
    for rank, i in enumerate(reversed(idx), start=1):
        # reversed ranks from n..1 to enforce monotonicity
        r = n - rank + 1
        val = pvals[i] * n / r
        prev = min(prev, val)
        q[i] = min(prev, 1.0)
    return q


def log_choose(n: int, k: int) -> float:
    # log(n choose k) using lgamma for stability
    if k < 0 or k > n:
        return float("-inf")
    return math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)


def hypergeom_sf(x: int, N: int, K: int, n: int) -> float:
    """
    Survival function P[X >= x] for hypergeometric.
    X ~ Hypergeom(N population, K successes in population, n draws)
    """
    # Sum_{i=x..min(K,n)} [C(K,i) C(N-K, n-i)] / C(N,n)
    max_i = min(K, n)
    if x > max_i:
        return 0.0
    denom = log_choose(N, n)
    # log-sum-exp
    logs = []
    for i in range(x, max_i + 1):
        logs.append(log_choose(K, i) + log_choose(N - K, n - i) - denom)
    m = max(logs)
    return float(sum(math.exp(li - m) for li in logs) * math.exp(m))


def ora_enrichment(
    gene_list: Set[str],
    background: Set[str],
    gene_sets: Dict[str, Set[str]],
    min_set_size: int = 10,
    max_set_size: int = 5000,
) -> pd.DataFrame:
    """
    Over-representation analysis (hypergeometric test) for gene_list against each gene set.

    background = universe genes you consider "measurable" in your experiment.
    """
    bg = set(background)
    genes = set(gene_list) & bg

    N = len(bg)
    n = len(genes)

    rows = []
    for pathway, members in gene_sets.items():
        members_bg = set(members) & bg
        K = len(members_bg)
        if K < min_set_size or K > max_set_size:
            continue

        overlap = genes & members_bg
        x = len(overlap)
        if x == 0:
            p = 1.0
        else:
            p = hypergeom_sf(x, N=N, K=K, n=n)

        rows.append(
            {
                "pathway": pathway,
                "set_size": K,
                "list_size": n,
                "overlap": x,
                "overlap_genes": ";".join(sorted(overlap)),
                "p_value": p,
            }
        )

    df = pd.DataFrame(rows)
    if len(df) == 0:
        return df

    df["fdr_bh"] = bh_fdr(df["p_value"].tolist())
    df = df.sort_values(["fdr_bh", "p_value", "overlap"], ascending=[True, True, False]).reset_index(drop=True)
    return df


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--ranked", required=True, help="CSV containing at least gene_name column (and coef optional).")
    ap.add_argument("--gmt", required=True, help="GMT file (Hallmark).")
    ap.add_argument("--outdir", required=True, help="Output directory.")
    ap.add_argument("--gene-col", default="gene_name", help="Gene symbol column in ranked CSV.")
    ap.add_argument("--background", default=None, help="Optional CSV/TSV file with a gene_name column for universe.")
    ap.add_argument("--min-set-size", type=int, default=10)
    ap.add_argument("--max-set-size", type=int, default=5000)
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    ranked_df = pd.read_csv(args.ranked)
    if args.gene_col not in ranked_df.columns:
        raise ValueError(f"Missing gene column '{args.gene_col}' in {args.ranked}. Columns: {list(ranked_df.columns)}")

    gene_list = set(ranked_df[args.gene_col].dropna().astype(str).str.strip().tolist())

    # --- background universe (optional) ---
    background_genes = None
    if args.background:
        bg = pd.read_csv(args.background)

        # If it didn't load a header (e.g., columns are [0]), retry as 1-col file.
        if args.gene_col not in bg.columns:
            bg = pd.read_csv(args.background, header=None, names=[args.gene_col])

        if args.gene_col not in bg.columns:
            raise ValueError(
                f"Missing gene column '{args.gene_col}' in background file. Columns: {list(bg.columns)}"
            )

        background_genes = set(bg[args.gene_col].dropna().astype(str).str.strip().tolist())

    gene_sets = read_gmt(args.gmt)

    results = ora_enrichment(
        gene_list=gene_list,
        background=background_genes,
        gene_sets=gene_sets,
        min_set_size=args.min_set_size,
        max_set_size=args.max_set_size,
    )
    # Save outputs
    results_csv = outdir / "ora_results.csv"
    results.to_csv(results_csv, index=False)

    meta = {
        "ranked": args.ranked,
        "gmt": args.gmt,
        "outdir": str(outdir),
        "gene_col": args.gene_col,
        "min_set_size": args.min_set_size,
        "max_set_size": args.max_set_size,
        "gene_list_size": len(gene_list),
        "background_size": len(background_genes),
        "n_tested_sets": int(len(results)),
    }
    (outdir / "run_meta.json").write_text(json.dumps(meta, indent=2))

    # Print top hits
    if len(results) == 0:
        print("[WARN] No results (check inputs / set size filters).")
        return

    print("\nTop pathways by FDR:")
    cols = ["pathway", "overlap", "set_size", "p_value", "fdr_bh"]
    print(results[cols].head(15).to_string(index=False))


if __name__ == "__main__":
    main()
