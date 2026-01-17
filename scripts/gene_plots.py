import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

OUTDIR = Path("artifacts/stage2_deliverable")
OUTDIR.mkdir(parents=True, exist_ok=True)

# -------------------------
# Plot 1: Top genes barplot
# -------------------------
coefs_path = "artifacts/model_runs/cox_cv/2026-01-16_195531/tables/final_model_coefs_named.csv"
df = pd.read_csv(coefs_path)

df["abs_coef"] = df["coef"].abs()
top = df.sort_values("abs_coef", ascending=False).head(20).copy()
colors = ["red" if c > 0 else "blue" for c in top["coef"]]

plt.figure(figsize=(6, 8))
plt.barh(top["gene_name"], top["coef"], color=colors)
plt.axvline(0, color="black", linewidth=0.8)
plt.xlabel("Cox Coefficient (log hazard ratio)")
plt.title("Top Survival-Associated Genes (LUAD)")
plt.gca().invert_yaxis()
plt.tight_layout()
plt.savefig(OUTDIR / "top_genes_barplot.png", dpi=200)
plt.close()

# -------------------------
# Plot 2: Risk pathway dotplot
# -------------------------
ora_path = "artifacts/model_runs/cox_cv/2026-01-16_195531/enrichment/hallmark_risk_ora_bg/ora_results.csv"
ora = pd.read_csv(ora_path)

ora = ora.sort_values("fdr_bh").head(10).copy()
# avoid -log10(0) just in case
ora["neglogFDR"] = -np.log10(np.clip(ora["fdr_bh"], 1e-300, 1.0))

plt.figure(figsize=(8, 5))
plt.scatter(
    ora["neglogFDR"],
    ora["pathway"],
    s=ora["overlap"] * 20,
    color="crimson"
)
plt.xlabel("-log10(FDR)")
plt.title("Enriched Hallmark Pathways (Risk Genes)")
plt.tight_layout()
plt.savefig(OUTDIR / "risk_pathway_dotplot.png", dpi=200)
plt.close()

print(f"[DONE] Wrote plots to: {OUTDIR}")