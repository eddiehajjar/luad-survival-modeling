Interpretable Survival Modeling in Lung Adenocarcinoma (TCGA-LUAD)

This project implements an end-to-end survival analysis pipeline using tumor RNA-seq data from TCGA lung adenocarcinoma (LUAD). The focus is on correct time-to-event modeling, stability in high-dimensional data, and biological interpretability, rather than black-box prediction.

⸻

Problem

Most cancer ML models reduce outcomes to binary labels and ignore censoring, time-dependence, and interpretability. In real oncology settings, this leads to brittle models that are difficult to trust or reason about biologically.

The goal here was not just to predict survival, but to:
	•	model time-to-event correctly
	•	handle ~60k genes with ~500 patients
	•	extract interpretable, biologically meaningful signal

⸻

Approach
	•	Data
	•	TCGA-LUAD tumor RNA-seq and curated overall survival endpoints
	•	Explicit censoring logic; no proxy labels
	•	Modeling
	•	Penalized Cox proportional hazards models
	•	Cross-validation to control overfitting
	•	L2 regularization for coefficient stability
	•	Top-variance gene selection performed within each fold to prevent leakage
	•	Validation
	•	Continuous risk scores derived from gene expression
	•	Kaplan–Meier stratification of high- vs low-risk patients
	•	Log-rank tests for time-to-event separation
	•	Interpretation
	•	Extraction of risk-increasing and protective genes from the multivariate model
	•	Background-aware pathway enrichment using MSigDB Hallmark gene sets
	•	Validation that learned signal reflects coherent oncogenic programs

⸻

Results
	•	Cross-validated C-index ≈ 0.63
	•	Clear survival separation between risk groups
	•	Enrichment of known aggressive LUAD biology, including:
	•	Epithelial–mesenchymal transition (EMT)
	•	Hypoxia and glycolysis
	•	KRAS signaling
	•	TNFα / NF-κB inflammatory signaling
	•	p53 pathway disruption
