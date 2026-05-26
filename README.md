# Selection-pressure diagnostics for per-protocol causal inference

This repository contains the reproducibility code for the Journal of Causal Inference submission:

**Selection-pressure diagnostics for per-protocol causal inference**

The project evaluates how selection pressure, measured by a time-resolved prognosis-dependent deviation index `SPD(t)`, relates to the operating characteristics of per-protocol estimators. The JCI version uses several safeguards added after code and estimand review:

1. **Aligned estimand for IPCW and DR-AIPW.** Standardized IPCW and DR-AIPW both target the same marginal fixed-horizon per-protocol risk contrast.
2. **Intervention-based Monte Carlo truth.** The non-null truth is computed by direct intervention on treatment assignment under the known data-generating mechanism, rather than by applying an evaluated estimator to a large simulated dataset.
3. **Risk-scale primary estimands.** The main analysis reports fixed-horizon log risk ratio and risk difference. Cox hazard-ratio outputs are kept as a supplementary consistency check only.
4. **Time-updated IPCW sensitivity analysis.** A supplementary IPCW estimator fits the deviation/adherence model on person-period data with observed time-updated prognosis and flexible time terms.

## Repository structure

```text
R/                         Core simulation and estimator functions
scripts/                   Main run, table, figure, and QA scripts
run_jci_final.R            R-only runner for smoke/full main reproduction
run_jci_timeupdated_only.R R-only runner for time-updated sensitivity rerun
docs/                      Method notes and guardrails
results/figures/           Final manuscript figures
results/tables/            Final main and supplementary tables
results/audit/             QA reports and session information
results/logs/              Full-run logs
precomputed/               Complete final output archive and SHA256 checksum
manuscript/                Manuscript and supplement source/PDF snapshot
install_packages.R         Minimal package installer
```

## Requirements

The final run was executed with R 4.2.2 on Debian GNU/Linux. Required CRAN packages:

- `data.table`
- `ggplot2`
- `survival`

Base R packages such as `stats`, `utils`, `parallel`, and `splines` are also used.

Install dependencies:

```bash
Rscript install_packages.R
```

## Quick smoke test

From the repository root:

```bash
Rscript run_jci_final.R smoke
```

This writes a small test run to `output_jci_smoke_final/`. Smoke-test outputs are **not** manuscript results.

## Full reproduction

From the repository root:

```bash
export N_CORES=12   # adjust to your machine
Rscript run_jci_final.R full
```

This generates:

```text
output_jci_final/raw/
output_jci_final/tables/
output_jci_final/figures/
output_jci_final/audit/
```

The full final run used:

```text
B=200
N=2000
TRUTH_LEVELS=null,non_null
MIS_SPEC_LEVELS=0,1
N_TRUTH=200000
B_TRUTH=5
RUN_DR=1
RUN_TIME_UPDATED_SENS=1
IPCW_ESTIMATOR=hajek
```

## Time-updated sensitivity rerun only

After a full main run exists in `output_jci_final/`, the time-updated IPCW sensitivity can be rerun alone:

```bash
export N_CORES=12
Rscript run_jci_timeupdated_only.R full
```

The final v7 time-updated sensitivity QA showed that the actual matrix-based time-updated GLM was used in all rows, with no fallback rows.

## Precomputed outputs

The exact final output archive used for the submission package is provided in:

```text
precomputed/output_jci_final_for_review_v7_timeupdated.zip
```

Verify integrity:

```bash
sha256sum -c precomputed/output_jci_final_for_review_v7_timeupdated.zip.sha256
```

## QA status of the final run

The final v7 outputs passed row-count and nonfinite-value checks. The time-updated sensitivity patch also passed model-status checks:

```text
proportion glm_timeupdated status: 1.0000
fallback rows: 0
```

See:

```text
results/audit/JCI_v7_full_output_QA_report.md
results/audit/jci_output_qa_report.txt
results/audit/timeupdated_sensitivity_patch_qa.txt
```

## Notes for public GitHub release

Before making the repository public, update:

1. `CITATION.cff`, especially `repository-code`.
2. Any author/contributor information you want visible.
3. The manuscript citation once the article has a DOI.

Do not mix older EJE/Paper B output files with this repository. The final JCI manuscript should use only the v7 outputs documented here.
