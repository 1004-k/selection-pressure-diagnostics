# JCI v7 full output QA report

## Overall verdict
The v7 output package is usable for the JCI manuscript. The main standardized IPCW and DR-AIPW full run is complete and internally consistent. The patched time-updated IPCW sensitivity rerun also succeeded: all 14,400 expected rows were generated and 100% of rows used the matrix-based time-updated GLM rather than the previous fallback model.

## Run integrity
- `replicate_results_rescue.csv`: 14,400 rows, expected 14,400.
- `ipcw_risk_rescue.csv`: 14,400 rows, expected 14,400.
- `dr_risk_rescue.csv`: 14,400 rows, expected 14,400.
- `ipcw_timeupdated_risk_rescue.csv`: 14,400 rows, expected 14,400.
- `spd_curves_rescue.csv`: 72,000 rows, expected 72,000.
- `weight_diagnostics_rescue.csv`: 288,000 rows, expected 288,000.
- `mc_truth_risk.csv`: 36 rows, expected 36.
- Nonfinite core estimates: 0 across main IPCW, DR-AIPW, SPD, weight diagnostics, MC truth, and time-updated sensitivity outputs.

## Time-updated sensitivity status
The previous v5/v6 problem is fixed.

- `proportion glm_timeupdated status`: 1.0000.
- `glm_timeupdated_richer`: 7,200 rows.
- `glm_timeupdated_reduced`: 7,200 rows.
- No fallback rows remain.

The time-updated IPCW estimator should be presented as a supplementary sensitivity analysis, not as a new primary estimator.

## Main operating characteristics, median across scenarios

| Nuisance regime | Type I IPCW | Coverage IPCW | Power IPCW | Sign error IPCW | Type I DR | Coverage DR | Power DR | Sign error DR |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Richer baseline set | 0.0425 | 0.9575 | 0.3775 | 0.0500 | 0.0425 | 0.9575 | 0.3875 | 0.0525 |
| Reduced baseline set | 0.0850 | 0.9150 | 0.2300 | 0.1225 | 0.0875 | 0.9125 | 0.2400 | 0.1250 |

Interpretation: after standardized IPCW correction, standardized IPCW and DR-AIPW behave very similarly. The manuscript should not claim that DR-AIPW broadly rescues IPCW. The safer claim is that SPD(t) indexes regimes in which both estimators remain stable versus regimes in which both degrade under reduced nuisance information.

## Time-updated IPCW sensitivity, median across scenarios

| Nuisance regime | Type I time-updated IPCW | Coverage time-updated IPCW | Power time-updated IPCW | Sign error time-updated IPCW | Median min c-hat |
|---:|---:|---:|---:|---:|---:|
| Richer baseline set | 0.0800 | 0.9200 | 0.1825 | 0.1325 | 0.0851 |
| Reduced baseline set | 0.0975 | 0.9025 | 0.1800 | 0.1425 | 0.1237 |

Interpretation: the time-updated sensitivity did not improve performance in this simulation setting. It was more conservative/less powerful and showed higher sign error than the primary standardized baseline/horizon-level IPCW. This is not a fatal issue if it is framed correctly: it shows that simply adding a time-updated prognostic score to the adherence model does not automatically rescue per-protocol inference under noisy selection pressure.

## Figures
- Figure 1 is technically acceptable: 2055 x 1470 pixels at 300 dpi. It is readable and should be retained.
- Figure 2 is technically acceptable but visually weak in its current labelled version because several scenario labels overlap. A cleaner unlabelled version has been generated:
  - `Figure2_operating_map_JCI_clean.png`
  - `Figure2_operating_map_JCI_clean.pdf`

Recommendation: use the clean unlabelled Figure 2 in the main manuscript and identify scenario IDs through Table 3 or supplementary source data.

## Remaining cautions before manuscript submission
1. Do not use the old v5 time-updated outputs. Use only v7 outputs.
2. Do not say “DR-AIPW rescued IPCW.” The corrected results show near-equivalence between standardized IPCW and DR-AIPW.
3. Use “richer baseline nuisance set” and “reduced baseline nuisance set,” not “well-specified” versus “misspecified,” unless the code implements the exact true time-varying nuisance model.
4. Mention that the time-updated sensitivity was supplementary and that it did not materially improve operating characteristics.
5. If publishing CSV outputs, consider renaming the `truth` value `null` to `causal_null` because some software imports the string `null` as missing by default.
6. Exclude `t=0` from any displayed weight-tail diagnostic summaries if possible, because baseline tail-share can be mechanically uninformative.

## Recommended manuscript framing
The main contribution should be framed as:

> We evaluated standardized per-protocol estimators targeting the same marginal fixed-horizon risk contrast. After aligning the target parameter for IPCW and DR-AIPW, their operating characteristics were similar under a richer baseline nuisance set. Under reduced nuisance information, both estimators showed inflated Type I error, reduced power, poorer coverage, and increased sign error. SPD(t) therefore functions as a diagnostic severity axis for interpreting estimator stability and for identifying settings in which neither estimator should be treated as decision-grade without stronger design or modelling information.

