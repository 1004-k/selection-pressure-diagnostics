# Precomputed full-run outputs

This folder contains the complete output archive from the final v7 run:

- `output_jci_final_for_review_v7_timeupdated.zip`
- SHA256: `7e0431c2e89c716f6aba8a78d861d1b5c7628ca1b1a710877168569ace778cb9`

The archive includes raw replicate-level simulation outputs, performance summaries, figures, tables, logs, and QA reports. It is provided so reviewers can inspect the exact outputs used to create the manuscript without rerunning the full simulation.

To regenerate these outputs from code, run from the repository root:

```bash
Rscript run_jci_final.R full
Rscript run_jci_timeupdated_only.R full
```

The full run used `B=200`, `N=2000`, `N_TRUTH=200000`, `B_TRUTH=5`, and `N_CORES=12` on a Debian/GCP VM. See `results/audit/sessionInfo_rescue.txt` and `results/logs/` for details.
