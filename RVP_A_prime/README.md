# Attention sensitivity (A prime): normative modeling

Uses an extended-support beta distribution with observations censored to the interval from 0 to 1. Exact values of 0 and 1 are supported.

The measure is `RVP_A_prime` in the input CSV (0–1). Higher scores mean better performance.

## Prepare and run

Read the [input guide](../docs/DATA.md) ([日本語](../docs/DATA_JA.md)), complete [environment setup](../environment/README.md), and use the [execution/settings/result guide](../docs/USAGE.md) ([日本語](../docs/USAGE_JA.md)).

First create the completed input CSV using the input guide. From `project/working_directory/` (the repository root):

```r
source("run_analysis.R")
input <- "../input/CANTAB_connect_RS6_rawscores.csv"
run_cantab("RVP_A_prime", input, "../output/RVP_A_prime_check", mode = "check")
run_cantab("RVP_A_prime", input, "../output/RVP_A_prime_fit", mode = "fresh")
# Only if fitting was interrupted, with unchanged inputs, code and settings:
run_cantab("RVP_A_prime", input, "../output/RVP_A_prime_fit", mode = "resume")
```

`check` validates input/model design without sampling. `fresh` fits the models and calculates scores in a new directory. `resume` continues an interrupted fitting run.

## Results

Open **`normative_scores.csv`** in the specified output directory. Its **`normative_z`** column is the normative z score for each assessment. `raw_score` retains the input value; `score_status` identifies rows with a usable score. `score_dictionary.csv` explains all columns, and `analysis_summary.txt` gives the model and score direction.

Eligible healthy controls are scored with participant-held-out cross-validation. Other participants use the model fitted to all eligible healthy controls. All selected analysis rows are preserved, including those without a score. With the default `reference_dir = NULL`, these are all original input rows; an optional cohort selection excludes unlisted participants before scoring. Use the [diagnostic guide](../docs/DIAGNOSTICS.md) ([日本語](../docs/DIAGNOSTICS_JA.md)) for exact log/result paths and interpretation of the model checks.

## Default model settings

- Model: Extended-support beta model for attention sensitivity.
- Four chains; all-control fit: 5000 total iterations including 1500 warmup per chain.
- Five-fold cross-validation; each fit: 5000 total iterations including 1500 warmup per chain.
- Sampler target acceptance (`adapt_delta`): 0.9995; maximum tree depth: 20.

[Model metadata](provenance/model.json) describes the defaults; it is documentation, not a user configuration file. The [model settings guide](../docs/SETTINGS.md) ([日本語](../docs/SETTINGS_JA.md)) explains the formulas, priors, data-dependent rules and recovery attempts.

Numerical scoring integrates the new-participant random effects using 32, 64, 128 and 256 quadrature nodes, stopping when successive z scores differ by at most 0.002. Scores must also have posterior Monte Carlo standard error at most 0.03. The direct prior check uses 40,000 independent parameter draws; calibration uses 2,000 participant-cluster bootstrap replicates.

This folder implements one normative model workflow, fitted to all eligible healthy controls and to five participant-held-out training sets. Each training set must contain at least 20 participants with two or more eligible assessments. In each training set, age/time use smooth terms when there are at least 6/5 distinct values, otherwise linear terms. FSIQ is included when it varies; gender, platform and occasion terms depend on represented levels. See the model settings guide for the exact rules. The model uses a fixed seed and four sampling chains. An initial fit rejected solely for divergent transitions may receive one retry with 2,500 warmup iterations and target acceptance 0.9999, preserving 3,500 retained iterations per chain. Every accepted fit is checked for complete draws, convergence, effective sample size, tree depth and energy diagnostics.

Implementation reference:

- [R/model_spec.R](R/model_spec.R)

Alternatively, source this folder's `run.R` and call `run_outcome()` with the same arguments except `outcome`. Keep the repository structure for shared input and result helpers.
