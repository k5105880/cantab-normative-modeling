# Spatial working-memory strategy: normative modeling

Combines additional probability at the minimum with a discretized beta distribution. RS6 scores are divided by four inside the model; the different score grids are retained.

The measure is `SWM_Strategy` in the input CSV (score). Higher scores mean poorer performance.

## Prepare and run

Read the [input guide](../docs/DATA.md), complete [environment setup](../environment/README.md), and use the [execution/settings/result guide](../docs/USAGE.md).

First create the completed input CSV using the input guide. From `project/working_directory/` (the repository root):

```r
source("run_analysis.R")
input <- "../input/CANTAB_connect_RS6_rawscores.csv"
run_cantab("SWM_Strategy", input, "../output/SWM_Strategy_check", mode = "check")
run_cantab("SWM_Strategy", input, "../output/SWM_Strategy_fit", mode = "fresh")
# Only if fitting was interrupted, with unchanged inputs, code and settings:
run_cantab("SWM_Strategy", input, "../output/SWM_Strategy_fit", mode = "resume")
```

`check` validates input/model design without sampling. `fresh` fits the normative model and calculates scores in a new directory. `resume` continues an interrupted fitting run.

## Results

Open **`normative_scores.csv`** in the specified output directory. Its **`normative_z`** column is the normative z score for each assessment. `raw_score` retains the input value; `score_status` identifies rows with a usable score. `score_dictionary.csv` explains all columns, and `analysis_summary.txt` gives the model and score direction.

Eligible healthy controls are scored with participant-held-out cross-validation. Other participants use the model fitted to all eligible healthy controls. All selected analysis rows are preserved, including those without a score. With the default `reference_dir = NULL`, these are all original input rows; an optional cohort selection excludes unlisted participants before scoring. Use the [diagnostic guide](../docs/DIAGNOSTICS.md) for exact log/result paths and interpretation of the model checks.

The pipeline contains one normative model. It exports assessment-level scores and calibration diagnostics.

## Default model settings

- Model: Minimum-inflated discretized beta model for strategy scores.
- Four chains; all-control fit: 4000 total iterations including 1500 warmup per chain.
- Five-fold cross-validation; each fit: 4000 total iterations including 1500 warmup per chain.
- Sampler target acceptance (`adapt_delta`): 0.999; maximum tree depth: 12.

[Model metadata](provenance/model.json) describes the defaults; it is documentation, not a user configuration file. The [model settings guide](../docs/SETTINGS.md) explains the formulas, priors, data-dependent rules and recovery attempts.

Advanced implementation references:

- [R/swm_core.R](R/swm_core.R)
- [R/swm_discrete_family.R](R/swm_discrete_family.R)
- [stan/swm_discrete_functions.stan](stan/swm_discrete_functions.stan)

Alternatively, source this folder's `run.R` and call `run_outcome()` with the same arguments except `outcome`. Keep the repository structure for shared input and result helpers.
