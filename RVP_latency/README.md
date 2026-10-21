# Response latency: normative modeling

Uses a Gamma distribution directly on positive milliseconds, with a log link for the mean and platform/occasion effects on its shape.

The measure is `RVP_Median_latency` in the input CSV (ms). Higher scores mean poorer performance.

## Prepare and run

Read the [input guide](../docs/DATA.md), complete [environment setup](../environment/README.md), and use the [execution/settings/result guide](../docs/USAGE.md).

First create the completed input CSV using the input guide. From `project/working_directory/` (the repository root):

```r
source("run_analysis.R")
input <- "../input/CANTAB_connect_RS6_rawscores.csv"
run_cantab("RVP_latency", input, "../output/RVP_latency_check", mode = "check")
run_cantab("RVP_latency", input, "../output/RVP_latency_fit", mode = "fresh")
# Only if fitting was interrupted, with unchanged inputs, code and settings:
run_cantab("RVP_latency", input, "../output/RVP_latency_fit", mode = "resume")
```

`check` validates input/model design without sampling. `fresh` fits the models and calculates scores in a new directory. `resume` continues an interrupted fitting run.

## Results

Open **`normative_scores.csv`** in the specified output directory. Its **`normative_z`** column is the normative z score for each assessment. `raw_score` retains the input value; `score_status` identifies rows with a usable score. `score_dictionary.csv` explains all columns, and `analysis_summary.txt` gives the model and score direction.

Eligible healthy controls are scored with participant-held-out cross-validation. Other participants use the model fitted to all eligible healthy controls. All selected analysis rows are preserved, including those without a score. With the default `reference_dir = NULL`, these are all original input rows; an optional cohort selection excludes unlisted participants before scoring. Use the [diagnostic guide](../docs/DIAGNOSTICS.md) for exact log/result paths and interpretation of the model checks.

## Default model settings

- Model: Gamma model for response latency.
- Four chains; all-control fit: 5000 total iterations including 1500 warmup per chain.
- Five-fold cross-validation; each fit: 3000 total iterations including 1000 warmup per chain.
- Sampler target acceptance (`adapt_delta`): 0.999; maximum tree depth: 20.

[Model metadata](provenance/model.json) describes the defaults; it is documentation, not a user configuration file. The [model settings guide](../docs/SETTINGS.md) explains the formulas, priors, data-dependent rules and recovery attempts.

Advanced implementation references:

- [model_functions.R](model_functions.R) and [model_settings.R](model_settings.R)
- [identifiability_runtime.R](identifiability_runtime.R)

Alternatively, source this folder's `run.R` and call `run_outcome()` with the same arguments except `outcome`. Keep the repository structure for shared input and result helpers.

The normative workflow fits the full eligible healthy-control sample and five participant-held-out training sets. It uses a participant random slope when at least 20 participants have repeated eligible assessments in the full training data; otherwise it uses a random intercept. A non-backend/non-compilation error during the full random-slope fit can also trigger the defined random-intercept retry. A failed convergence gate is not that retry trigger. Cross-validation refits the resulting full-model structure; see the settings guide for details. New-participant random effects are integrated numerically; scores are withheld if integration or posterior Monte Carlo precision fails.
