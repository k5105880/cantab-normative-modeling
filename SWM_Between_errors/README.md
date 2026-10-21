# Spatial working-memory errors: normative modeling

Separates zero errors from positive counts and uses a discretized normal distribution for positive counts. Relative trial exposure is fixed at RS6 = 4 and Connect = 1. Confirm that this 4:1 ratio matches your test protocol. Input counts are not divided by four.

The measure is `SWM_Between_errors` in the input CSV (count). Higher scores mean poorer performance.

## Prepare and run

Read the [input guide](../docs/DATA.md), complete [environment setup](../environment/README.md), and use the [execution/settings/result guide](../docs/USAGE.md).

First create the completed input CSV using the input guide. From `project/working_directory/` (the repository root):

```r
source("run_analysis.R")
input <- "../input/CANTAB_connect_RS6_rawscores.csv"
run_cantab("SWM_Between_errors", input, "../output/SWM_Between_errors_check", mode = "check")
run_cantab("SWM_Between_errors", input, "../output/SWM_Between_errors_fit", mode = "fresh")
# Only if fitting was interrupted, with unchanged inputs, code and settings:
run_cantab("SWM_Between_errors", input, "../output/SWM_Between_errors_fit", mode = "resume")
```

`check` validates input/model design without sampling. `fresh` fits the normative model and calculates scores in a new directory. `resume` continues an interrupted fitting run.

## Results

Open **`normative_scores.csv`** in the specified output directory. Its **`normative_z`** column is the normative z score for each assessment. `raw_score` retains the input value; `score_status` identifies rows with a usable score. `score_dictionary.csv` explains all columns, and `analysis_summary.txt` gives the model and score direction.

Eligible healthy controls are scored with participant-held-out cross-validation. Other participants use the model fitted to all eligible healthy controls. All selected analysis rows are preserved, including those without a score. With the default `reference_dir = NULL`, these are all original input rows; an optional cohort selection excludes unlisted participants before scoring. Use the [diagnostic guide](../docs/DIAGNOSTICS.md) for exact log/result paths and interpretation of the model checks.

## Default model settings

- Model: Hurdle discretized normal model for error counts with trial-count adjustment.
- Four chains; all-control fit: 3000 total iterations including 1000 warmup per chain.
- Five-fold cross-validation; each fit: 3000 total iterations including 1000 warmup per chain.
- Sampler target acceptance (`adapt_delta`): 0.995; maximum tree depth: 18.

[Model metadata](provenance/model.json) describes the defaults; it is documentation, not a user configuration file. The [model settings guide](../docs/SETTINGS.md) explains the formulas, priors, data-dependent rules and recovery attempts.

## Source files

The workflow fits one normative model. Its definitions are organized by role:

| File | Role |
|---|---|
| [vendor/common_model_definitions.R](vendor/common_model_definitions.R) | Mean formula and shared priors |
| [vendor/final_model_definition.R](vendor/final_model_definition.R) | Complete normative model, including scale and zero-error probability |
| [vendor/count_distribution.R](vendor/count_distribution.R) | Count probabilities, likelihood and predictive random-number generation |
| [R/models.R](R/models.R) | Fitting, diagnostics and prediction dispatch |
| [vendor/probability_scores.R](vendor/probability_scores.R) | Conversion of probabilities to z scores and calibration moments |
| [vendor/normative_prediction.R](vendor/normative_prediction.R) | Normative predictions and predictive-interval coverage |

The model uses `normative` as its internal name. The main score column is `normative_z` in both the intermediate score table and the common output.

The files separate shared model definitions, the fitted model specification, probability calculations, fitting and prediction, and score conversion. The workflow fits one model specification.

Alternatively, source this folder's `run.R` and call `run_outcome()` with the same arguments except `outcome`. Keep the repository structure for shared input and result helpers.
