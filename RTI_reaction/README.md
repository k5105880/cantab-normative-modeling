# Reaction time: normative modeling

Uses a Student-t distribution for log(1 + milliseconds). Participant-specific intercepts can differ by platform, with a common participant time slope.

The measure is `RTI_Median_five-choice_reaction_time` in the input CSV (ms). Higher scores mean poorer performance.

## Prepare and run

Read the [input guide](../docs/DATA.md), complete [environment setup](../environment/README.md), and use the [execution/settings/result guide](../docs/USAGE.md).

First create the completed input CSV using the input guide. From `project/working_directory/` (the repository root):

```r
source("run_analysis.R")
input <- "../input/CANTAB_connect_RS6_rawscores.csv"
run_cantab("RTI_reaction", input, "../output/RTI_reaction_check", mode = "check")
run_cantab("RTI_reaction", input, "../output/RTI_reaction_fit", mode = "fresh")
# Only if fitting was interrupted, with unchanged inputs, code and settings:
run_cantab("RTI_reaction", input, "../output/RTI_reaction_fit", mode = "resume")
```

`check` validates input/model design without sampling. `fresh` fits the normative model and calculates scores in a new directory. `resume` continues an interrupted fitting run.

## Results

Open **`normative_scores.csv`** in the specified output directory. Its **`normative_z`** column is the normative z score for each assessment. `raw_score` retains the input value; `score_status` identifies rows with a usable score. `score_dictionary.csv` explains all columns, and `analysis_summary.txt` gives the model and score direction.

Eligible healthy controls are scored with participant-held-out cross-validation. Other participants use the model fitted to all eligible healthy controls. All selected analysis rows are preserved, including those without a score. With the default `reference_dir = NULL`, these are all original input rows; an optional cohort selection excludes unlisted participants before scoring. Use the [diagnostic guide](../docs/DIAGNOSTICS.md) for exact log/result paths and interpretation of the model checks.

## Default model settings

- Model: Student-t model for log-transformed reaction time with platform-specific participant intercepts.
- Four chains; all-control fit: 5000 total iterations including 1500 warmup per chain.
- Five-fold cross-validation; each fit: 3000 total iterations including 1000 warmup per chain.
- Sampler target acceptance (`adapt_delta`): 0.999; maximum tree depth: 20.
- A failed convergence check triggers one recovery attempt with twice the iterations and warmup, target acceptance 0.99999, and the next seed. Scores require successful convergence, numerical integration, and posterior precision checks.

The mean includes nonlinear baseline age and elapsed time, gender, baseline IQ, practice, platform, and platform interactions with time and practice. The time smooth excludes its linear term so that platform-specific linear slopes are identifiable. Residual scale depends on platform interactions with practice and baseline age. Participant effects are correlated platform-specific intercepts and a time slope.

Normative percentiles integrate over the effects of a new participant. They use all ordered posterior draws and preserve log-scale tail probabilities. Numerical integration checks the z score and log density at increasing quadrature orders; the posterior Monte Carlo standard error of z must not exceed 0.03.

[Model metadata](provenance/model.json) describes the defaults; it is documentation, not a user configuration file. The [model settings guide](../docs/SETTINGS.md) explains the formulas, priors, data-dependent rules and recovery attempts.

Advanced implementation references:

- [rti_core.R](rti_core.R)
- [rti_cdf.R](rti_cdf.R)

Alternatively, source this folder's `run.R` and call `run_outcome()` with the same arguments except `outcome`. Keep the repository structure for shared input and result helpers.

From the repository root in R, `source("scripts/run_tests.R")` runs the public software tests, including reaction-time numerical integration and precision, without production model fitting. No saved fit or private input is required. The optional `RTI_TEST_FIT` environment variable is only for advanced saved-fit adapter checks; leave it unset for the standard tests.
