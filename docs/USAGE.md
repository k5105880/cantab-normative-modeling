# Running an analysis and reading the results

## 1. Prepare folders and software

The examples assume the following layout, with `working_directory/` as the R working directory.

```text
project/
  working_directory/    code; set this as the R working directory
  input/                participant and assessment tables, prepared input CSV
  output/               result folders for individual runs
```

Follow [R/Stan setup](../environment/README.md) once, then [prepare the CSV](DATA.md). The code needs healthy-control data that support the model in the full sample and all five training partitions. In particular, A′ requires at least 20 participants with two or more eligible assessments in each training partition. This is an implementation requirement, not a sufficient sample-size recommendation. A successful CSV check alone does not guarantee that the model can be estimated. See the [model settings](SETTINGS.md) for the design-dependent rules.

## 2. Choose the user settings

Source `run_analysis.R`, then call `run_cantab()`. Alternatively, edit the four settings at the top of [examples/run_one_outcome.R](../examples/run_one_outcome.R) and source that file from the repository root.

| Argument | Meaning | Default / allowed values |
|---|---|---|
| `outcome` | One outcome folder name | Required; use the six names in the README |
| `input_file` | Prepared input CSV | Required; an existing file |
| `output_dir` | Directory for this run's results | Required; outside the code repository |
| `mode` | What to do | `"check"` by default; `"fresh"` or `"resume"` |
| `reference_dir` | Optional separate participant-selection settings | `NULL`; usually leave it unspecified |

## 3. Check, fit, or resume

```r
source("run_analysis.R")
input_file <- "../input/CANTAB_connect_RS6_rawscores.csv"

# Input, all training-fold model designs, dependencies and Stan syntax.
run_cantab("SWM_Strategy", input_file, "../output/strategy_check",
           mode = "check")

# Fit the model and calculate normative scores.
run_cantab("SWM_Strategy", input_file, "../output/strategy_fit",
           mode = "fresh")

# Only after interruption: same fit directory, input, code, settings and software.
run_cantab("SWM_Strategy", input_file, "../output/strategy_fit",
           mode = "resume")
```

`check` creates diagnostic files but **does not create normative scores**. `fresh` needs a new output directory; do not reuse the check directory. `resume` reuses verified completed stages from an interrupted fitting run. Changes to input values/order, code, reference tables, package versions or settings require a new directory and a fresh run. Sourcing `run_analysis.R` only loads functions.

Every outcome folder also has `run.R`; source it and call `run_outcome()` with the same arguments except `outcome`. Keep the complete repository because these entry points use shared helpers.

## 4. Open the results

After a successful fit, the output directory contains:

| File | Purpose |
|---|---|
| `normative_scores.csv` | Main result, one row per assessment in the selected analysis input |
| `score_dictionary.csv` | Definition of every result column |
| `analysis_summary.txt` | Model description, score direction and number of available scores |

These names are the same for all outcomes. Use a different output directory for each outcome. Read the [diagnostic guide](DIAGNOSTICS.md) for the exact status/log filenames, convergence checks, score precision and held-out calibration files to inspect before interpreting the scores.

```r
scores <- read.csv("../output/strategy_fit/normative_scores.csv",
                   colClasses = c(Subject_ID = "character"))
scores[, c("Subject_ID", "Timepoint", "raw_score", "normative_z", "score_status")]
```

| Result column | Meaning |
|---|---|
| `Subject_ID`, `Diagnosis`, `Timepoint`, `Platform` | Participant and assessment identity |
| `Gender`, `FSIQ`, `Age_baseline_years`, `Time_since_baseline_years` | Model input covariates; age/time are in years |
| `outcome` | The measure being scored |
| `raw_score` | Original outcome value; no sign reversal or extra transformation |
| `normative_z` | Normative z score for this assessment |
| `score_status` | `available`: usable score; `not_available`: row retained without a usable score |
| `reference_model` | `cross_validation`: this healthy participant was excluded from training; `all_healthy_controls`: all eligible healthy controls were used |
| `score_direction` | `higher_is_better` for A′; `higher_is_worse` for the other five outcomes |

Within the selected analysis input, rows without a score remain in the file, with an empty `normative_z`. When `reference_dir` is omitted, this includes every original input row. When it selects a cohort, participants absent from its reference table are excluded before scoring; `COHORT_SELECTION.csv` records every source row and whether it was included, and `INPUT_AUDIT_R.json` records the counts. Do not replace missing scores with 0. The reference-model field is empty on these rows. Time outcomes remain in milliseconds in `raw_score`; A′ remains on 0–1; SWM values retain their input scale.

## 5. Interpret a normative z score

A normative model estimates the distribution of scores expected for the supplied covariates, including age, sex/gender, IQ, assessment timing and platform. Each observed value is located within that predictive distribution and expressed on the standard normal scale.

- `0`: cumulative probability 0.5; for continuous distributions this is the predicted median.
- About `+1.645`: cumulative probability 0.95; about `−1.645`: cumulative probability 0.05.
- Higher A′ is better. Higher reaction/movement time, latency, error count or strategy score is worse. The software preserves these directions.

For discrete scores and boundary point masses, the probability is the midpoint of the interval assigned to that observed value: `P(Y < y) + 0.5 * P(Y = y)`. A raw median therefore need not receive a score of zero, and some z scores may be unattainable.

A normative score is not simply `(value − sample mean) / sample SD`. It is a per-assessment score, not a longitudinal change score or an individual diagnosis. Check convergence and held-out calibration before interpreting a fitted model, using the [diagnostic files and explanations](DIAGNOSTICS.md). A successful software check is not evidence of model suitability.

## 6. Models and fixed estimation settings

- **Reaction time:** Uses a Student-t distribution for log(1 + milliseconds). Participant-specific intercepts can differ by platform, with a common participant time slope.
- **Movement time:** Uses a Student-t distribution for log(1 + milliseconds). Residual variability can differ by platform and assessment occasion.
- **Attention sensitivity (A prime):** Uses an extended-support beta distribution with observations censored to the interval from 0 to 1. Exact values of 0 and 1 are supported.
- **Response latency:** Uses a Gamma distribution directly on positive milliseconds, with a log link for the mean and platform/occasion effects on its shape.
- **Spatial working-memory errors:** Separates zero errors from positive counts and uses a discretized normal distribution for positive counts. Trial exposure is fixed at RS6 = 4 and Connect = 1. Confirm that this 4:1 ratio matches the test protocol underlying your exported counts; input counts are not divided by four.
- **Spatial working-memory strategy:** Combines additional probability at the minimum with a discretized beta distribution. RS6 scores are divided by four inside the model; the different score grids are retained.

Each reference model uses healthy-control data (`HC`). Healthy participants are divided into five groups **by participant**: all visits from one person stay together. Each group is scored with the model trained on the other four. Non-control participants use the model fitted to all eligible healthy controls.

The following are initial default settings. Iterations are **per chain**, and the first number **includes** warmup; warmup is adaptation and is not retained for inference. All models use four chains. The second model column applies to each of the five cross-validation fits.

| Outcome | All-control fit: total / warmup | Each cross-validation fit: total / warmup | `adapt_delta` | `max_treedepth` |
|---|---|---|---|---|
| `RTI_reaction` | 5000 / 1500 | 3000 / 1000 | 0.999 | 20 |
| `RTI_movement` | 5000 / 1500 | 6000 / 2000 | 0.999 | 20 |
| `RVP_A_prime` | 5000 / 1500 | 5000 / 1500 | 0.9995 | 20 |
| `RVP_latency` | 5000 / 1500 | 3000 / 1000 | 0.999 | 20 |
| `SWM_Between_errors` | 3000 / 1000 | 3000 / 1000 | 0.995 | 18 |
| `SWM_Strategy` | 4000 / 1500 | 4000 / 1500 | 0.999 | 12 |

`adapt_delta` is a target acceptance setting for the sampler; `max_treedepth` limits its trajectory-building depth. Prespecified recovery attempts may increase some settings if a fit has sampling problems. The [settings guide](SETTINGS.md) explains these rules and where the run records are stored; diagnostic CSVs do not all contain the same settings fields.

For formulas, priors, automatic changes when data have limited variation, and sampling controls, read the [model settings guide](SETTINGS.md). Movement time and latency can use a participant random intercept instead of a random slope under their defined fallback; A′ builds age/time and covariate terms from the variation in each training set. The per-outcome `provenance/model.json` and [outcomes.csv](outcomes.csv) describe the defaults. They are **reference documentation, not editable run configuration**. The supported user interface selects the outcome, paths and mode. Changing model formulas, priors or sampling settings is an advanced source-code modification requiring fresh validation and a new run.

## 7. Terms and common errors

| Term | Meaning here |
|---|---|
| Healthy controls / HC | Reference participants labelled `HC` |
| Covariate | An input used to describe the expected score distribution |
| Cross-validation fold | One of five participant groups held out from model training |
| Held-out / out-of-fold | The scored participant's data were excluded from that reference fit |
| MCMC / chain | Sampling used to estimate model parameters; a chain is one sampling sequence |
| Warmup | Initial sampler adaptation, excluded from retained estimates |
| Convergence / calibration | Whether sampling is reliable / predictions behave appropriately in held-out controls |
| Stan / CmdStan | Software used to compile and sample the statistical models |

| Message or situation | Action |
|---|---|
| Missing columns | Use the input template; preserve hyphens with `check.names = FALSE` |
| Output directory already exists | Use a new directory for a new run; use `resume` only for an interrupted fit |
| Input/code/settings changed | Start a new run rather than combining incompatible results |
| Locked output directory | Confirm the previous R process has stopped before removing its stale lock |
| Missing package / CmdStan / compiler | Complete the environment setup and rerun `check` |
| Insufficient variation / nonidentifiable model design | Inspect healthy-control coverage across covariates, platform and visits; do not bypass the check |
| `check` succeeded but no score CSV | Expected: use `mode = "fresh"` to fit |
| Score status is `not_available` | Check the outcome and required covariates for that assessment; use the [diagnostic guide](DIAGNOSTICS.md) to find eligibility/precision details. Retain the missing value |

The optional `reference_dir` is explained at the end of the [input guide](DATA.md). It is unnecessary for the default workflow.
