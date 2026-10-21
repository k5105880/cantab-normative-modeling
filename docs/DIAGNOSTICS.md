# Reading run status and model diagnostics

Start with the [execution guide](USAGE.md). This page identifies the files to inspect after `check`, during a long run, or before interpreting fitted scores.

Every path below is relative to the `output_dir` supplied to `run_cantab()`. Files appear as their stages finish: a `check` run has no posterior convergence results or normative scores, and an interrupted run may not yet have an aggregate diagnostic table.

## 1. Find the run status and log

| Outcome | Status file | Successful check / successful fitting run | Log |
|---|---|---|---|
| `RTI_reaction` | `PIPELINE_STATUS.txt` | `CHECK_PASSED_NO_MCMC` / `NORMATIVE_COMPLETE` | `logs/` contains timestamped `*_model_*.log` |
| `RTI_movement` | `run_status.txt` | `CHECK_COMPLETE: no MCMC or normative scores.` / `COMPLETE` | `console_<mode>_<timestamp>_<pid>.log` |
| `RVP_A_prime` | `run_status.json` | `CHECK_PASSED_NO_MCMC` / `COMPLETE` in the `status` field | `check.log`, `prepare.log`, `fit.log`, `score.log`, `prior.log`, `calibrate.log`, `verify_models.log` |
| `RVP_latency` | `run_status.txt` | `CHECK_COMPLETE: no MCMC or normative scores.` / `COMPLETE` | `console_<mode>_<timestamp>_<pid>.log` |
| `SWM_Between_errors` | `PORTABLE_RUN_STATUS.txt` | `CHECK_COMPLETE: no MCMC or final scores.` / `NORMATIVE_COMPLETE` | `console_<mode>_<timestamp>_<pid>.log`; fitting details in `logs/normative.log` |
| `SWM_Strategy` | `PORTABLE_RUN_STATUS.txt` | `CHECK_COMPLETE: no MCMC or final scores.` / `NORMATIVE_COMPLETE` | `console_<mode>_<timestamp>_<pid>.log` |

After successful SWM Between errors scoring, `diagnostic_summaries.zip` contains the files from that run's `aggregate/` directory. It is a convenience copy of diagnostic reports, not the public code release. It stays in the output directory; review reports and run metadata before sharing it.

`RUNNING...` records that execution started. After an abrupt termination it can be stale; confirm whether the R process is still active. `FAILED`, `FAILED_CHECK`, `FAILED_OR_INTERRUPTED` or `INTERRUPTED` mean that the launcher did not finish successfully. Read the named log and the error message before choosing a next step. If validation failed before a run was initialized, some of these files may not exist.

The table uses the outer launcher's status file. Inner files such as `model/STATUS.txt`, `model_run_status.txt` or `RUN_STATUS.txt` describe individual stages. An inner fitting stage can finish before score export finishes.

For example, from another R session while a run is active, or after it returns:

```r
run_dir <- "../output/strategy_fit"
readLines(file.path(run_dir, "PORTABLE_RUN_STATUS.txt"))
logs <- list.files(run_dir, pattern = "^console_.*[.]log$",
                   full.names = TRUE)
logs # Choose the file corresponding to the run you are inspecting.
if (length(logs)) {
  latest <- logs[which.max(file.info(logs)$mtime)]
  cat(tail(readLines(latest, warn = FALSE), 40), sep = "\n")
}
```

## 2. Read the input and design checks

All outcomes save `analysis_input.csv`, the selected analysis data, and `INPUT_AUDIT_R.json`, with source and selected counts. `COHORT_SELECTION.csv` identifies each original input row and its `included` value. Without `reference_dir`, all source rows are selected. Selection is separate from outcome eligibility: a selected row can still lack the outcome or covariates required to calculate a score.

| Outcome | Files produced by the input/design stage |
|---|---|
| `RTI_reaction` | `model/input_counts.csv`, `model/design_checks.csv`, `model/design_audit_full.csv` |
| `RTI_movement`, `RVP_latency` | `preflight/primary_subject_folds.csv`, `preflight/standata_rank.csv` |
| `RVP_A_prime` | `input/prepared_counts.csv`, `input/model_design_audit.csv`, `input/design_audit_full.csv` |
| `SWM_Between_errors` | `aggregate/input_group_counts.csv`, `aggregate/design_ranks.csv`, `aggregate/fold_balance.csv`, `aggregate/PREFLIGHT_PASS.txt` |
| `SWM_Strategy` | `input_quality_control.csv`, `design_ranks.csv`, `healthy_control_folds.csv`, `fold_balance_audit.csv` |

Counts distinguish assessments from participants. A fold is a group of participants; all their visits stay together. Design-rank checks ask whether the model's predictor columns contain enough independent information in the full sample and each training fold. A failed design check calls for reviewing the data and the [model's requirements](SETTINGS.md).

## 3. Check sampling and score precision

Convergence concerns estimation of model parameters. Numerical integration concerns calculating a predictive probability from an accepted fit. Posterior Monte Carlo standard error (MCSE) measures remaining simulation error in that calculated z score; it is not a participant's confidence interval or a raw-score standard deviation.

| Outcome | Sampling diagnostics | Score-precision details |
|---|---|---|
| `RTI_reaction` | `model/sampling_diagnostics.csv`; earlier per-attempt files in `model/fits/` | `model/numerical_quality.csv`, `model/normative_prediction_details.csv` |
| `RTI_movement` | `full_fit_diagnostics.csv`, `cv_diagnostics_*.csv` | `RTI_movement_student_log1p_full_population_state_cdf_diagnostics.csv`, `RTI_movement_student_log1p_oof_population_state_cdf_diagnostics.csv` |
| `RVP_A_prime` | `models/model_diagnostics.csv`, `models/actual_verification_diagnostics.csv`; actual sampling settings in `models/model_manifest.csv` | `results/normative_score_precision.csv`; a failed scoring batch can leave `score_cache/*_UNRESOLVED.csv` |
| `RVP_latency` | `full_fit_diagnostics.csv`, `cv_diagnostics_*.csv` | `RVP_latency_gamma_raw_full_population_state_cdf_diagnostics.csv`, `RVP_latency_gamma_raw_oof_population_state_cdf_diagnostics.csv` |
| `SWM_Between_errors` | `aggregate/all_model_diagnostics.csv` | `private/normative_fold1_state_precision.csv` through `private/normative_fold5_state_precision.csv`; `private/normative_full_state_precision.csv` when non-HC rows are scored |
| `SWM_Strategy` | `all_fit_diagnostics.csv`; per-fit initial/retry files also remain | `main_fold1_state_precision.csv` through `main_fold5_state_precision.csv`; `main_other_participants_chunk1_state_precision.csv` and subsequent chunks when non-HC rows are scored |

Here `*` means a filename component generated by the run. For movement/latency, `cv_diagnostics_*.csv` includes the role and seed in its filename. `full_population` contains scores from the full-HC reference model; `oof_population` contains scores for held-out HC.

The principal numeric sampler gates implemented by this release are:

| Quantity | Meaning | Implemented requirement for accepted posterior fits |
|---|---|---|
| `max_rhat` | Largest across-chain convergence statistic among monitored parameters | `< 1.01` |
| `min_ess_bulk`, `min_ess_tail` | Smallest effective sample sizes for the body and tails of monitored parameter distributions | Both `>= 400` |
| `n_divergent` or `divergences` | Divergent sampler transitions | `0` |
| `n_max_treedepth` or `treedepth_hits` | Transitions reaching the permitted tree depth | `0` |
| `min_bfmi`, `min_ebfmi` or `min_EBFMI` | Minimum chain energy-mixing diagnostic | `> 0.3` for A′ and Between errors; `>= 0.3` for the other four |

Required diagnostic values must also be finite. A′ verifies complete chain draws and sampler records; the SWM workflows also report parameters without finite diagnostics. Use each workflow's computed `pass` or `validation_pass` field where present. A′ and Between errors enforce their gates without a common `pass` column in the aggregate table. Reaction-time attempt tables can contain a rejected initial attempt followed by an accepted retry; inspect the attempt identifier. [Settings and recovery rules](SETTINGS.md) explain which retries are permitted.

Scoring additionally requires successful integration and finite posterior MCSE(z). Its limit is `0.02` for movement time and latency, and `0.03` for the other four outcomes. The corresponding fields are `numerical_ok` and `posterior_precision_ok`, reaction time's `posterior_mc_ok`, or A′'s `precision_status = "ok"`. A′ also requires successive quadrature z estimates to differ by at most `0.002`; the other scorers check their integration tolerances separately. Do not infer numerical success from sampler convergence alone.

A scoring failure can stop final export. Keep its precision files and log, identify whether the issue is integration or posterior sampling, and review the settings before any new fit.

## 4. Examine held-out prediction, not only successful computation

| Outcome | Calibration files to open | What to examine |
|---|---|---|
| `RTI_reaction` | `model/HC_OOF_calibration.csv` | Counts, z mean/SD, tail proportions and `central95`, by platform/visit/age; `weighting` distinguishes equal assessments from equal participants |
| `RTI_movement`, `RVP_latency` | `OOF_calibration.csv`, `OOF_calibration_by_platform.csv`, `OOF_calibration_by_practice.csv`, `OOF_calibration_by_covariate.csv`, `OOF_interval_coverage.csv` | Held-out HC z distributions overall and within subgroups; nominal interval level versus observed raw-score predictive coverage |
| `RVP_A_prime` | `results/OOF_calibration.csv`, `results/calibration_methods.txt` | Integrated randomized-PIT moments/tails and participant-cluster bootstrap intervals |
| `SWM_Between_errors` | `aggregate/calibration_by_stratum.csv`, `aggregate/platform_gaps_bootstrap.csv` | Integrated randomized-PIT moments, observed/predicted zero proportions and interval coverage; Connect-minus-RS6 difference with bootstrap intervals |
| `SWM_Strategy` | `held_out_calibration.csv` | Mid-P z moments, randomized-PIT moments, observed/predicted minimum-score proportions and `coverage95`, with group counts |

PIT means *probability integral transform*: it locates observations in their predictive distributions. With discrete values or probability masses at a boundary, randomized-PIT summaries account for the entire probability interval of a value. They can differ from the mean/SD of the exported deterministic mid-P z scores.

These summaries describe held-out HC from the supplied cohort. Read `n`/`n_rows` and `subjects`/`n_subjects` alongside subgroup estimates; an empty or small group gives limited information. A′'s `PIT_central95` and reaction time's `central95` describe central probability/score mass; they are not raw-score predictive-interval coverage. Between-errors columns named `within_illustrative_margin_*` use illustrative margins, not an established acceptance criterion. Bootstrap intervals condition on the fitted predictions; they do not establish validity in a new population.

There is no universal calibration pass threshold implemented across these six workflows. Their sampler/precision gates protect computational reliability; they do not establish that the chosen model is suitable for your data. Interpret calibration in relation to the test protocol, covariate coverage and intended use.

## 5. First troubleshooting steps

| Situation | First action |
|---|---|
| `check` completed but scores are absent | Expected. Use `fresh` with a new output directory to fit. |
| Input or design error | Read the first error in the log. Compare units, column names, missingness and eligible HC/fold counts with the [input guide](DATA.md). |
| `Malformed CSV` | Check the source export's delimiter, header count and quotation marks. Re-export a complete UTF-8 CSV; do not suppress the error or fill missing fields automatically. |
| Old `brms`, missing CmdStan or compiler | Follow [environment setup](../environment/README.md), then rerun `check` in a new directory. |
| A completed result contains `not_available` | Match `Subject_ID` and `Timepoint` to the selected input. Check the required outcome/covariates; Strategy additionally reports `eligible` and the score-range reason in `input_quality_control.csv`. Keep the missing score. |
| Sampler or score-precision failure | Open the files in section 3. Preserve the evidence; do not overwrite failed outputs or change gates merely to obtain scores. |
| Interrupted run | Confirm the previous R process stopped. Resume the same fitting directory only with unchanged input, code, reference tables and software. |
| Locked output | The error names the lock. Verify that no process is using the directory before removing a stale lock. A status file alone is insufficient to establish this. |
| Input, model settings or software changed | Use a new fitting directory and validate again. |

For deeper inspection, the common result's source tables are `analysis/RTI_analysis_ready.csv`, `RTI_movement_scores.csv`, `results/RVP_A_prime_scores.csv`, `RVP_latency_scores.csv`, `private/SWM_scores.csv`, and `SWM_Strategy_scores.csv`, respectively. They retain implementation-specific eligibility or precision fields. The common `normative_scores.csv` remains the main score table.
