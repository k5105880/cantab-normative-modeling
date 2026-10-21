# Settings and source-file roles

[Run an analysis](USAGE.md) · [Inspect diagnostics](DIAGNOSTICS.md)

## Start with these four settings

For ordinary use, keep the model defaults and edit only the four assignments in [examples/run_one_outcome.R](../examples/run_one_outcome.R), or pass them to `run_cantab()`:

| Setting | What to enter |
|---|---|
| `outcome` | Exactly one of the six folder names in the README |
| `input_file` | Your prepared CSV, for example `../input/CANTAB_connect_RS6_rawscores.csv` |
| `output_dir` | A separate directory for this outcome and run, for example `../output/reaction_fit` |
| `mode` | `"check"` to validate; `"fresh"` to fit; `"resume"` only to continue an interrupted fit |

`source("examples/run_one_outcome.R")` runs the example in check mode. Before fitting, set `mode = "fresh"` and choose a new `output_dir`.

The optional `reference_dir` filters participants and, for two outcomes, can supply fold assignments. It does not provide a trained reference model. Leave it as `NULL` unless using the [documented cohort settings](DATA.md#optional-cohort-and-fold-settings).

Set `CMDSTAN` only when selecting a custom CmdStan installation as described in the [setup guide](../environment/README.md). Internal environment variables are set automatically by the launcher.

## Which files do what?

| Files | Role |
|---|---|
| `run_analysis.R`, `examples/run_one_outcome.R` | Public interface and editable execution example |
| Each outcome's `run.R` | Sourceable entry point for that single outcome; call `run_outcome()` after sourcing |
| `RTI_reaction/config.R`, `SWM_Between_errors/config.R`, `SWM_Strategy/config.R` | Active model, sampling and diagnostic constants |
| `RTI_movement/model_settings.R`, `RVP_latency/model_settings.R` | Active model specification and sampling constants |
| `RVP_A_prime/config.json` | Active A′ configuration; validated by `R/config.R` |
| `state_cdf.R`, `R/state_precision.R` | Predictive probability integration and score-precision constants |
| `core.R`, `model_functions.R`, `R/model_spec.R`, `vendor/*definition*.R` | Data adaptation, model formulas and prior distributions; exact roles are listed in the outcome READMEs |
| `main.R`, `run_model.R`, outcome-internal `run_analysis.R` | Implementation scripts invoked by the launcher in its execution process; not additional setup steps |
| Each outcome's `provenance/model.json`, `docs/outcomes.csv` | Descriptions of the supplied defaults; editing these does not change calculations |
| `environment/recorded_R_packages.csv` | Recorded package versions; not an installer or dependency lockfile |
| `tests/`, `test_*.R` | Software checks with synthetic inputs; start them with `source("scripts/run_tests.R")` |

The remaining sections explain the active constants for readers inspecting or adapting the implementation. They are **not additional inputs needed to start an analysis**. Changing formulas, priors, numerical tolerances, seeds or sampling plans creates a modified model implementation. Such a change needs appropriate validation and a new output directory; it cannot be applied with `resume`. Loading a configuration in the interactive console does not override the copy loaded in the child R session.

## Terms used in the settings

| Term | Meaning and units |
|---|---|
| `iter`, `warmup` | Iterations per chain; `iter` includes warmup. Warmup is discarded. Retained draws per chain are `(iter - warmup) / thin`. |
| `adapt_delta`, `max_treedepth` | Sampler target acceptance probability and maximum tree depth; these are not score cutoffs. |
| `seed`, `fold_seed` | Reproducible random-number seeds, not dates or participant-selection criteria. Different stages use fixed offsets. |
| `folds`, `cv_K` | Five participant groups for cross-validation; every visit from a person stays in one group. |
| Prior / posterior | Model parameter distribution before / after fitting the supplied healthy-control data. |
| `draws`, `ndraws` | Number of parameter or simulated-observation draws for the specified step; the purpose differs by field below. |
| Quadrature / `nodes` | Deterministic numerical integration over new-participant effects; not the number of people or MCMC draws. |
| MCSE | Monte Carlo standard error of the estimated z score; computational precision, not a participant's clinical uncertainty. |
| `chunk`, `batch_size` | Number of rows or draws processed together to manage memory. See each field for its unit. |
| `state` in numerical helper names | A score at one assessment, obtained from a new-participant predictive distribution. The public column is `normative_z`. |
| `sd`, `rate`, `shape`, `eta` | Distribution parameters. An exponential/Gamma rate is the inverse of scale; `sd` is a standard deviation. |

The initial sampling plans are listed in [USAGE.md](USAGE.md#6-models-and-fixed-estimation-settings). The settings below also describe recovery attempts and diagnostic calculations.

## RTI reaction

Active object: `rti_config` in [RTI_reaction/config.R](../RTI_reaction/config.R).

| Field(s) | Supplied value | Meaning |
|---|---|---|
| `version` | `reaction_normative_v1` | Cache/format identity; does not select a different model |
| `seed`; `folds` | `20260921`; `1:5` | Base random seed; required fold IDs. Fold construction uses seed 1334 and 500 balanced candidates in `run.R`. |
| `chains` | 4 | Number of sampling chains |
| `full_iter`, `full_warmup` | 5000, 1500 | All-control fit |
| `cv_iter`, `cv_warmup` | 3000, 1000 | Each held-out fit |
| `adapt_delta`, `max_treedepth` | 0.999, 20 | Initial sampler settings |
| `retry_adapt_delta`, `retry_iteration_multiplier` | 0.99999, 2 | One failed-diagnostic retry doubles both total iterations and warmup; depth stays 20 |
| `max_rhat`, `min_ess`, `min_bfmi` | 1.01, 400, 0.3 | Acceptance requires R-hat **below** 1.01, bulk/tail ESS at least 400 and BFMI at least 0.3; zero divergences/depth hits are also required |
| `posterior_mcse_z_limit` | 0.03 | Maximum accepted score MCSE. Final export independently enforces a ceiling of 0.03; increasing this value cannot relax that ceiling. |
| `quad_nodes`, `check_nodes`, `quad_max_nodes` | 40, 80, 320 | Initial, comparison and maximum integration orders |
| `quad_z_tolerance`, `quad_log_score_tolerance` | 0.005, 0.0001 | Allowed change in z and log predictive density between integration resolutions |
| `cdf_chunk_size` | 40 | Assessment rows per integration batch |
| `run_prior_predictive`; `prior_draws` | `TRUE`; 1000 | Run the prior predictive check; simulate from 1000 draws. Its separate prior-only sampler uses 2 chains, 1500 total/750 warmup. |

## RTI movement and RVP latency

Active files: [movement model_settings.R](../RTI_movement/model_settings.R) and [latency model_settings.R](../RVP_latency/model_settings.R).

| Field(s) | Supplied value | Meaning |
|---|---|---|
| `cv_K`, `cv_joint_likelihood` | 5, `"group"` | Five folds; joint predictive likelihood is evaluated by participant |
| `cv_balance_folds`, `cv_balance_candidates` | `TRUE`, 500 | Select a balanced partition from 500 candidate allocations |
| `cv_balance_include_fsiq` | `FALSE` | IQ is not used in the partition-balance objective; this does **not** remove IQ from the model |
| `seed_main` | 1234 | Base seed for fitting, folds, predictive checks and calibration streams |
| `chains_fit`, `iter_fit`, `warmup_fit` | 4, 5000, 1500 | All-control sampling |
| `cv_chains`, `cv_iter`, `cv_warmup` | Movement: 4, 6000, 2000; latency: 4, 3000, 1000 | Cross-validation sampling. Movement's iteration values must match its fixed recovery plan below. |
| `adapt_delta_fit`, `max_treedepth_fit` | 0.999, 20 | Full, held-out and prior-only sampler settings |
| `file_refit_mode` | `"on_change"` | Internal brms cache behavior; outer run checks still require a fresh directory when inputs/code/settings change |
| `prior_pred_chains`, `prior_pred_iter`, `prior_pred_warmup`, `prior_pred_draws` | 2, 1200, 600, 1000 | Separate prior-only fit and the number of draws used in its predictive summaries |
| `reference_age_years`, `age_scale_years` | 50, 10 | Model age = `(baseline age - 50) / 10`; not eligibility limits |
| `reference_fsiq`, `fsiq_scale_points` | 100, 15 | Model IQ = `(FSIQ - 100) / 15`; input the original IQ |
| `age_term_type`, `time_term_type` | `"smooth"`, `"smooth"` | Smooth age/time terms when training data have at least 6 distinct ages/5 distinct times; otherwise linear terms |
| `include_platform_mean` | `TRUE` | Platform mean term when two platform levels are present; also enables applicable scale/interaction terms |
| `include_platform_time_interaction` | `TRUE` | Platform-specific linear time contribution |
| `include_platform_practice_interaction` | `TRUE` | Platform-by-occasion term when both factors vary |
| `include_scale_platform`, `include_scale_practice` | `TRUE`, `TRUE` | Platform/occasion terms for residual sigma (movement) or Gamma shape (latency), where levels permit |
| `random_effect_mode` | `"random_slope"` | Request correlated participant intercept/time slope when supported |
| `min_subjects_with_2plus_visits_for_random_slope` | 20 | Minimum repeated participants in the training data for the slope term; otherwise use a random intercept |
| `allow_random_slope_fallback` | `TRUE` | A failed full-model fitting call may retry with random intercepts, except compilation/backend failures; convergence failures do not bypass their checks |
| `prior_cor_eta` | 2 | LKJ prior parameter for intercept/slope correlation |

Score integration uses **all retained posterior draws**. Its single batch-size setting is `state_cdf_settings$chunk_size` (25 rows), documented below.

Movement's held-out plan in `rti_execution_helpers.R` has one recovery attempt after a CV convergence-gate failure: 6000 total/2000 warmup/thin 1, then 10000/2000/thin 2. Sampling, design and other validation errors stop immediately. Both attempts retain 4000 draws per chain for the cross-validation prediction matrix. Changing only `cv_iter` or `cv_warmup` conflicts with that plan and stops execution. Latency has no corresponding longer cross-validation retry.

Movement prior constants (the `RTI_MT_` prefix means RTI movement time):

| Field(s) after `RTI_MT_prior_` | Value | Prior role |
|---|---|---|
| `intercept_mean`, `intercept_sd` | `log1p(300)`, 0.2 | Normal location/SD for the mean intercept on the log(1 + milliseconds) scale |
| `b_sd` | 0.15 | SD of zero-centered Normal mean coefficients |
| `sds_sd` | 0.15 | Positive-constrained Normal scale for smooth-effect SDs |
| `subject_intercept_sd`, `subject_slope_sd` | 0.2, 0.1 | Positive-constrained Normal scales for participant intercept/slope SDs |
| `sigma_intercept_mean`, `sigma_intercept_sd` | `log(0.2)`, 0.3 | Normal location/SD for log residual scale |
| `sigma_b_sd` | 0.2 | Normal coefficient SD for log residual scale |
| `nu_shape`, `nu_rate` | 5, 0.5 | Gamma shape/rate prior for Student-t degrees of freedom |

Latency prior constants (prefix `rvp_latency_prior_`):

| Field(s) after the prefix | Value | Prior role |
|---|---|---|
| `intercept_mean_ms`, `intercept_log_sd` | 500, 0.35 | Normal location `log(500)` and SD 0.35 for the log mean intercept; 500 is not an arithmetic prior mean after exponentiation |
| `b_sd`, `aux_b_sd` | 0.25, 0.25 | Normal coefficient SDs for log mean and log shape |
| `sd_rate`, `sds_rate` | 5, 5 | Exponential rates for participant-effect and smooth-effect SDs |
| `shape_median`, `shape_log_sd` | 25, 0.6 | Normal location `log(25)` and SD 0.6 for log shape; 25 is the median on the exponentiated scale |

## RVP A prime

Active file: [RVP_A_prime/config.json](../RVP_A_prime/config.json). A dotted name below means a nested JSON field, e.g. `sampling.iter`.

Use exactly the documented field names and JSON types. Unknown, duplicate or missing fields are rejected before analysis; numeric values must be numbers, not quoted strings. The `sampling_recovery` object is optional: omitting it or setting it to `null` disables the retry. All other listed fields are required.

| Field(s) | Value | Meaning |
|---|---|---|
| `schema_version`; `analysis` | 1; `RVP_A_prime_normative` | Configuration-format version and descriptive identity; not model selectors |
| `seed` | 20260922 | Base seed |
| `sampling.chains`, `.iter`, `.warmup` | 4, 5000, 1500 | Each full/held-out fit |
| `sampling.adapt_delta`, `.max_treedepth`, `.refresh` | 0.9995, 20, 200 | Sampler settings and progress-message interval |
| `sampling_recovery.trigger`, `.max_retries` | `"divergence_only"`, 1 | One retry only for the specified divergence-only failure; other diagnostic failures stop |
| `sampling_recovery.warmup`, `.adapt_delta` | 2500, 0.9999 | Retry settings; total iterations increase to 6000 to preserve 3500 retained draws per chain |
| `prior.draws`, `.batch_size` | 40000, 500 | Direct independent prior draws in batches of 500; no prior-only MCMC for this check |
| `calibration_bootstraps` | 2000 | Participant-level bootstrap repetitions for held-out calibration |
| `scoring.nodes` | `[32, 64, 128, 256]` | Required integration-resolution sequence |
| `scoring.tolerance`, `.mcse_limit` | 0.002, 0.03 | Maximum z change between resolutions and score MCSE |
| `scoring.chunk` | 8 | Assessment rows per integration task |

The validator enforces four chains, this integration sequence, at least 40000 prior draws, at least 2000 bootstraps, and precision thresholds no looser than those above. These restrictions implement the supplied calculation plan. Formula construction in `R/model_spec.R` adapts terms to training-data variation and requires at least 20 participants with repeated eligible visits in every full/fold training set; see the [data guide](DATA.md).

## SWM Between errors

Active object: `CONFIG` in [SWM_Between_errors/config.R](../SWM_Between_errors/config.R).

| Field(s) | Value | Meaning |
|---|---|---|
| `control_label`, `groups` | `"HC"`, `c("HC","MDD","BD")` | Training group and supported group order; match the shared input contract |
| `trials` | `c(RS6=4, Connect=1)` | Trial-exposure mapping in the statistical model. Confirm that the source protocol matches; supply original counts. |
| `models` | `"normative"` | Fixed implementation identity; other model IDs are rejected |
| `folds`, `fold_seed`, `fold_candidates` | 5, 1334, 500 | Number of folds, allocation seed and balanced candidates |
| `chains`, `iter`, `warmup` | 4, 3000, 1000 | Initial sampling |
| `adapt_delta`, `max_treedepth` | 0.995, 18 | Initial sampler settings |
| `retry_iter`, `retry_warmup`, `retry_adapt_delta`, `retry_max_treedepth` | 5000, 1500, 0.999, 20 | Second attempt after failed posterior diagnostics |
| `prediction_draws` | 8000 | Maximum draws for simulated interval-coverage diagnostics; normative z integration uses all retained draws |
| `bootstrap_reps` | 2000 | Participant-level bootstrap repetitions for platform-gap summaries |
| `prior_predictive` | `TRUE` | Separate prior-only fit: 1200 total/600 warmup, with up to 1200 draws for predictive summaries |
| `seed` | 202609210 | Base seed; stages use fixed offsets |

Definitions are in `vendor/common_model_definitions.R` and `vendor/final_model_definition.R`; the positive-count probability distribution is in `vendor/count_distribution.R`. The `vendor` directory contains the supporting source shipped with this workflow; it is not a separate installation step.

## SWM Strategy

Active object: `NORMATIVE_CONFIG` in [SWM_Strategy/config.R](../SWM_Strategy/config.R).

| Field(s) | Value | Meaning |
|---|---|---|
| `chains`, `iter`, `warmup` | 4, 4000, 1500 | Each full/held-out fit |
| `adapt_delta`, `max_treedepth` | 0.999, 12 | Initial sampler settings |
| `init_radius` | 0.1 | Sampler random-initialization radius; not a prior width |
| `refresh` | 25 | Iterations between progress messages |
| `ndraws` | 8000 | Sampled probability-mass/diagnostic draws, with replacement if necessary; normative z uses all retained draws deterministically |
| `seed`; `fold_seed`, `fold_candidates` | 20260921; 1334, 500 | Base seed; fold seed and balanced candidate count. Randomized calibration probabilities use a separate fixed seed 20260917. |
| `prediction_chunk_subjects` | 25 | Non-control participants per prediction batch; changing this may alter sampled diagnostic streams |

There is one divergence-only recovery attempt, with warmup raised to at least 2500 and `adapt_delta` to at least 0.9999. Retained draws per chain and maximum tree depth stay unchanged. Exact detection rules and saved diagnostics are in `R/swm_runtime.R` and `R/swm_diagnostics.R`.

## Numerical integration settings outside the main configuration

These constants govern score precision. They are read from the named source files automatically. The sampled diagnostic draw counts above do not reduce the posterior draws used for normative z.

| Constants | Movement / latency: `state_cdf_settings` in `state_cdf.R` | Both SWM: `swm_state_settings` in `R/state_precision.R` |
|---|---|---|
| `quad_nodes`, `check_nodes`, `max_nodes` | 32, 64, 1024 | 16, 32, 1024 |
| `z_tolerance` | 0.0001 | 0.0001 |
| `tail_relative_tolerance` | 0.0001 | 0.0001 |
| `draw_probability_tolerance` | 0.000001 | 0.000001 |
| `posterior_mcse_z_tolerance` | 0.02 | 0.03 |
| `chunk_size` | 25 assessment rows | 20 assessment rows |

Integration increases its resolution for unresolved rows. The three integration tolerances check absolute z change, the maximum relative tail change `abs(coarse / fine - 1)` over both tails, and the maximum absolute per-draw probability change over both tails. SWM uses midpoint tail probabilities for these comparisons. Each tolerance and the MCSE limit uses `<=`. These helpers also require the integration step relative to distribution scale (`quad_step_to_scale`) to be at most 1, and probability-normalization error to be strictly below `1e-10` for movement/latency or `1e-9` for SWM. Reaction has a separate normalization gate below `1e-10`.

The MCSE check separately assesses finite-posterior-sample error. Failing checks do not authorize treating an unreliable score as available. Use the [diagnostic guide](DIAGNOSTICS.md) to locate the corresponding results.
