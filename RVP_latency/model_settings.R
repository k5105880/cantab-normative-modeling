# Defaults for the normative model. Parallelism does not change the sampling plan.
# See ../docs/SETTINGS.md for every setting, units and fixed constraints.
# Participant-held-out CV: keep all visits together; balance uses baseline covariates.
cv_K <- 5
# Joint likelihood by participant; not a separate model-selection option.
cv_joint_likelihood <- "group"
cv_balance_folds <- TRUE
cv_balance_candidates <- 500
# This affects fold balance only; FSIQ still enters the model when it varies.
cv_balance_include_fsiq <- FALSE
# Separate prior-only sampler. Iteration counts include discarded warmup.
prior_pred_chains <- 2
prior_pred_iter <- 1200
prior_pred_warmup <- 600
# Selected prior draws for predictive summaries, capped at available draws.
prior_pred_draws <- 1000
# Full healthy-control fit; iteration counts are per chain and include warmup.
chains_fit <- 4
iter_fit <- 5000
warmup_fit <- 1500
# core.R uses one worker; all sampling chains and iterations are preserved.
# Each held-out fit uses the following chain/iteration/warmup settings.
cv_chains <- 4
cv_iter <- 3000
cv_warmup <- 1000
# NUTS controls shared by full/CV/prior fits; depth also enters diagnostic gates.
adapt_delta_fit <- 0.999
max_treedepth_fit <- 20
# brms cache policy; the outer run identity still rejects changed runs on resume.
file_refit_mode <- "on_change"
# All posterior draws score normative z; numerical tolerances/chunks: state_cdf.R.
# Fixed role-specific offsets are added to this seed.
seed_main <- 1234
# Center and scale age/IQ; changing these alters coefficient/prior interpretation.
reference_age_years <- 50
age_scale_years <- 10
reference_fsiq <- 100
fsiq_scale_points <- 15
# Smooth age/time require >=6/>=5 unique values respectively; otherwise linear.
age_term_type <- "smooth"
time_term_type <- "smooth"
# Terms are also conditional on levels present. Platform mean gates its interactions.
include_platform_mean <- TRUE
include_platform_time_interaction <- TRUE
include_platform_practice_interaction <- TRUE
# Scale terms model Gamma shape, not the mean response.
include_scale_platform <- TRUE
include_scale_practice <- TRUE
# Random time slopes need enough repeated participants; otherwise intercept only.
random_effect_mode <- "random_slope"
# Retry after a model fitting error, not a relaxation of convergence diagnostics.
allow_random_slope_fallback <- TRUE
min_subjects_with_2plus_visits_for_random_slope <- 20
# LKJ prior for participant intercept/slope correlation, when slopes are present.
prior_cor_eta <- 2
# log(mean_ms) is the Normal prior location for log mean response, not its raw mean.
rvp_latency_prior_intercept_mean_ms <- 500
rvp_latency_prior_intercept_log_sd <- 0.35
# Zero-centered Normal SDs for mean-response and log-shape coefficients.
rvp_latency_prior_b_sd <- 0.25
rvp_latency_prior_aux_b_sd <- 0.25
# Exponential RATE for participant/smooth SDs; higher rate favors smaller effects.
rvp_latency_prior_sd_rate <- 5
rvp_latency_prior_sds_rate <- 5
# Normal prior for log Gamma shape, centered at log(shape_median).
rvp_latency_prior_shape_median <- 25
rvp_latency_prior_shape_log_sd <- 0.6
