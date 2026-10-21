# Synthetic workflow regression: folds, seeds, output routing and cache reuse.
# Fitting is replaced by a recorder; this test performs no MCMC.
args <- commandArgs(FALSE)
self <- sub("^--file=", "", args[grepl("^--file=", args)][1])
root <- dirname(dirname(normalizePath(gsub("~+~", " ", self, fixed = TRUE))))
for (file in c("config.R", "R/data.R", "R/models.R", "R/pipeline.R",
               "vendor/calibration.R", "R/reporting.R"))
    source(file.path(root, file))
stopifnot(identical(CONFIG$models, "normative"))
raw <- data.frame(
    Subject_ID = c(rep(sprintf("person_%02d", 1:5), each = 2), "person_06", "person_07"),
    Diagnosis = c(rep("HC", 10), "MDD", "BD"),
    Timepoint = c(rep(1:2, 5), 1, 1),
    Platform = c(rep(c("RS6", "Connect"), 5), "Connect", "RS6"),
    Gender = rep(c("F", "M"), 6), FSIQ = 105,
    Age_baseline_years = c(rep(seq(30, 70, 10), each = 2), 42, 55),
    Time_since_baseline_years = c(rep(c(0, 1), 5), 0, 0),
    SWM_Between_errors = c(0:10, NA))
p <- prepare_raw_data(raw, CONFIG)
p$hc$CV_fold <- rep(1:5, each = 2)
p$all$CV_fold <- c(rep(1:5, each = 2), NA, NA)
p$valid$CV_fold <- p$all$CV_fold[p$all$score_eligible]

# Empty non-control cohorts must return without entering the predictor.
stopifnot(nrow(predict_state(NULL,p$valid[FALSE,],"normative",1L,8000L))==0L)

trace <- list()
fit_fixed <- function(id, data, tag, seed, ctx, prior_only = FALSE) {
    trace[[length(trace) + 1L]] <<- list(stage = "fit", seed = seed,
        prior = prior_only, rows = data$row_id_global)
    list(fit = list(training_IDs = as.character(data$Subject_ID)),
         diagnostics = data.frame(pass = TRUE), fit_signature = "synthetic")
}
diagnostic_pass <- function(d) isTRUE(d$pass)
predict_state <- function(fit, data, id, seed, draws, coverage, diagnostics_file) {
    if(!nrow(data))return(data.frame())
    stopifnot(!length(intersect(fit$training_IDs, as.character(data$Subject_ID))))
    trace[[length(trace) + 1L]] <<- list(stage = "predict", seed = seed, rows = data$row_id_global)
    result <- data.frame(row_id_global = data$row_id_global, Subject_ID = as.character(data$Subject_ID),
        z = data$row_id_global / 10, posterior_precision_ok = TRUE,
        numerical_ok = TRUE, posterior_mcse_z = 0.001)
    result[c("Age_baseline_years", "Platform", "Practice_bin", "NM_response")] <-
        data[c("Age_baseline_years", "Platform", "Practice_bin", "NM_response")]
    result$pit_mean <- result$z
    result$pit_second <- result$z^2 + 0.1
    result$pit_low <- result$pit_high <- 0.025
    result$pred_zero <- 0.1
    result$log_pred_density <- -1
    for (level in c(50, 80, 95)) {
        result[[paste0("coverage_", level)]] <- TRUE
        result[[paste0("expected_coverage_", level)]] <- level / 100
    }
    result
}
summarize_prior_predictions <- function(fit, d, seed) {
    trace[[length(trace) + 1L]] <<- list(stage = "prior_summary", seed = seed)
    data.frame(zero_fraction = 0.1)
}
run <- tempfile("swm_workflow_")
dir.create(file.path(run, "private"), recursive = TRUE)
ctx <- list(run = run, signature = "synthetic", config = CONFIG)
job <- normative_job("normative", ctx, p)
stopifnot(length(trace) == 14L)
for (k in 0:5) {
    stopifnot(trace[[2*k+1]]$seed == CONFIG$seed + 1000 + k*10,
              trace[[2*k+2]]$seed == CONFIG$seed + 1001 + k*10)
}
stopifnot(trace[[13]]$seed == CONFIG$seed + 1500, trace[[13]]$prior,
          trace[[14]]$seed == CONFIG$seed + 1501,
          identical(sort(job$oof$row_id_global), 1:10),
          identical(job$full$row_id_global, 11L))
scores <- assemble_scores(p, list(normative = job), CONFIG)
stopifnot(identical(as.character(scores$Subject_ID), raw$Subject_ID),
          identical(scores$normative_z, c((1:11)/10, NA_real_)),
          all(scores$normative_score_source[1:10] == "HC_subject_held_out"),
          all(scores$normative_score_source[11:12] == "nonHC_full_HC"),
          is.na(scores$posterior_precision_ok[12]))
before_trace <- trace
stopifnot(identical(normative_job("normative", ctx, p), job),
          identical(trace, before_trace))
dir.create(file.path(run, "aggregate"))
ctx$config$bootstrap_reps <- 20L
report_normative(p, list(normative = job), scores, ctx)
expected_reports <- c("all_model_diagnostics.csv", "prior_predictive.csv",
    "calibration_by_stratum.csv", "platform_gaps_bootstrap.csv",
    "score_availability.csv", "score_dictionary.csv", "sample_counts.csv", "RESULTS_GUIDE.txt")
stopifnot(setequal(list.files(file.path(run, "aggregate")), expected_reports))
dictionary <- read.csv(file.path(run, "aggregate/score_dictionary.csv"))
stopifnot(identical(dictionary$column, "normative_z"))
fails <- function(expr) inherits(tryCatch({force(expr); NULL}, error = identity), "error")
bad <- job
bad$full$posterior_precision_ok[1] <- FALSE
stopifnot(fails(assemble_scores(p, list(normative = bad), CONFIG)))
bad$full$posterior_precision_ok[1] <- NA
stopifnot(fails(assemble_scores(p, list(normative = bad), CONFIG)))
bad <- job
bad$full$Subject_ID[1] <- "wrong-person"
stopifnot(fails(assemble_scores(p, list(normative = bad), CONFIG)),
          fails(fixed_definition("unsupported", p$hc)))
# A control-only input completes all five held-out predictions and scoring.
control_only <- prepare_raw_data(raw[1:10,], CONFIG)
control_only$hc$CV_fold <- control_only$all$CV_fold <- control_only$valid$CV_fold <- rep(1:5,each=2)
control_run <- tempfile("swm_controls_only_")
dir.create(file.path(control_run,"private"),recursive=TRUE)
control_ctx <- list(run=control_run,signature="synthetic_controls",config=CONFIG)
control_job <- normative_job("normative",control_ctx,control_only)
control_scores <- assemble_scores(control_only,list(normative=control_job),CONFIG)
stopifnot(nrow(control_job$full)==0L,nrow(control_job$oof)==10L,
          identical(control_scores$normative_z,(1:10)/10))
cat("PASS: normative workflow, five held-out folds, fixed seeds, missing rows, reports, precision gates and cache reuse.\n")
