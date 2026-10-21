# Synthetic checks for the public score tables; no model fitting.
local({
  root <- normalizePath(".", mustWork = TRUE)
  source(file.path(root, "R/input_data.R"), local = TRUE)
  stopifnot("Shared result export is available" = file.exists(file.path(root, "R/results.R")))
  source(file.path(root, "R/results.R"), local = TRUE)
  work <- tempfile("cantab_readable_results_"); dir.create(work)
  on.exit(unlink(work, recursive = TRUE), add = TRUE)
  raw <- data.frame(Subject_ID = c("001", "9007199254740993", "0002"),
    Diagnosis = c("HC", "MDD", "HC"), Timepoint = 1, Platform = c("RS6", "Connect", "RS6"),
    Gender = c("F", "M", "F"), FSIQ = c(110, 95, NA),
    Age_baseline_years = c(40, 52, 30), Time_since_baseline_years = 0,
    SWM_Between_errors = c(4, 0, NA), SWM_Strategy = c(32, 2, NA),
    RVP_A_prime = c(.9, 1, NA), RVP_Median_latency = c(400, 450, NA),
    check.names = FALSE)
  raw[["RTI_Median_five-choice_reaction_time"]] <- c(310, 350, NA)
  raw[["RTI_Median_five-choice_movement_time"]] <- c(220, 250, NA)
  z <- c("-1.2345678901234567", "2.3456789012345678", "")
  cases <- list(
    RTI_reaction = list(path = "analysis/RTI_analysis_ready.csv", score = "normative_z",
      raw = "RTI_Median_five-choice_reaction_time", source = "normative_source",
      values = c("HC_subject_OOF", "patient_full_HC_new_subject", ""),
      checks = list(normative_usable = c("TRUE", "TRUE", "FALSE"))),
    RTI_movement = list(path = "RTI_movement_scores.csv",
      score = "normative_z", source = "score_source",
      values = c("HC_OOF", "full_HC_new_subject", ""),
      raw = "RTI_Median_five-choice_movement_time",
      checks = list(precision_status = c("ok", "ok", ""))),
    RVP_A_prime = list(path = "results/RVP_A_prime_scores.csv",
      score = "normative_z", raw = "RVP_A_prime", source = "score_source",
      values = c("HC_SUBJECT_OOF", "NONHC_FULL_HC_MODEL", ""),
      checks = list(precision_status = c("ok", "ok", "not_scored"))),
    RVP_latency = list(path = "RVP_latency_scores.csv", score = "normative_z",
      source = "score_source", values = c("HC_OOF", "full_HC_new_subject", ""),
      raw = "RVP_Median_latency",
      checks = list(precision_status = c("ok", "ok", ""))),
    SWM_Between_errors = list(path = "private/SWM_scores.csv", score = "normative_z",
      raw = "SWM_Between_errors", source = "normative_score_source",
      values = c("HC_subject_held_out", "nonHC_full_HC", ""),
      checks = list(posterior_precision_ok = c("TRUE", "TRUE", ""),
                    numerical_ok = c("TRUE", "TRUE", ""))),
    SWM_Strategy = list(path = "SWM_Strategy_scores.csv", score = "normative_z",
      raw = "SWM_Strategy", source = "normative_score_source",
      values = c("participant_OOF", "full_HC_new_person", "excluded"),
      checks = list(precision_status = c("passed", "passed", "not_scored")))
  )
  fails <- function(expr, text) {
    err <- tryCatch({force(expr); NULL}, error = identity)
    stopifnot(inherits(err, "error"), grepl(text, conditionMessage(err), fixed = TRUE))
  }
  for (outcome in names(cases)) {
    case <- cases[[outcome]]; out <- file.path(work, outcome)
    dir.create(out)
    input <- file.path(out, "analysis_input.csv")
    write.csv(raw, input, row.names = FALSE, na = "")
    native <- raw
    native[[case$score]] <- z
    if (outcome == "RTI_reaction") native$NM_raw <- raw[[case$raw]]
    if (!is.null(case$source)) native[[case$source]] <- case$values
    for (name in names(case$checks)) native[[name]] <- case$checks[[name]]
    if (outcome %in% c("RTI_movement", "RVP_latency")) {
      native$normative_z_full_reference <- c("9.99", z[2], "")
      native$normative_z_cross_validation <- c(z[1], "", "")
    }
    path <- file.path(out, case$path)
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    save_native <- function(x) write.csv(x, path, row.names = FALSE, na = "")
    save_native(native)
    exported <- cantab_export_results(outcome, input, out)
    result <- read.csv(exported, colClasses = "character", check.names = FALSE, na.strings = NULL)
    stopifnot(identical(result$Subject_ID, raw$Subject_ID),
      identical(result$normative_z, z),
      identical(result$score_status, c("available", "available", "not_available")),
      identical(result$reference_model, c("cross_validation", "all_healthy_controls", "")),
      identical(as.numeric(result$raw_score[1:2]), raw[[case$raw]][1:2]),
      all(result$outcome == outcome),
      all(result$score_direction == if (outcome == "RVP_A_prime") "higher_is_better" else "higher_is_worse"),
      identical(names(result), c("Subject_ID", "Diagnosis", "Timepoint", "Platform", "Gender",
        "FSIQ", "Age_baseline_years", "Time_since_baseline_years", "outcome", "raw_score",
        "normative_z", "score_status", "reference_model", "score_direction")),
      identical(exported, cantab_export_results(outcome, input, out)),
      file.exists(file.path(out, "analysis_summary.txt")),
      file.exists(file.path(out, "score_dictionary.csv")))
    changed <- native[c(2, 1, 3), ]; save_native(changed)
    fails(cantab_export_results(outcome, input, out), "identity")
    changed <- native; changed[[case$score]][1] <- "Inf"; save_native(changed)
    fails(cantab_export_results(outcome, input, out), "finite")
    changed <- native; changed[[names(case$checks)[1]]][1] <- "failed"; save_native(changed)
    fails(cantab_export_results(outcome, input, out), "precision")
    changed <- native; changed[[case$source]][1] <- case$values[2]; save_native(changed)
    fails(cantab_export_results(outcome, input, out), "reference")
    if (outcome %in% c("RTI_movement", "RVP_latency")) {
      changed <- native; changed$normative_z_cross_validation[1] <- "9.99"; save_native(changed)
      fails(cantab_export_results(outcome, input, out), "reference")
    }
    save_native(native)
    writeLines("edited externally", exported)
    fails(cantab_export_results(outcome, input, out), "changed")
  }
  cat("PASS: six common score schemas, exact score text/IDs, reference routing, precision, row identity and overwrite protection.\n")
})
