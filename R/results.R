# Common user-facing results. Internal model tables remain available for auditing.
cantab_result_spec <- function(outcome) {
  specs <- list(
    RTI_reaction = list(file = "analysis/RTI_analysis_ready.csv", score = "normative_z",
      raw = "RTI_Median_five-choice_reaction_time", native_raw = "NM_raw",
      model = "Student-t model for log-transformed reaction time with platform-specific participant intercepts",
      source = "normative_source", sources = c("HC_subject_OOF", "patient_full_HC_new_subject"),
      checks = c(normative_usable = "TRUE")),
    RTI_movement = list(file = "RTI_movement_scores.csv",
      score = "normative_z",
      raw = "RTI_Median_five-choice_movement_time",
      model = "Student-t model for log-transformed movement time",
      source = "score_source", sources = c("HC_OOF", "full_HC_new_subject"),
      reference_full = "normative_z_full_reference", reference_cv = "normative_z_cross_validation",
      checks = c(precision_status = "ok")),
    RVP_A_prime = list(file = "results/RVP_A_prime_scores.csv",
      score = "normative_z", raw = "RVP_A_prime",
      model = "Extended-support beta model for attention sensitivity",
      source = "score_source", sources = c("HC_SUBJECT_OOF", "NONHC_FULL_HC_MODEL"),
      checks = c(precision_status = "ok")),
    RVP_latency = list(file = "RVP_latency_scores.csv",
      score = "normative_z", raw = "RVP_Median_latency",
      model = "Gamma model for response latency",
      source = "score_source", sources = c("HC_OOF", "full_HC_new_subject"),
      reference_full = "normative_z_full_reference", reference_cv = "normative_z_cross_validation",
      checks = c(precision_status = "ok")),
    SWM_Between_errors = list(file = "private/SWM_scores.csv",
      score = "normative_z", raw = "SWM_Between_errors",
      model = "Hurdle discretized normal model for error counts with trial-count adjustment",
      source = "normative_score_source", sources = c("HC_subject_held_out", "nonHC_full_HC"),
      checks = c(posterior_precision_ok = "TRUE", numerical_ok = "TRUE")),
    SWM_Strategy = list(file = "SWM_Strategy_scores.csv",
      score = "normative_z", raw = "SWM_Strategy",
      model = "Minimum-inflated discretized beta model for strategy scores",
      source = "normative_score_source", sources = c("participant_OOF", "full_HC_new_person"),
      checks = c(precision_status = "passed"))
  )
  if (length(outcome) != 1L || !outcome %in% names(specs)) stop("Unknown outcome.")
  specs[[outcome]]
}

cantab_export_results <- function(outcome, input_file, output_dir) {
  spec <- cantab_result_spec(outcome)
  output_dir <- normalizePath(output_dir, winslash = "/", mustWork = TRUE)
  input <- read_cantab_input(input_file)
  raw <- attr(input, "raw_characters")
  path <- file.path(output_dir, spec$file)
  if (length(path) != 1L || !file.exists(path)) stop("Expected one completed score table for ", outcome)
  native <- cantab_read_characters(path)
  native_raw <- if (is.null(spec$native_raw)) spec$raw else spec$native_raw
  required <- c("Subject_ID", "Diagnosis", "Platform", "Timepoint", spec$score,
                native_raw, names(spec$checks), spec$source, spec$reference_full, spec$reference_cv)
  if (!all(required %in% names(native))) stop("Completed score table is missing required fields.")
  numeric_cells <- function(x) {
    n <- suppressWarnings(as.numeric(x))
    if (any(!is.na(x) & (!is.finite(n) | is.na(n))))
      stop("Scores and numeric fields must be finite numbers or missing.")
    n
  }
  if (nrow(native) != nrow(input) ||
      any(!vapply(c("Subject_ID", "Diagnosis", "Platform"),
        function(nm) identical(native[[nm]], raw[[nm]]), logical(1))) ||
      !identical(numeric_cells(native$Timepoint), input$Timepoint))
    stop("Score/input row identity differs.")
  if (!identical(numeric_cells(native[[native_raw]]), input[[spec$raw]]))
    stop("Score/input raw value identity differs.")
  score <- numeric_cells(native[[spec$score]])
  available <- is.finite(score)
  for (nm in names(spec$checks))
    if (anyNA(native[[nm]][available]) ||
        any(native[[nm]][available] != spec$checks[[nm]]))
      stop("Score precision checks have not passed.")
  healthy <- input$Diagnosis == "HC"
  expected <- ifelse(healthy, spec$sources[1], spec$sources[2])
  if (anyNA(native[[spec$source]][available]) ||
      any(native[[spec$source]][available] != expected[available]))
    stop("Incorrect score reference model.")
  if (!is.null(spec$reference_full)) {
    expected <- numeric_cells(native[[spec$reference_full]])
    expected[healthy] <- numeric_cells(native[[spec$reference_cv]])[healthy]
    if (!identical(score, expected)) stop("Incorrect score reference model.")
  }
  columns <- c("Subject_ID", "Diagnosis", "Timepoint", "Platform", "Gender", "FSIQ",
               "Age_baseline_years", "Time_since_baseline_years")
  result <- input[columns]
  attr(result, "raw_characters") <- NULL
  # Preserve the original score text, including all digits and missing rows.
  result$outcome <- outcome
  result$raw_score <- raw[[spec$raw]]
  result$normative_z <- native[[spec$score]]
  result$score_status <- ifelse(available, "available", "not_available")
  result$reference_model <- ifelse(available,
    ifelse(healthy, "cross_validation", "all_healthy_controls"), "")
  result$score_direction <- if (outcome == "RVP_A_prime") "higher_is_better" else "higher_is_worse"
  meanings <- c(
    Subject_ID = "Participant identifier; import as text to preserve leading zeros.",
    Diagnosis = "HC: healthy control; MDD: major depressive disorder; BD: bipolar disorder.",
    Timepoint = "Assessment number: 1 is baseline, then 2, 3 and so on.",
    Platform = "CANTAB platform: RS6 or Connect.",
    Gender = "Input sex/gender code: F or M.",
    FSIQ = "Baseline estimated full-scale IQ (mean 100, SD 15), e.g. from NART or a validated language adaptation such as JART; use the IQ estimate, not the reading-test raw score. Fixed across visits.",
    Age_baseline_years = "Age at baseline in years.",
    Time_since_baseline_years = "Elapsed time since baseline in years.",
    outcome = "Outcome folder name identifying the fitted measure.",
    raw_score = "Original outcome value; time measures are in milliseconds and A prime is on 0 to 1.",
    normative_z = "Predictive cumulative probability expressed on the standard normal scale. Zero corresponds to probability 0.5. Discrete values and boundary masses use the midpoint of their probability interval, so a raw median need not score zero.",
    score_status = "available: a score is present; not_available: the row was retained without a usable score.",
    reference_model = "cross_validation: reference model excludes this healthy participant; all_healthy_controls: model uses all eligible healthy controls.",
    score_direction = "higher_is_better for A prime; higher_is_worse for the other five measures. Signs are not reversed."
  )
  dictionary <- data.frame(column = names(meanings), meaning = unname(meanings))
  summary <- c(
    "CANTAB normative modeling results",
    paste("Outcome:", outcome), paste("Model:", spec$model),
    paste("Assessment rows:", nrow(result)), paste("Available scores:", sum(available)),
    paste("Rows without a score:", sum(!available)),
    paste("Direction:", result$score_direction[1]),
    "Open normative_scores.csv for the result table and score_dictionary.csv for column definitions.",
    "The score is for one assessment relative to its model-predicted reference distribution.",
    "It is not a change score and it is not an individual diagnosis.",
    "Healthy participants are scored with participant-grouped five-fold cross-validation.",
    "Other participants are scored using the model fitted to all eligible healthy controls.",
    "Selected analysis rows, including rows without scores, remain in their input order. Normative scores are copied without rounding.",
    "Other files in this run directory are model objects, intermediate results and diagnostics."
  )
  staging <- tempfile("result_export_", tmpdir = output_dir)
  dir.create(staging)
  on.exit(unlink(staging, recursive = TRUE), add = TRUE)
  files <- c("normative_scores.csv", "score_dictionary.csv", "analysis_summary.txt")
  write.csv(result, file.path(staging, files[1]), row.names = FALSE, na = "", fileEncoding = "UTF-8")
  write.csv(dictionary, file.path(staging, files[2]), row.names = FALSE, fileEncoding = "UTF-8")
  writeLines(summary, file.path(staging, files[3]), useBytes = TRUE)
  bytes <- function(p) readBin(p, "raw", n = file.info(p)$size)
  for (name in files) {
    target <- file.path(output_dir, name)
    if (file.exists(target) && !identical(bytes(target), bytes(file.path(staging, name))))
      stop("Exported result has changed: ", target, ". Use a new output directory.")
  }
  for (name in files) {
    target <- file.path(output_dir, name)
    if (!file.exists(target) && !file.copy(file.path(staging, name), target))
      stop("Could not write result: ", target)
  }
  target <- file.path(output_dir, files[1])
  message("Normative scores: ", target)
  invisible(target)
}
