# Public diagnosis labels must select the correct reference model and remain unchanged.
local({
  source("R/input_data.R", local = TRUE)
  source("R/results.R", local = TRUE)
  work <- tempfile("diagnosis_labels_")
  dir.create(work)
  on.exit(unlink(work, recursive = TRUE), add = TRUE)
  labels <- c("HC", "MDD", "BD")
  participants <- data.frame(Subject_ID = c("001", "002", "003"), Diagnosis = labels,
    Gender = c("F", "M", "F"), FSIQ = c(110, 95, 100), Age_baseline_years = c(40, 52, 45),
    Baseline_date = "2020-01-01")
  assessments <- data.frame(Subject_ID = participants$Subject_ID,
    Assessment_date = "2020-01-01", Timepoint = 1, Platform = "Connect",
    SWM_Between_errors = c(0, 3, 2), SWM_Strategy = c(2, 4, 3),
    RVP_A_prime = c(.95, .8, .85), RVP_Median_latency = c(400, 450, 420),
    check.names = FALSE)
  assessments[["RTI_Median_five-choice_reaction_time"]] <- c(300, 350, 320)
  assessments[["RTI_Median_five-choice_movement_time"]] <- c(200, 250, 220)
  input <- file.path(work, "input.csv")
  create_cantab_input(participants, assessments, input)
  dat <- read_cantab_input(input)
  stopifnot(identical(dat$Diagnosis, labels))
  refs <- file.path(work, "reference"); dir.create(refs)
  write.csv(participants[c(1, 3), c("Subject_ID", "Diagnosis")],
    file.path(refs, "diagnosis_reference.csv"), row.names = FALSE)
  subset <- cantab_prepare_input(input, file.path(work, "subset"), refs, "RTI_reaction")
  stopifnot(identical(read_cantab_input(subset)$Diagnosis, c("HC", "BD")))
  for (bad_label in c("unknown", "hc", "H", "U", "B", " H ")) {
    bad <- participants; bad$Diagnosis[1] <- bad_label
    error <- tryCatch(create_cantab_input(bad, assessments, file.path(work, "invalid.csv")),
      error = identity)
    stopifnot(inherits(error, "error"), !file.exists(file.path(work, "invalid.csv")))
  }
  for (outcome in c("RTI_reaction", "RTI_movement", "RVP_A_prime", "RVP_latency",
                    "SWM_Between_errors", "SWM_Strategy")) {
    spec <- cantab_result_spec(outcome)
    out <- file.path(work, outcome)
    native_file <- file.path(out, spec$file)
    dir.create(dirname(native_file), recursive = TRUE)
    native <- dat[c("Subject_ID", "Diagnosis", "Platform", "Timepoint")]
    native_raw <- if (is.null(spec$native_raw)) spec$raw else spec$native_raw
    native[[native_raw]] <- dat[[spec$raw]]
    native[[spec$score]] <- c(-.2, .3, .4)
    native[[spec$source]] <- spec$sources[c(1, 2, 2)]
    for (field in names(spec$checks)) native[[field]] <- spec$checks[[field]]
    if (!is.null(spec$reference_full)) {
      native[[spec$reference_full]] <- c(9, .3, .4)
      native[[spec$reference_cv]] <- c(-.2, 9, 9)
    }
    write.csv(native, native_file, row.names = FALSE)
    path <- cantab_export_results(outcome, input, out)
    result <- cantab_read_characters(path)
    stopifnot(identical(result$Diagnosis, labels),
      identical(result$reference_model, c("cross_validation", "all_healthy_controls", "all_healthy_controls")),
      identical(result$score_status, rep("available", 3L)))
    dictionary <- cantab_read_characters(file.path(out, "score_dictionary.csv"))
    meaning <- dictionary$meaning[dictionary$column == "Diagnosis"]
    stopifnot(grepl("HC: healthy control", meaning, fixed = TRUE),
      grepl("MDD: major depressive disorder", meaning, fixed = TRUE),
      grepl("BD: bipolar disorder", meaning, fixed = TRUE))
  }
  cat("PASS: HC/MDD/BD input creation, unchanged cohort labels, invalid-label rejection, six exports and correct control/patient reference routing.\n")
})
