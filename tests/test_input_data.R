# Synthetic records only: tests for the public input preparation contract.
root <- normalizePath(if (file.exists("R/input_data.R")) "." else "..", mustWork = TRUE)
stopifnot("R input helper has been implemented" = file.exists(file.path(root, "R/input_data.R")))
source(file.path(root, "R/input_data.R"))
expect_error <- function(expr, pattern) {
  error <- tryCatch({ force(expr); NULL }, error = identity)
  stopifnot(inherits(error, "error"), grepl(pattern, conditionMessage(error)))
}
work <- tempfile("cantab_input_test_"); dir.create(work)
participants <- data.frame(Subject_ID = c("S001", "S002"), Diagnosis = c("HC", "MDD"),
  Gender = c("F", "M"), FSIQ = c(110, 95), Age_baseline_years = c(40, 52),
  Baseline_date = c("2020-01-01", "2021-02-03"))
assessments <- data.frame(Subject_ID = c("S002", "S001", "S001"),
  Assessment_date = c("2021-02-03", "2020-01-01", "2021-01-01"),
  Timepoint = c(1, 1, 2), Platform = c("Connect", "RS6", "Connect"),
  SWM_Between_errors = c(0, 4, NA), SWM_Strategy = c(2, 32, 8),
  RVP_A_prime = c(1, .9, .95), RVP_Median_latency = c(450, 400, 390),
  check.names = FALSE)
assessments[["RTI_Median_five-choice_reaction_time"]] <- c(350, 310, 300)
assessments[["RTI_Median_five-choice_movement_time"]] <- c(250, 220, 210)
raw <- file.path(work, "raw.csv")
create_cantab_input(participants, assessments, raw)
dat <- read_cantab_input(raw)
stopifnot(identical(dat$Subject_ID, assessments$Subject_ID),
          identical(dat$SWM_Strategy, c(2, 32, 8)), is.na(dat$SWM_Between_errors[3]),
          identical(dat$Time_since_baseline_years[1:2], c(0, 0)),
          abs(dat$Time_since_baseline_years[3] - 366 / 365.25) < 1e-12,
          abs(dat$Age_years[3] - 40 - dat$Time_since_baseline_years[3]) < 1e-12)
# Extra export metadata must neither be copied to the input schema nor affect validation.
metadata_names <- c("SWM_Problem_reached", "SWM Problem reached", "SWM_Problem_Reached",
  "SWM_ProblemReached", "SWM Problem Reached", "SWM_Problem_reached_QC",
  "HRSD17", "Onset_age", "illness_duration", "Export_note")
for (column in metadata_names) {
  extra <- assessments; extra[[column]] <- c("0", "15", "not recorded")
  extra_file <- file.path(work, "extra.csv")
  create_cantab_input(participants, extra, extra_file, overwrite = TRUE)
  stopifnot("Unrelated export fields must not enter the generated CSV" =
              identical(readLines(raw), readLines(extra_file)))
  extra <- cantab_read_characters(raw); extra[[column]] <- c("0", "15", "not recorded")
  write.csv(extra, extra_file, row.names = FALSE, na = "")
  checked <- read_cantab_input(extra_file)
  stopifnot(all(vapply(cantab_required_columns,
    function(nm) identical(checked[[nm]], dat[[nm]]), logical(1))))
}
expect_error(create_cantab_input(participants, assessments, raw), "exists")
bad <- participants; bad$Subject_ID[2] <- "S001"
expect_error(create_cantab_input(bad, assessments, file.path(work, "bad.csv")), "duplicate|Duplicate")
bad <- assessments; bad$Subject_ID[1] <- "UNKNOWN"
expect_error(create_cantab_input(participants, bad, file.path(work, "bad.csv")), "participant")
write_bad <- function(x) {
  p <- file.path(work, "bad.csv"); write.csv(x, p, row.names = FALSE, na = "")
  p
}
bad <- dat; bad$Timepoint[3] <- 1
expect_error(read_cantab_input(write_bad(bad)), "Duplicate")
bad <- dat; bad$Diagnosis[3] <- "BD"
expect_error(read_cantab_input(write_bad(bad)), "Diagnosis changes")
bad <- dat; bad$Time_since_baseline_years[3] <- 366
expect_error(read_cantab_input(write_bad(bad)), "years")
bad <- dat; bad[["RTI_Median_five-choice_reaction_time"]][1] <- 0
expect_error(read_cantab_input(write_bad(bad)), "positive")
bad <- dat; bad$FSIQ[3] <- 105
expect_error(read_cantab_input(write_bad(bad)), "FSIQ")
bad <- dat; bad$Platform[1] <- "connect"
expect_error(read_cantab_input(write_bad(bad)), "Platform")
bad <- dat; bad$Subject_ID[1] <- " S002"
expect_error(read_cantab_input(write_bad(bad)), "Subject_ID")
for (column in c("Diagnosis", "Platform", "Gender")) {
  bad <- dat; bad[[column]][1] <- paste0(bad[[column]][1], " ")
  expect_error(read_cantab_input(write_bad(bad)), "spaces")
}
bad <- dat; bad$RVP_A_prime[1] <- Inf
expect_error(read_cantab_input(write_bad(bad)), "finite")
out <- file.path(work, "out")
analysis <- cantab_prepare_input(raw, out, outcome = "RTI_reaction")
stopifnot(nrow(read_cantab_input(analysis)) == 3,
          identical(analysis, cantab_prepare_input(raw, out, outcome = "RTI_reaction")))
refs <- file.path(work, "refs"); dir.create(refs)
write.csv(participants[1, c("Subject_ID", "Diagnosis")],
          file.path(refs, "diagnosis_reference.csv"), row.names = FALSE)
sub <- cantab_prepare_input(raw, file.path(work, "subset"), refs, "RVP_A_prime")
stopifnot(identical(read_cantab_input(sub)$Subject_ID, c("S001", "S001")))
expect_error(cantab_prepare_input(raw, out, refs, "RTI_reaction"), "changed")
bad <- dat; bad$RVP_A_prime[1] <- .8
expect_error(cantab_prepare_input(write_bad(bad), out, outcome = "RTI_reaction"), "changed")
cat("PASS: CSV creation, row/score preservation, validation, reference cohort and audit identity\n")
