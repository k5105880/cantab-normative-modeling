# Preparing CANTAB_connect_RS6_rawscores.csv

The input is a wide table with **one row per participant and assessment occasion**. Repeated assessments have additional rows with the same `Subject_ID`. The pair `Subject_ID,Timepoint` must be unique, even when platforms differ. Use actual assessment order. The conventional filename may be changed when passing `input_file`.

The examples use the [documented folder layout](USAGE.md): working directory `project/working_directory/`, input files in `project/input/`.

Save tables as **comma-separated UTF-8 CSV**; a UTF-8 byte-order mark (BOM) is accepted. Include one header for every column and do not export R row names. Quoted commas and line breaks inside cells are supported, but unclosed quotes, blank records and inconsistent field counts cause an error.

## Source records and column mapping

Combine CANTAB task exports with a participant master table. Export column names vary with software and export settings: explicitly map the relevant raw measures to the names below. Do not substitute means for medians or other reaction-time conditions for five-choice scores. Read identifiers as character strings so that leading zeros remain intact. A person must have the same pseudonymous ID across platforms.

Prepare two CSVs using [participants_header.csv](../data/templates/participants_header.csv) and [assessments_header.csv](../data/templates/assessments_header.csv). These are header-only templates.

### Participant master: one row per person

| Column | Definition |
|---|---|
| `Subject_ID` | Unique, nonempty pseudonymous text ID; no surrounding spaces |
| `Diagnosis` | `HC` healthy control, `MDD` major depressive disorder, `BD` bipolar disorder; use the resolved participant labels |
| `Gender` | `F` or `M`, the codes supported by this implementation; blank if unavailable |
| `FSIQ` | Baseline estimated full-scale IQ on the mean-100, SD-15 scale, e.g. `105`; for premorbid IQ, use NART or a validated language adaptation such as JART. Enter the IQ estimate, not a reading-test raw score, and keep it constant across visits |
| `Age_baseline_years` | Baseline age in years, e.g. 42.5; constant across visits |
| `Baseline_date` | Baseline assessment date, strictly `YYYY-MM-DD` |

`FSIQ` names the numerical IQ covariate, not a particular test. NART (National Adult Reading Test) is an English reading test; JART (Japanese Adult Reading Test) is its Japanese adaptation. Use an instrument and IQ-conversion equation validated for the participant's language and relevant population. Record the instrument, language/version and conversion method with your dataset. Apply a consistent definition within an analysis; the software does not convert reading errors to IQ or automatically harmonize estimates from different instruments. The model transformation remains `(FSIQ - 100) / 15`.

NART and JART both have published equations for estimating full-scale IQ: [NART validation](https://pubmed.ncbi.nlm.nih.gov/27624393/) and [JART development](https://onlinelibrary.wiley.com/doi/10.1111/j.1440-1819.2006.01510.x).

### Assessment table: one row per person/occasion

| Column | Definition and units |
|---|---|
| `Subject_ID` | Exact ID from the master table |
| `Assessment_date` | Assessment date, `YYYY-MM-DD` |
| `Timepoint` | Positive integer, baseline = 1; subsequent occasions 2, 3, 4, ... |
| `Platform` | Exactly `RS6` or `Connect` |
| `SWM_Between_errors` | Raw nonnegative integer count; do not divide RS6 by four |
| `SWM_Strategy` | Raw integer score: RS6 8–56, Connect 2–14; the model performs the RS6 division by four |
| `RVP_A_prime` | Raw A′ in [0,1], including 0 and 1; not a percentage |
| `RVP_Median_latency` | Positive median latency, milliseconds |
| `RTI_Median_five-choice_reaction_time` | Positive five-choice median reaction time, milliseconds |
| `RTI_Median_five-choice_movement_time` | Positive five-choice median movement time, milliseconds |

Keep the hyphen in `five-choice`. When reading or constructing an R data frame use `check.names = FALSE`. Convert seconds to milliseconds (`*1000`) only if the source units are verified as seconds. Do not log-transform, standardize or platform-correct the raw scores. Fourth and later occasions are grouped internally by the models.

The Between-errors model fixes relative trial exposure at **RS6 = 4 and Connect = 1**. Confirm that this 4:1 ratio matches the task protocol represented by your exported raw counts. It is a model assumption, not a value inferred from the CSV. Keep RS6 counts undivided; a different protocol requires review of the [model settings](SETTINGS.md).

## Missing values and joining exports

Use blank or `NA` (also accepted: `NaN`) for missing cells. Translate proprietary missing-value codes only according to the export specification. Never use 0, 999 or −9 as generic missing codes. Zero is valid for Between errors and A′, but invalid for the three time measures. A missing outcome does not require deleting other valid outcomes on that row; eligibility is outcome-specific. If you analyze only one outcome, keep all six outcome column headers and leave the unused outcome values blank. Strategy scores are checked against the supported range. Missing required outcomes or covariates are not imputed: an otherwise structurally valid row with missing `Gender`, `FSIQ`, baseline age or elapsed years is excluded from that outcome's model training and scoring. It can remain in the result with `score_status = "not_available"`. Some structural checks still require values, including valid dates for the CSV builder and elapsed time zero for a baseline row; fix a reported structural error before fitting.

Map export columns explicitly. The names on the right below are illustrative placeholders to replace with your actual source names:

```r
source("run_analysis.R")
export <- cantab_read_characters("../input/export.csv")
column_map <- c(
  Subject_ID = "participant_code", Assessment_date = "test_date",
  Timepoint = "visit_number", Platform = "test_platform",
  SWM_Between_errors = "swm_between_errors_raw", SWM_Strategy = "swm_strategy_raw",
  RVP_A_prime = "rvp_aprime_raw", RVP_Median_latency = "rvp_latency_median_ms",
  "RTI_Median_five-choice_reaction_time" = "rti_five_choice_reaction_median_ms",
  "RTI_Median_five-choice_movement_time" = "rti_five_choice_movement_median_ms"
)
stopifnot(all(unname(column_map) %in% names(export)))
assessments <- export[, unname(column_map), drop = FALSE]
names(assessments) <- names(column_map)
write.csv(assessments, "../input/assessments.csv", row.names = FALSE, na = "",
          fileEncoding = "UTF-8")
```

`cantab_read_characters()` checks CSV structure before reading and keeps all columns as text, including participant IDs and decimal score strings. It accepts export-specific column names; the mapping above supplies the required names. Empty cells, `NA` and `NaN` become missing values.

If tasks arrive in separate tables, normalize ID, date, platform and occasion first; then join on those keys, never on row position. Keep missing tasks with a full join and verify the keys before and after merging:

```r
keys <- c("Subject_ID", "Timepoint", "Assessment_date", "Platform")
stopifnot(!anyDuplicated(swm[keys]), !anyDuplicated(rvp[keys]), !anyDuplicated(rti[keys]))
assessments <- Reduce(function(x, y) merge(x, y, by = keys, all = TRUE, sort = FALSE),
                      list(swm, rvp, rti))
```

Here `swm`, `rvp`, and `rti` are already harmonized data frames containing the key columns and their respective outcome columns. Check row counts, duplicate assessments and mismatched dates/platforms. After finalizing the input, keep its values and row order unchanged when resuming a run.

## Build and validate the final CSV

```r
source("run_analysis.R")
create_cantab_input(
  participants = "../input/participants.csv",
  assessments = "../input/assessments.csv",
  output_file = "../input/CANTAB_connect_RS6_rawscores.csv"
)
input <- "../input/CANTAB_connect_RS6_rawscores.csv"
x <- read_cantab_input(input)
nrow(x)
length(unique(x$Subject_ID))
table(x$Diagnosis, x$Platform)
colSums(is.na(x[, cantab_required_columns]))
run_cantab("RTI_reaction", input, "../output/reaction_check", mode = "check")
```

The builder matches by ID, preserves assessment row order and raw score values, and computes:

- `Time_since_baseline_years = (Assessment_date - Baseline_date) / 365.25`.
- `Age_years = Age_baseline_years + Time_since_baseline_years`.

Dates themselves are omitted from the output. The output contains the 14 required columns in [primary_raw_header.csv](../data/templates/primary_raw_header.csv), plus `Age_years`. Baseline elapsed time must be zero; later visits must have positive, strictly increasing elapsed time. Occasion numbers are supplied explicitly. An existing output file is not overwritten unless `overwrite = TRUE` is specified.

If you already have validated baseline age and elapsed years, directly prepare the final-column template and validate it with `read_cantab_input()`; no date recalculation is necessary. Preserve source precision. Alternative day columns are accepted by the common reader only when consistent with years; using the explicit year columns is the portable contract for all six outcomes.

Age checks of 5–110 years and elapsed-time checks of 0–50 years detect likely unit errors; they do not define model validity limits. Structural validation does not guarantee sufficient HC participants, platform/occasion overlap, or identifiable designs in every training fold. The outcome `check` stage evaluates these model-specific requirements without MCMC. A′ specifically requires at least 20 participants with two or more eligible assessments in the full HC set and in each held-out training set. Having 20 such participants in the full CSV is not enough to establish that every fold meets the condition; see the [settings guide](SETTINGS.md).

[participants.csv](../data/example/participants.csv) and [assessments.csv](../data/example/assessments.csv) are wholly synthetic examples. `source("data/example/create_example.R")` creates a final CSV in a temporary directory. The example is intentionally small and is not a model-fitting dataset.

## Optional cohort and fold settings

A fold is one of five participant groups used for cross-validation; all visits from a person stay together. Normally omit `reference_dir`. All supplied participants are then considered, subject to outcome-specific missingness and score-range checks. Models are trained on `Diagnosis = HC`, and participant-grouped folds are generated automatically.

To select a predefined participant set, place `diagnosis_reference.csv` with columns `Subject_ID,Diagnosis` in a separate directory and pass that directory as `reference_dir`. IDs must be unique and diagnoses must match the input. Participants absent from the table are excluded at all visits, including from `normative_scores.csv`; labels are never rewritten. The run saves a source-row selection map in `COHORT_SELECTION.csv` and source/selected counts in `INPUT_AUDIT_R.json`. The two-column CSV is sufficient; no JSON file needs to be created. If an existing `diagnosis_reference.json` is supplied, its saved identity and hashes are also checked.

RTI reaction and RVP A′ also accept a supplied fold table within that directory:

| Outcome | Fold filename |
|---|---|
| RTI reaction | `frozen_subject_folds.csv` or `input/frozen_subject_folds.csv` |
| RVP A′ | `primary_subject_folds.csv` or `reference/primary_subject_folds.csv` |

Fold columns are `Subject_ID,CV_fold`, with values 1–5. Include exactly the eligible HC participants and keep all observations of a participant together. Use only one of the two listed locations for each outcome; multiple fold files cause an error rather than selecting one silently. If no fold table is present, a partition is generated. The other four outcomes use their own automatic fold procedures.

```r
run_cantab("RTI_reaction", input, "../output/reaction_selected_check",
           mode = "check", reference_dir = "../input/cohort_settings")
```

Next, follow the [execution, results and settings guide](USAGE.md).
