# Input conversion and fixed predictor scaling for this outcome.
prepare_normative_data <- function(infile) {
  data <- read.csv(
    infile,
    na.strings = c("", "NA", "NaN"),
    colClasses = c(Subject_ID = "character"),
    fileEncoding = "UTF-8-BOM", stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  data <- rvp_input_units(data)
  
  required_cols <- c(
    "Subject_ID", "Diagnosis", "Gender", "Platform", "Timepoint",
    "Age", "Age_baseline", "Time_since_baseline", "FSIQ",
    "RVP_Median_latency"
  )
  missing_cols <- setdiff(required_cols, names(data))
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }
  
  numeric_cols <- c(
    "Timepoint", "Age", "Age_baseline", "Time_since_baseline", "FSIQ",
    "RVP_Median_latency"
  )
  
  data <- data %>%
    mutate(
      row_id_global = row_number(),
      across(all_of(numeric_cols), ~ rvp_strict_numeric(.x, cur_column())),
      Platform = factor(Platform),
      Subject_ID = factor(Subject_ID),
      Diagnosis = factor(Diagnosis),
      Gender = factor(Gender),
      Timepoint = as.numeric(Timepoint),
      timepoint_f = factor(Timepoint),
      Practice_order = Timepoint - 1,
      Practice_bin = case_when(
        Practice_order == 0 ~ "0_baseline",
        Practice_order == 1 ~ "1_second",
        Practice_order == 2 ~ "2_third",
        Practice_order >= 3 ~ "3_fourth_or_later",
        TRUE ~ NA_character_
      ),
      Practice_bin = factor(
        Practice_bin,
        levels = c("0_baseline", "1_second", "2_third", "3_fourth_or_later")
      ),
      Age_years = Age / 365.25,
      Age_baseline_years = Age_baseline / 365.25,
      Time_since_baseline_years = Time_since_baseline / 365.25,
      Age_baseline_10y = (Age_baseline_years - reference_age_years) / age_scale_years,
      Time_since_baseline_1y = Time_since_baseline_years,
      FSIQ_15 = (FSIQ - reference_fsiq) / fsiq_scale_points
    )
  
  rvp_validate_input(data)
  
  control_label <- if ("HC" %in% levels(data$Diagnosis)) {
    "HC"
  } else {
    stop("Diagnosis must contain HC for the normative reference.")
  }
  message("Healthy-control label: ", control_label)
  
  if ("RS6" %in% levels(data$Platform)) {
    data$Platform <- relevel(data$Platform, ref = "RS6")
  }
  
  list(data=data,control_label=control_label)
}
