# Mean formula and shared priors for the SWM error-count normative model.
swm_mean_formula <- function() {
    terms <- c(
        "s(Age_baseline_10y, k = 5)",
        "s(Time_since_baseline_1y, k = 4, bs = 'tp', m = c(2, 0))",
        "Gender", "FSIQ_15", "Practice_bin", "Platform",
        "(1 + Time_since_baseline_1y | Subject_ID)",
        "Platform:Time_since_baseline_1y", "Practice_bin:Platform",
        "offset(log_trial_exposure)"
    )
    as.formula(paste("NM_response ~", paste(terms, collapse = " + ")), env = environment(swm_mean_formula))
}

swm_common_priors <- function(formula, data, family) {
    priors <- c(
        set_prior(paste0("normal(", log(5), ",0.7)"), class = "Intercept"),
        set_prior("normal(0,0.3)", class = "b"),
        set_prior("exponential(3.5)", class = "sd", group = "Subject_ID", coef = "Intercept"),
        set_prior("exponential(15)", class = "sd", group = "Subject_ID", coef = "Time_since_baseline_1y"),
        set_prior("exponential(3.5)", class = "sds"),
        set_prior("lkj(2)", class = "cor"),
        set_prior(paste0("normal(", qlogis(0.1), ",0.9)"), class = "Intercept", dpar = "hu"),
        set_prior("normal(0,0.3)", class = "b", dpar = "hu"),
        set_prior("normal(0,0.7)", class = "b", dpar = "hu", coef = "Age_baseline_10y")
    )
    available <- get_prior(formula, data = data, family = family)
    terms <- available$coef[available$class == "b" & available$dpar == "" & nzchar(available$coef)]
    for (term in terms) {
        sd <- if (term == "PlatformConnect") 1 else if (
            grepl("PlatformConnect", term) && grepl("Practice_bin", term)
        ) 0.5 else NA_real_
        if (is.finite(sd))
            priors <- c(priors, set_prior(paste0("normal(0,", sd, ")"), class = "b", coef = term))
    }
    priors <- c(priors,
        set_prior("normal(0,1.5)", class = "b", dpar = "hu", coef = "PlatformConnect"),
        set_prior("normal(0,0.8)", class = "b", dpar = "hu", coef = "Connect_followup")
    )
    age_smooth <- available$coef[available$class == "sds" & available$dpar == "" &
        grepl("Age_baseline_10y", available$coef, fixed = TRUE)]
    age_linear <- available$coef[available$class == "b" & available$dpar == "" &
        grepl("sAge_baseline_10y_", available$coef, fixed = TRUE)]
    stopifnot(length(age_smooth) == 1, length(age_linear) == 1)
    c(priors,
        set_prior("exponential(1)", class = "sds", coef = age_smooth),
        set_prior("normal(0,0.7)", class = "b", coef = age_linear)
    )
}
