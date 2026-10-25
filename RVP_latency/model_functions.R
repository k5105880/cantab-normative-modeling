# Normative model construction, grouped validation and predictive calibration.
family_support_check <- function(x, family_id, mode = c("fit", "prediction")) {
    if(!identical(family_id,"gamma_raw"))stop("Only the prespecified normative family is supported")
    mode <- match.arg(mode)
    x <- x[!is.na(x)]
    if (length(x) == 0) {
        return(list(valid = FALSE, reason = "no non-missing observations"))
    }
    ok <- all(x > 0)
    return(list(valid = ok, reason = if (ok) "ok" else "positive family requires x > 0 for fit and prediction"))
}

make_response <- function(x, family_id) {
    if(!identical(family_id,"gamma_raw"))stop("Only the prespecified normative family is supported")
    x
}

inverse_response_draws <- function(pp, family_id) {
    if(!identical(family_id,"gamma_raw"))stop("Only the prespecified normative family is supported")
    pp
}

brms_family_from_id <- function(family_id) {
    if(!identical(family_id,"gamma_raw"))stop("Only the prespecified normative family is supported")
    Gamma(link = "log")
}

prepare_model_data <- function(dat, spec, family_id, control_only = TRUE) {
    if(!identical(family_id,"gamma_raw"))stop("Only the prespecified normative family is supported")
    d <- dat
    if (control_only) 
        d <- d %>% filter(Diagnosis == control_label)
    support_mode <- if (control_only) 
        "fit"
    else "prediction"
    vars_needed <- c(spec$raw_col, "Age_baseline_10y", "Time_since_baseline_1y", "Gender", "FSIQ_15", 
        "Platform", "Practice_bin", "Subject_ID", "row_id_global", "Timepoint", "Time_since_baseline_years", 
        "Diagnosis")
    d <- d %>% filter(if_all(all_of(vars_needed), ~!is.na(.x))) %>% droplevels()
    support <- family_support_check(d[[spec$raw_col]], family_id, mode = support_mode)
    if (!support$valid) {
        stop("Family ", family_id, " is invalid for ", spec$outcome_prefix, ": ", support$reason)
    }
    d$NM_raw <- d[[spec$raw_col]]
    d$NM_response <- make_response(d[[spec$raw_col]], family_id)
    d
}

make_model_priors <- function(family_id, use_smooth, use_random_slope, aux_dpars, outcome_prefix = NULL) {
    if(!identical(family_id,"gamma_raw"))stop("Only the prespecified normative family is supported")
    b_sd_use <- rvp_latency_prior_b_sd
    aux_b_sd_use <- rvp_latency_prior_aux_b_sd
    sd_rate_use <- rvp_latency_prior_sd_rate
    sds_rate_use <- rvp_latency_prior_sds_rate
    priors <- c(set_prior(paste0("normal(0, ", b_sd_use, ")"), class = "b"), set_prior(paste0("exponential(", 
        sd_rate_use, ")"), class = "sd"))
    if (use_smooth) {
        priors <- c(priors, set_prior(paste0("exponential(", sds_rate_use, ")"), class = "sds"))
    }
    if (use_random_slope) {
        priors <- c(priors, set_prior(paste0("lkj(", prior_cor_eta, ")"), class = "cor"))
    }
    priors <- c(priors, set_prior(sprintf("normal(%.8f, %.4f)", log(rvp_latency_prior_intercept_mean_ms), 
        rvp_latency_prior_intercept_log_sd), class = "Intercept"), set_prior(sprintf("normal(%.8f, %.4f)", 
        log(rvp_latency_prior_shape_median), rvp_latency_prior_shape_log_sd), class = "Intercept", dpar = "shape"))
    for (dpar in aux_dpars) {
        priors <- c(priors, set_prior(paste0("normal(0, ", aux_b_sd_use, ")"), class = "b", dpar = dpar))
    }
    priors
}

build_model_definition <- function(model_dat, family_id, force_random_intercept = FALSE, outcome_prefix = NULL) {
    if(!identical(family_id,"gamma_raw"))stop("Only the prespecified normative family is supported")
    use_age_smooth <- identical(age_term_type, "smooth") && n_distinct(model_dat$Age_baseline_10y) >= 
        6
    use_time_smooth <- identical(time_term_type, "smooth") && n_distinct(model_dat$Time_since_baseline_1y) >= 
        5
    age_term <- if (use_age_smooth) {
        "s(Age_baseline_10y, k = 5)"
    }
    else {
        "Age_baseline_10y"
    }
    time_term <- cantab_time_term(use_time_smooth, include_platform_mean && nlevels(model_dat$Platform) >= 
        2 && include_platform_time_interaction)
    use_gender <- nlevels(model_dat$Gender) >= 2
    use_practice <- nlevels(model_dat$Practice_bin) >= 2
    use_platform <- include_platform_mean && nlevels(model_dat$Platform) >= 2
    use_fsiq <- is.finite(sd(model_dat$FSIQ_15)) && sd(model_dat$FSIQ_15) > 0
    subject_visit_counts <- model_dat %>% count(Subject_ID, name = "n_visits")
    n_subjects_2plus <- sum(subject_visit_counts$n_visits >= 2)
    use_random_slope <- !force_random_intercept && identical(random_effect_mode, "random_slope") && n_subjects_2plus >= 
        min_subjects_with_2plus_visits_for_random_slope
    random_term <- if (use_random_slope) {
        "(1 + Time_since_baseline_1y | Subject_ID)"
    }
    else {
        "(1 | Subject_ID)"
    }
    rhs_terms <- c(age_term, time_term)
    if (use_gender) 
        rhs_terms <- c(rhs_terms, "Gender")
    if (use_fsiq) 
        rhs_terms <- c(rhs_terms, "FSIQ_15")
    if (use_practice) 
        rhs_terms <- c(rhs_terms, "Practice_bin")
    if (use_platform) 
        rhs_terms <- c(rhs_terms, "Platform")
    if (use_platform && include_platform_time_interaction) {
        rhs_terms <- c(rhs_terms, "Platform:Time_since_baseline_1y")
    }
    if (use_platform && use_practice && include_platform_practice_interaction) {
        rhs_terms <- c(rhs_terms, "Platform:Practice_bin")
    }
    rhs_terms <- c(rhs_terms, random_term)
    main_formula <- as.formula(paste("NM_response ~", paste(rhs_terms, collapse = " + ")))
    scale_terms <- character(0)
    if (use_platform && include_scale_platform) 
        scale_terms <- c(scale_terms, "Platform")
    if (use_practice && include_scale_practice) 
        scale_terms <- c(scale_terms, "Practice_bin")
    scale_rhs <- if (length(scale_terms)) 
        paste(scale_terms, collapse = " + ")
    else "1"
    aux_formulas <- list()
    aux_dpars <- character(0)
    aux_b_dpars <- character(0)
    aux_formulas <- c(aux_formulas, list(as.formula(paste("shape ~", scale_rhs))))
    aux_dpars <- c(aux_dpars, "shape")
    if (scale_rhs != "1") 
        aux_b_dpars <- c(aux_b_dpars, "shape")
    model_formula <- do.call(brms::bf, c(list(main_formula), aux_formulas))
    family_obj <- brms_family_from_id(family_id)
    priors <- make_model_priors(family_id = family_id, use_smooth = use_age_smooth || use_time_smooth, use_random_slope = use_random_slope, 
        aux_dpars = unique(aux_b_dpars), outcome_prefix = outcome_prefix)
    list(formula = model_formula, family = family_obj, priors = priors, use_random_slope = use_random_slope, 
        model_label = if (use_random_slope) "random_slope" else "random_intercept", n_subjects_2plus = n_subjects_2plus, 
        use_age_smooth = use_age_smooth, use_time_smooth = use_time_smooth, use_platform = use_platform, 
        use_practice = use_practice, aux_dpars = unique(aux_dpars), aux_b_dpars = unique(aux_b_dpars))
}

is_stan_backend_or_compilation_error <- function(e) {
    msg <- tolower(conditionMessage(e))
    patterns <- c("compil", "toolchain", "cmdstan", "rstan", "stanheaders", "rtools", "mingw", "g\\+\\+", 
        "makefile", "dll", "connection", "processx", "executable")
    any(vapply(patterns, grepl, logical(1), x = msg, fixed = FALSE))
}

fit_normative_model <- function(model_dat, outcome_prefix, family_id, seed, cache_tag = NULL) {
    if (nrow(model_dat) < 20) {
        stop("Too few complete HC observations for ", outcome_prefix, ": ", nrow(model_dat))
    }
    def <- build_model_definition(model_dat, family_id, force_random_intercept = FALSE, outcome_prefix = outcome_prefix)
    cache_suffix <- if (is.null(cache_tag) || !nzchar(cache_tag)) 
        ""
    else paste0("_", cache_tag)
    cache_base <- paste0("brms_cache_", outcome_prefix, "_", family_id, "_", model_cache_version, "_", 
        def$model_label, cache_suffix)
    message("--------------------------------------------------")
    message("Fitting outcome: ", outcome_prefix)
    message("Family: ", family_id)
    message("Random effects: ", def$model_label)
    message("HC rows: ", nrow(model_dat), "; subjects >=2 visits: ", def$n_subjects_2plus)
    if (verbose_checks) {
        message("Formula:")
        print(def$formula)
    }
    fit <- tryCatch(cantab_brm(def$formula, data = model_dat, family = def$family, prior = def$priors, 
        chains = chains_fit, iter = iter_fit, warmup = warmup_fit, cores = cores_fit, seed = seed, init = 0, 
        backend = "cmdstanr", save_pars = save_pars(all = TRUE), refresh = 50, file = cache_base, file_refit = file_refit_mode, 
        control = list(adapt_delta = adapt_delta_fit, max_treedepth = max_treedepth_fit)), error = function(e) e)
    if (inherits(fit, "error")) {
        if (is_stan_backend_or_compilation_error(fit)) 
            stop(fit)
        if (def$use_random_slope && allow_random_slope_fallback) {
            message("Random-slope fit failed; retrying with random intercept only.")
            message("Original error: ", conditionMessage(fit))
            def <- build_model_definition(model_dat, family_id, force_random_intercept = TRUE, outcome_prefix = outcome_prefix)
            cache_base <- paste0("brms_cache_", outcome_prefix, "_", family_id, "_", model_cache_version, 
                "_random_intercept_retry", cache_suffix)
            fit <- cantab_brm(def$formula, data = model_dat, family = def$family, prior = def$priors, 
                chains = chains_fit, iter = iter_fit, warmup = warmup_fit, cores = cores_fit, seed = seed, 
                init = 0, backend = "cmdstanr", save_pars = save_pars(all = TRUE), refresh = 50, file = cache_base, 
                file_refit = file_refit_mode, control = list(adapt_delta = adapt_delta_fit, max_treedepth = max_treedepth_fit))
            def$model_label <- "random_intercept_retry"
        }
        else {
            stop(fit)
        }
    }
    dg <- rvp_add_bfmi(diagnose_brms_fit(fit, outcome_prefix, family_id, def$model_label), fit)
    write.csv(dg, paste0(cache_base, "_diagnostics.csv"), row.names = FALSE)
    if (!all(dg$validation_pass)) 
        stop("Full-model convergence gate failed: ", cache_base)
    list(fit = fit, definition = def)
}

diagnose_brms_fit <- function(fit, outcome_prefix, family_id, model_label) {
    np <- tryCatch(brms::nuts_params(fit), error = function(e) NULL)
    n_divergent <- NA_integer_
    n_max_treedepth <- NA_integer_
    if (!is.null(np)) {
        n_divergent <- sum(np$Parameter == "divergent__" & np$Value == 1, na.rm = TRUE)
        n_max_treedepth <- sum(np$Parameter == "treedepth__" & np$Value >= max_treedepth_fit, na.rm = TRUE)
    }
    draws <- tryCatch(posterior::as_draws_df(fit), error = function(e) NULL)
    max_rhat <- NA_real_
    min_ess_bulk <- NA_real_
    min_ess_tail <- NA_real_
    if (!is.null(draws)) {
        # Only the first Cholesky diagonal and upper triangle are deterministic.
        # Undefined diagnostics for any other parameter must fail the gate.
        variables <- posterior::variables(draws)
        cholesky <- grep("^L_[0-9]+\\[[0-9]+,[0-9]+\\]$", variables, value = TRUE)
        structural <- vapply(cholesky, function(variable) {
            ij <- as.integer(strsplit(gsub(".*\\[|\\]", "", variable), ",")[[1L]])
            ij[1L] < ij[2L] || all(ij == 1L)
        }, logical(1))
        if (any(structural))
            draws <- posterior::subset_draws(draws, variable = setdiff(variables, cholesky[structural]))
        conv <- tryCatch(posterior::summarise_draws(draws, "rhat", "ess_bulk", "ess_tail"), error = function(e) NULL)
        if (!is.null(conv) && nrow(conv) > 0L &&
            all(is.finite(as.matrix(conv[c("rhat", "ess_bulk", "ess_tail")])))) {
            max_rhat <- max(conv$rhat)
            min_ess_bulk <- min(conv$ess_bulk)
            min_ess_tail <- min(conv$ess_tail)
        }
    }
    data.frame(outcome = outcome_prefix, family_id = family_id, model_label = model_label, n_divergent = n_divergent, 
        n_max_treedepth = n_max_treedepth, max_rhat = max_rhat, min_ess_bulk = min_ess_bulk, min_ess_tail = min_ess_tail, 
        stringsAsFactors = FALSE)
}

clamp_probability <- function(p) {
    pmin(pmax(p, 1e-06), 1 - 1e-06)
}

percentile_from_predictive_matrix <- function(pp, y, family_id) {
    if (ncol(pp) != length(y)) 
        stop("pp/y dimension mismatch")
    p <- rep(NA_real_, length(y))
    for (j in seq_along(y)) {
        draws_j <- pp[, j]
        draws_j <- draws_j[is.finite(draws_j)]
        if (!length(draws_j) || !is.finite(y[j])) 
            next
        p[j] <- mean(draws_j <= y[j])
    }
    p
}

randomized_pit_from_predictive_matrix <- function(pp, y, family_id, seed) {
    if (ncol(pp) != length(y)) 
        stop("pp/y dimension mismatch")
    set.seed(seed)
    u <- rep(NA_real_, length(y))
    for (j in seq_along(y)) {
        draws_j <- pp[, j]
        draws_j <- draws_j[is.finite(draws_j)]
        if (!length(draws_j) || !is.finite(y[j])) 
            next
        lower <- mean(draws_j < y[j])
        mass <- mean(draws_j == y[j])
        if (FALSE || mass > 0) {
            u[j] <- lower + runif(1) * mass
        }
        else {
            u[j] <- mean(draws_j <= y[j])
        }
    }
    clamp_probability(u)
}

sample_draw_ids <- function(fit, n_draw, seed) {
    draws_df <- posterior::as_draws_df(fit)
    set.seed(seed)
    sample(seq_len(nrow(draws_df)), size = min(n_draw, nrow(draws_df)), replace = FALSE)
}

score_reference_population <- function(fit, pred_dat, family_id, chunk_size = state_cdf_settings$chunk_size, score_label = "full_population") {
    scored <- state_scoring_guard({
        scored <- state_predict_marginal(fit, pred_dat, family_id, chunk_size = chunk_size)
        state_write_diagnostics(scored, output_prefix, family_id, score_label)
        scored
    })
    return(scored$z)
}

sample_moments <- function(x) {
    x <- x[is.finite(x)]
    n <- length(x)
    if (n < 4) {
        return(c(mean = NA_real_, sd = NA_real_, skewness = NA_real_, excess_kurtosis = NA_real_))
    }
    m <- mean(x)
    s <- sd(x)
    if (!is.finite(s) || s == 0) {
        return(c(mean = m, sd = s, skewness = NA_real_, excess_kurtosis = NA_real_))
    }
    z <- (x - m)/s
    c(mean = m, sd = s, skewness = mean(z^3), excess_kurtosis = mean(z^4) - 3)
}

coverage_from_yrep <- function(yrep, y, probs = c(0.5, 0.8, 0.9, 0.95)) {
    out <- lapply(probs, function(level) {
        alpha <- (1 - level)/2
        lower <- apply(yrep, 2, quantile, probs = alpha, na.rm = TRUE)
        upper <- apply(yrep, 2, quantile, probs = 1 - alpha, na.rm = TRUE)
        data.frame(interval = level, coverage = mean(y >= lower & y <= upper, na.rm = TRUE), stringsAsFactors = FALSE)
    })
    bind_rows(out)
}

make_calibration_by_covariate <- function(oof_hc, calibration_z_col, outcome, family_id, calibration_method) {
    if (!nrow(oof_hc)) 
        return(data.frame())
    d <- oof_hc %>% mutate(Age_quartile = factor(dplyr::ntile(Age_baseline_years, 4), levels = 1:4, labels = paste0("Q", 
        1:4)), FSIQ_quartile = factor(dplyr::ntile(FSIQ, 4), levels = 1:4, labels = paste0("Q", 1:4)), 
        Time_bin = cut(Time_since_baseline_years, breaks = c(-Inf, 1e-10, 1, 3, Inf), labels = c("baseline", 
            "0_to_1y", "1_to_3y", "gt_3y"), include.lowest = TRUE, right = TRUE))
    group_vars <- c("Platform", "Practice_bin", "Gender", "Age_quartile", "FSIQ_quartile", "Time_bin")
    bind_rows(lapply(group_vars, function(g) {
        tmp <- d %>% filter(!is.na(.data[[g]]), is.finite(.data[[calibration_z_col]])) %>% group_by(.data[[g]]) %>% 
            summarise(n = n(), n_subjects = n_distinct(Subject_ID), mean_z = mean(.data[[calibration_z_col]], 
                na.rm = TRUE), sd_z = sd(.data[[calibration_z_col]], na.rm = TRUE), .groups = "drop")
        names(tmp)[1] <- "level"
        tmp %>% mutate(outcome = outcome, family_id = family_id, calibration_method = calibration_method, 
            covariate = g, level = as.character(level), .before = 1)
    }))
}

make_subject_cv_profile <- function(model_dat) {
    model_dat %>% mutate(Subject_ID_chr = as.character(Subject_ID), Platform_chr = as.character(Platform), 
        Gender_chr = as.character(Gender)) %>% group_by(Subject_ID_chr) %>% summarise(n_rows = n(), n_visits = n_distinct(Timepoint), 
        longitudinal = n_distinct(Timepoint) >= 2, Age_baseline_10y = median(Age_baseline_10y, na.rm = TRUE), 
        FSIQ_15 = median(FSIQ_15, na.rm = TRUE), Gender = first(Gender_chr), platform_pattern = paste(sort(unique(Platform_chr)), 
            collapse = "+"), .groups = "drop") %>% rename(Subject_ID = Subject_ID_chr)
}

score_grouped_fold_balance <- function(model_dat, folds) {
    profile <- make_subject_cv_profile(model_dat)
    subj_fold <- data.frame(Subject_ID = as.character(model_dat$Subject_ID), fold = as.integer(folds), 
        stringsAsFactors = FALSE) %>% distinct()
    if (anyDuplicated(subj_fold$Subject_ID)) {
        stop("Grouped CV leakage detected while scoring folds.")
    }
    profile <- left_join(profile, subj_fold, by = "Subject_ID")
    if (anyNA(profile$fold)) 
        stop("Fold balance scoring failed: missing subject fold.")
    numeric_dispersion <- function(x) {
        ok <- is.finite(x)
        if (sum(ok) < 2) 
            return(0)
        x <- x[ok]
        f <- profile$fold[ok]
        sdx <- sd(x)
        if (!is.finite(sdx) || sdx < 1e-08) 
            return(0)
        fold_means <- tapply(x, f, mean)
        mean(((fold_means - mean(x))/sdx)^2)
    }
    categorical_dispersion <- function(x) {
        x <- as.character(x)
        levs <- sort(unique(x[!is.na(x)]))
        if (length(levs) <= 1) 
            return(0)
        sum(vapply(levs, function(lv) {
            p0 <- mean(x == lv, na.rm = TRUE)
            pf <- tapply(x == lv, profile$fold, mean, na.rm = TRUE)
            denom <- max(p0 * (1 - p0), 0.05)
            mean((pf - p0)^2, na.rm = TRUE)/denom
        }, numeric(1)))
    }
    n_subj <- as.numeric(table(factor(profile$fold, levels = seq_len(max(profile$fold)))))
    rows_by_fold <- as.numeric(tapply(profile$n_rows, factor(profile$fold, levels = seq_len(max(profile$fold))), 
        sum))
    size_subject_score <- mean(((n_subj - mean(n_subj))/max(mean(n_subj), 1))^2)
    size_row_score <- mean(((rows_by_fold - mean(rows_by_fold))/max(mean(rows_by_fold), 1))^2)
    score <- 2 * size_subject_score + 1 * size_row_score + 1 * numeric_dispersion(profile$Age_baseline_10y) + 
        0.5 * numeric_dispersion(profile$n_visits) + 1 * categorical_dispersion(profile$Gender) + 1 * 
        categorical_dispersion(profile$platform_pattern) + 1 * categorical_dispersion(profile$longitudinal)
    if (cv_balance_include_fsiq) {
        score <- score + 0.5 * numeric_dispersion(profile$FSIQ_15)
    }
    as.numeric(score)
}

make_fold_balance_audit <- function(model_dat, folds) {
    profile <- make_subject_cv_profile(model_dat)
    subj_fold <- data.frame(Subject_ID = as.character(model_dat$Subject_ID), fold = as.integer(folds), 
        stringsAsFactors = FALSE) %>% distinct()
    profile <- left_join(profile, subj_fold, by = "Subject_ID")
    numeric_metrics <- c("Age_baseline_10y", "FSIQ_15", "n_visits")
    out <- list()
    z <- 1L
    for (metric in numeric_metrics) {
        overall <- mean(profile[[metric]], na.rm = TRUE)
        for (k in sort(unique(profile$fold))) {
            value <- mean(profile[[metric]][profile$fold == k], na.rm = TRUE)
            out[[z]] <- data.frame(fold = k, metric = paste0("mean_", metric), level = "", value = value, 
                overall_value = overall, abs_diff = abs(value - overall), stringsAsFactors = FALSE)
            z <- z + 1L
        }
    }
    categorical_metrics <- c("Gender", "platform_pattern", "longitudinal")
    for (metric in categorical_metrics) {
        vals <- as.character(profile[[metric]])
        for (lv in sort(unique(vals[!is.na(vals)]))) {
            overall <- mean(vals == lv, na.rm = TRUE)
            for (k in sort(unique(profile$fold))) {
                idx <- profile$fold == k
                value <- mean(vals[idx] == lv, na.rm = TRUE)
                out[[z]] <- data.frame(fold = k, metric = paste0("prop_", metric), level = lv, value = value, 
                  overall_value = overall, abs_diff = abs(value - overall), stringsAsFactors = FALSE)
                z <- z + 1L
            }
        }
    }
    n_subjects_all <- nrow(profile)
    n_rows_all <- sum(profile$n_rows)
    for (k in sort(unique(profile$fold))) {
        idx <- profile$fold == k
        out[[z]] <- data.frame(fold = k, metric = "n_subjects", level = "", value = sum(idx), overall_value = n_subjects_all/length(unique(profile$fold)), 
            abs_diff = abs(sum(idx) - n_subjects_all/length(unique(profile$fold))), stringsAsFactors = FALSE)
        z <- z + 1L
        out[[z]] <- data.frame(fold = k, metric = "n_rows", level = "", value = sum(profile$n_rows[idx]), 
            overall_value = n_rows_all/length(unique(profile$fold)), abs_diff = abs(sum(profile$n_rows[idx]) - 
                n_rows_all/length(unique(profile$fold))), stringsAsFactors = FALSE)
        z <- z + 1L
    }
    bind_rows(out)
}

make_grouped_folds <- function(model_dat, K, seed, balanced = cv_balance_folds, n_candidates = cv_balance_candidates) {
    if (n_distinct(model_dat$Subject_ID) < K) {
        stop("Fewer HC subjects than requested CV folds.")
    }
    if (!balanced || n_candidates <= 1) {
        set.seed(seed)
        best_folds <- loo::kfold_split_grouped(K = K, x = if(is.factor(model_dat$Subject_ID))droplevels(model_dat$Subject_ID)else model_dat$Subject_ID)
        best_score <- score_grouped_fold_balance(model_dat, best_folds)
        best_candidate <- 1L
    }
    else {
        best_score <- Inf
        best_folds <- NULL
        best_candidate <- NA_integer_
        candidate_scores <- rep(NA_real_, n_candidates)
        for (cand in seq_len(n_candidates)) {
            set.seed(seed + cand - 1L)
            f <- loo::kfold_split_grouped(K = K, x = if(is.factor(model_dat$Subject_ID))droplevels(model_dat$Subject_ID)else model_dat$Subject_ID)
            sc <- score_grouped_fold_balance(model_dat, f)
            candidate_scores[cand] <- sc
            if (is.finite(sc) && sc < best_score) {
                best_score <- sc
                best_folds <- f
                best_candidate <- cand
            }
        }
        if (is.null(best_folds)) 
            stop("Unable to generate a valid balanced grouped fold partition.")
    }
    subj_fold <- data.frame(Subject_ID = as.character(model_dat$Subject_ID), fold = best_folds, stringsAsFactors = FALSE) %>% 
        distinct()
    bad <- subj_fold %>% count(Subject_ID) %>% filter(n > 1)
    if (nrow(bad) > 0) 
        stop("Grouped CV leakage detected: a subject appears in multiple folds.")
    audit <- make_fold_balance_audit(model_dat, best_folds) %>% mutate(balance_score = best_score, selected_candidate = best_candidate, 
        n_candidates = if (balanced) 
            n_candidates
        else 1L, K = K, .before = 1)
    list(folds = as.integer(best_folds), balance_score = best_score, selected_candidate = best_candidate, 
        audit = audit)
}

run_grouped_kfold <- function(fit, folds, save_fits, seed) {
    kf <- cantab_kfold(fit, K = max(folds), folds = folds, group = "Subject_ID", joint = cv_joint_likelihood, 
        save_fits = TRUE, compare = FALSE, chains = cv_chains, iter = cv_iter, warmup = cv_warmup, cores = cv_cores, 
        seed = seed, refresh = 100, file = NULL, control = list(adapt_delta = adapt_delta_fit, max_treedepth = max_treedepth_fit))
    rvp_validate_kfold(kf, fit$data, "RVP", seed, diagnose_brms_fit, max_treedepth_fit)
    if (!save_fits) 
        kf$fits <- NULL
    kf
}

extract_kfold_summary <- function(kf, outcome_prefix, family_id, prespecified) {
    est <- tryCatch(kf$estimates, error = function(e) NULL)
    if (!is.null(est)) 
        est <- tryCatch(as.matrix(est), error = function(e) NULL)
    get_est_flexible <- function(row_candidates, col = "Estimate") {
        if (is.null(est) || is.null(rownames(est)) || is.null(colnames(est)) || !col %in% colnames(est)) 
            return(NA_real_)
        hit <- row_candidates[row_candidates %in% rownames(est)]
        if (!length(hit)) {
            patt <- paste(row_candidates, collapse = "|")
            hit <- rownames(est)[grepl(patt, rownames(est), ignore.case = TRUE)]
        }
        if (!length(hit)) 
            return(NA_real_)
        val <- suppressWarnings(as.numeric(est[hit[1], col]))
        if (length(val) == 1 && is.finite(val)) 
            val
        else NA_real_
    }
    elpd <- get_est_flexible(c("elpd_kfold", "elpd"), "Estimate")
    se_elpd <- get_est_flexible(c("elpd_kfold", "elpd"), "SE")
    p_kfold <- get_est_flexible(c("p_kfold", "p_eff", "p"), "Estimate")
    kfoldic <- get_est_flexible(c("kfoldic", "kfold_ic", "ic"), "Estimate")
    if (!is.finite(elpd)) {
        pw <- tryCatch(kf$pointwise, error = function(e) NULL)
        if (!is.null(pw)) 
            pw <- tryCatch(as.matrix(pw), error = function(e) NULL)
        if (!is.null(pw) && !is.null(colnames(pw))) {
            elpd_cols <- intersect(c("elpd_kfold", "elpd"), colnames(pw))
            if (!length(elpd_cols)) {
                elpd_cols <- colnames(pw)[grepl("^elpd", colnames(pw), ignore.case = TRUE)]
            }
            if (length(elpd_cols)) {
                v <- suppressWarnings(as.numeric(pw[, elpd_cols[1]]))
                v <- v[is.finite(v)]
                if (length(v)) {
                  elpd <- sum(v)
                  se_elpd <- if (length(v) > 1) 
                    sqrt(length(v) * stats::var(v))
                  else NA_real_
                }
            }
        }
    }
    if (!is.finite(kfoldic) && is.finite(elpd)) 
        kfoldic <- -2 * elpd
    status <- if (is.finite(elpd)) 
        "ok"
    else "elpd_unavailable"
    reason <- if (is.finite(elpd)) 
        ""
    else paste0("Could not extract K-fold ELPD; estimate rows=", paste(if (is.null(est)) 
        character(0)
    else rownames(est), collapse = "|"), "; pointwise cols=", paste(tryCatch(colnames(as.matrix(kf$pointwise)), 
        error = function(e) character(0)), collapse = "|"))
    data.frame(outcome = outcome_prefix, family_id = family_id, prespecified = prespecified, status = status, 
        reason = reason, model_label = NA_character_, n_divergent_fullfit = NA_real_, n_max_treedepth_fullfit = NA_real_, 
        max_rhat_fullfit = NA_real_, min_ess_bulk_fullfit = NA_real_, min_ess_tail_fullfit = NA_real_, 
        elpd_kfold = elpd, se_elpd_kfold = se_elpd, p_kfold = p_kfold, kfoldic = kfoldic, stringsAsFactors = FALSE)
}

score_cross_validation <- function(kf_obj, family_id, model_dat, rpit_seed, score_label = "oof_population") {
    set.seed(rpit_seed + 101L)
    kp <- brms::kfold_predict(kf_obj, method = "posterior_predict", sample_new_levels = "gaussian")
    idx_model <- rvp_oof_index(kp$y, model_dat$NM_response)
    global_ids <- model_dat$row_id_global[idx_model]
    p_mid <- percentile_from_predictive_matrix(kp$yrep, kp$y, family_id)
    z <- qnorm(clamp_probability(p_mid))
    state_diagnostics <- NULL
    state_diagnostics <- state_scoring_guard({
        scored <- state_predict_oof(kf_obj, model_dat, family_id)
        state_write_diagnostics(scored, output_prefix, family_id, score_label)
        scored
    })
    pos <- match(global_ids, state_diagnostics$row_id_global)
    if (anyNA(pos)) 
        stop("Analytic OOF state score alignment failed")
    p_mid <- state_diagnostics$p_lower[pos]
    z <- state_diagnostics$z[pos]
    pit_randomized <- randomized_pit_from_predictive_matrix(kp$yrep, kp$y, family_id, seed = rpit_seed)
    z_rpit <- qnorm(clamp_probability(pit_randomized))
    list(row_id_global = global_ids, z = z, percentile_mid = p_mid, pit_randomized = pit_randomized, 
        z_rpit = z_rpit, state_diagnostics = state_diagnostics, y = kp$y, yrep = kp$yrep)
}

run_prior_predictive_check <- function(model_dat, definition, spec, selected_family, outcome_index) {
    cache_file <- paste0("brms_prior_cache_", spec$outcome_prefix, "_", selected_family, "_", analysis_version, 
        "_", definition$model_label)
    prior_fit <- tryCatch(cantab_brm(definition$formula, data = model_dat, family = definition$family, 
        prior = definition$priors, sample_prior = "only", chains = prior_pred_chains, iter = prior_pred_iter, 
        warmup = prior_pred_warmup, cores = min(cores_fit, prior_pred_chains), seed = seed_main + outcome_index * 
            4e+05 + 101, init = 0, backend = "cmdstanr", save_pars = save_pars(all = TRUE), refresh = 0, 
        file = cache_file, file_refit = file_refit_mode, control = list(adapt_delta = adapt_delta_fit, 
            max_treedepth = max_treedepth_fit)), error = function(e) e)
    if (inherits(prior_fit, "error")) {
        return(data.frame(outcome = spec$outcome_prefix, family_id = selected_family, status = "failed", 
            reason = conditionMessage(prior_fit), stringsAsFactors = FALSE))
    }
    draw_ids <- sample_draw_ids(prior_fit, prior_pred_draws, seed = seed_main + outcome_index * 4e+05 + 
        113)
    pp <- tryCatch(posterior_predict(prior_fit, draw_ids = draw_ids), error = function(e) e)
    if (inherits(pp, "error")) {
        return(data.frame(outcome = spec$outcome_prefix, family_id = selected_family, status = "prediction_failed", 
            reason = conditionMessage(pp), stringsAsFactors = FALSE))
    }
    pp_raw <- inverse_response_draws(pp, selected_family)
    vals <- as.numeric(pp_raw)
    vals <- vals[is.finite(vals)]
    obs <- model_dat$NM_raw[is.finite(model_dat$NM_raw)]
    if (!length(vals) || !length(obs)) {
        return(data.frame(outcome = spec$outcome_prefix, family_id = selected_family, status = "empty", 
            reason = "No finite prior predictive or observed values", stringsAsFactors = FALSE))
    }
    qs <- quantile(vals, probs = c(0.01, 0.05, 0.5, 0.95, 0.99), na.rm = TRUE)
    data.frame(outcome = spec$outcome_prefix, family_id = selected_family, status = "ok", reason = "", 
        n_prior_predictive_values = length(vals), observed_min = min(obs), observed_max = max(obs), prior_pred_q01 = unname(qs[1]), 
        prior_pred_q05 = unname(qs[2]), prior_pred_median = unname(qs[3]), prior_pred_q95 = unname(qs[4]), 
        prior_pred_q99 = unname(qs[5]), prior_pred_prop_below_observed_min = mean(vals < min(obs)), prior_pred_prop_above_observed_max = mean(vals > 
            max(obs)), stringsAsFactors = FALSE)
}
