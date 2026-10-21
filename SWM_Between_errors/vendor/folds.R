# Outcome-blind participant-grouped cross-validation folds.
make_subject_cv_profile <- function(model_dat) {
    model_dat %>% mutate(Subject_ID_chr = as.character(Subject_ID), Platform_chr = as.character(Platform), Gender_chr = as.character(Gender)) %>% 
        group_by(Subject_ID_chr) %>% summarise(n_rows = n(), n_visits = n_distinct(Timepoint), longitudinal = n_distinct(Timepoint) >= 
        2, Age_baseline_10y = median(Age_baseline_10y, na.rm = TRUE), FSIQ_15 = median(FSIQ_15, na.rm = TRUE), 
        Gender = first(Gender_chr), platform_pattern = paste(sort(unique(Platform_chr)), collapse = "+"), .groups = "drop") %>% 
        rename(Subject_ID = Subject_ID_chr)
}

score_grouped_fold_balance <- function(model_dat, folds) {
    profile <- make_subject_cv_profile(model_dat)
    subj_fold <- data.frame(Subject_ID = as.character(model_dat$Subject_ID), fold = as.integer(folds), stringsAsFactors = FALSE) %>% 
        distinct()
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
    score <- 2 * size_subject_score + 1 * size_row_score + 1 * numeric_dispersion(profile$Age_baseline_10y) + 0.5 * 
        numeric_dispersion(profile$n_visits) + 1 * categorical_dispersion(profile$Gender) + 1 * categorical_dispersion(profile$platform_pattern) + 
        1 * categorical_dispersion(profile$longitudinal)
    as.numeric(score)
}

make_fold_balance_audit <- function(model_dat, folds) {
    profile <- make_subject_cv_profile(model_dat)
    subj_fold <- data.frame(Subject_ID = as.character(model_dat$Subject_ID), fold = as.integer(folds), stringsAsFactors = FALSE) %>% 
        distinct()
    profile <- left_join(profile, subj_fold, by = "Subject_ID")
    numeric_metrics <- c("Age_baseline_10y", "FSIQ_15", "n_visits")
    out <- list()
    z <- 1L
    for (metric in numeric_metrics) {
        overall <- mean(profile[[metric]], na.rm = TRUE)
        for (k in sort(unique(profile$fold))) {
            value <- mean(profile[[metric]][profile$fold == k], na.rm = TRUE)
            out[[z]] <- data.frame(fold = k, metric = paste0("mean_", metric), level = "", value = value, overall_value = overall, 
                abs_diff = abs(value - overall), stringsAsFactors = FALSE)
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
                out[[z]] <- data.frame(fold = k, metric = paste0("prop_", metric), level = lv, value = value, overall_value = overall, 
                  abs_diff = abs(value - overall), stringsAsFactors = FALSE)
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
        out[[z]] <- data.frame(fold = k, metric = "n_rows", level = "", value = sum(profile$n_rows[idx]), overall_value = n_rows_all/length(unique(profile$fold)), 
            abs_diff = abs(sum(profile$n_rows[idx]) - n_rows_all/length(unique(profile$fold))), stringsAsFactors = FALSE)
        z <- z + 1L
    }
    bind_rows(out)
}

make_grouped_folds <- function(model_dat, K, seed, balanced = TRUE, n_candidates = 500L) {
    # Retain participant order while removing unused factor levels before loo.
    groups <- if (is.factor(model_dat$Subject_ID)) droplevels(model_dat$Subject_ID) else model_dat$Subject_ID
    if (n_distinct(model_dat$Subject_ID) < K) {
        stop("Fewer HC subjects than requested CV folds.")
    }
    if (!balanced || n_candidates <= 1) {
        set.seed(seed)
        best_folds <- loo::kfold_split_grouped(K = K, x = groups)
        best_score <- score_grouped_fold_balance(model_dat, best_folds)
        best_candidate <- 1L
    }
    else {
        best_score <- Inf
        best_folds <- NULL
        best_candidate <- NA_integer_
        for (cand in seq_len(n_candidates)) {
            set.seed(seed + cand - 1L)
            f <- loo::kfold_split_grouped(K = K, x = groups)
            sc <- score_grouped_fold_balance(model_dat, f)
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
    list(folds = as.integer(best_folds), balance_score = best_score, selected_candidate = best_candidate, audit = audit)
}

