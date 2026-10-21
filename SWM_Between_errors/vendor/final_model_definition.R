# Final normative model: hurdle discretized normal counts with exposure offsets.
# mu models positive-count location, sigma its scale, and hu the probability of zero.
swm_final_model <- function(data) {
    main <- swm_mean_formula()
    sigma <- sigma ~ Platform + Practice_bin + Age_baseline_10y + offset(log_trial_exposure)
    hu <- hu ~ Platform + Practice_bin + Age_baseline_10y + Connect_followup
    environment(main) <- environment(sigma) <- environment(hu) <- environment(swm_final_model)
    formula <- bf(main, sigma, hu)
    family <- hdg_family()
    priors <- c(swm_common_priors(formula, data, family),
        set_prior(paste0("normal(", log(5), ",0.7)"), class = "Intercept", dpar = "sigma"),
        set_prior("normal(0,0.3)", class = "b", dpar = "sigma"),
        set_prior("normal(0,0.7)", class = "b", dpar = "sigma", coef = "PlatformConnect"),
        set_prior("normal(0,0.4)", class = "b", dpar = "sigma", coef = "Age_baseline_10y")
    )
    list(id = "normative", formula = formula, priors = priors,
         family = family, stanvars = hdg_stanvars())
}
