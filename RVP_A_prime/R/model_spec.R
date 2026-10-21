local({
  paths <- vapply(sys.frames(), function(f) if(is.null(f$ofile)) '' else f$ofile, character(1))
  source(file.path(dirname(normalizePath(tail(paths[nzchar(paths)],1))), '../cantab_design_audit.R'), local=parent.env(environment()))
})
# Normative RVP A-prime model specification.
build_model <- function(data) {
  required <- c("NM_response","Age_baseline_10y","Time_since_baseline_1y",
    "Gender","FSIQ_15","Practice_bin","Platform","PlatformConnect_indicator","Subject_ID")
  if (!all(required %in% names(data))) stop("Model data missing: ",paste(setdiff(required,names(data)),collapse=", "))
  if (any(!is.finite(data$NM_response)) || any(data$NM_response<0 | data$NM_response>1))
    stop("RVP A-prime responses must lie in [0,1]")
  if (!"xbeta" %in% getNamespaceExports("brms")) stop("brms with xbeta support is required")
  unique_n <- function(x) length(unique(x[!is.na(x)]))
  age <- if(unique_n(data$Age_baseline_10y)>=6) "s(Age_baseline_10y, k = 5)" else "Age_baseline_10y"
  time <- if(unique_n(data$Time_since_baseline_1y)>=5) "s(Time_since_baseline_1y, k = 4)" else "Time_since_baseline_1y"
  terms <- c(age,time)
  if(unique_n(data$Gender)>=2) terms <- c(terms,"Gender")
  if(is.finite(sd(data$FSIQ_15)) && sd(data$FSIQ_15)>0) terms <- c(terms,"FSIQ_15")
  practice <- unique_n(data$Practice_bin)>=2
  platform <- unique_n(data$Platform)>=2
  if(practice) terms <- c(terms,"Practice_bin")
  if(platform) terms <- c(terms,"Platform")
  if(platform) terms <- c(terms,"PlatformConnect_indicator:Time_since_baseline_1y")
  if(platform && practice) terms <- c(terms,"Platform:Practice_bin")
  # The model specification has correlated subject intercepts and slopes.
  # Failure cannot silently change the scientific model to an intercept-only fit.
  if(sum(table(data$Subject_ID)>=2)<20) stop("Fewer than 20 subjects support the specified random slope")
  terms <- c(terms,"(1 + Time_since_baseline_1y | Subject_ID)")
  scale <- c(if(platform) "Platform",if(practice) "Practice_bin")
  phi <- if(length(scale)) paste(scale,collapse=" + ") else "1"
  formula <- brms::bf(as.formula(paste("NM_response ~",paste(terms,collapse=" + "))),
                      as.formula(paste("phi ~",phi)))
  prior <- c(brms::set_prior("normal(0, 0.5)",class="b"),
    brms::set_prior("exponential(1)",class="sd"),
    brms::set_prior("lkj(2)",class="cor"),
    brms::set_prior("student_t(3, 0, 2.5)",class="Intercept"),
    brms::set_prior("student_t(3, 0, 2.5)",class="Intercept",dpar="phi"),
    brms::set_prior("gamma(0.01, 0.01)",class="kappa"))
  if(any(grepl("^s\\(",terms))) prior <- c(prior,brms::set_prior("exponential(1)",class="sds"))
  if(length(scale)) prior <- c(prior,brms::set_prior("normal(0, 0.4)",class="b",dpar="phi"))
  list(formula=formula,prior=prior,family=brms::xbeta(),kind="primary")
}

model_design_audit <- function(model,data) {
  standata <- brms::make_standata(model$formula,data=data,family=model$family,prior=model$prior)
  full_audit <- cantab_audit_standata(standata,model$kind)
  x <- cbind(standata$X,standata$Xs)
  out <- data.frame(n=nrow(x),columns=ncol(x),rank=qr(x)$rank)
  out$full_rank <- out$columns==out$rank
  if(!out$full_rank) stop("Nonidentifiable fixed/spline model design")
  attr(out,"full_audit") <- full_audit
  out
}
