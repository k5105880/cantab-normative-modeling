# Truth checks are independent one-dimensional integrate() or closed forms.
here <- dirname(normalizePath(gsub('~+~',' ',sub('^--file=', '', grep('^--file=', commandArgs(FALSE), value=TRUE)[1]),fixed=TRUE)))
helper <- file.path(here, 'state_cdf.R')
stopifnot('analytic state CDF helper must exist' = file.exists(helper))
source(helper)
check <- function(ok, label) { if (!isTRUE(ok)) stop(label); cat('PASS', label, '\n') }
reject <- function(expr, pattern) {
  e <- tryCatch({force(expr); NULL}, error=identity)
  check(inherits(e,'error') && grepl(pattern,conditionMessage(e)), paste('reject',pattern))
}
meta <- data.frame(.chain=rep(1:4,each=100),.iteration=rep(1:100,4),.draw=1:400)
mu <- matrix(rep(c(-1,0,2),each=400),400)
student <- state_cdf_checked(c(-2,0,9),mu,matrix(.7,400,3),rep(4,400),matrix(0,400,3),'student_log1p',meta)
truth <- pt((c(-2,0,9)-c(-1,0,2))/.7,4)
check(max(abs(student$p_lower-truth))<1e-12,'Student exact no random effects')
check(max(abs(student$z-qnorm(truth)))<1e-10,'probability averaged before qnorm')
check(all(student$posterior_mcse_z==0),'constant posterior has zero Monte Carlo error')
gamma <- state_cdf_checked(c(1e-8,420,1500),matrix(log(400),400,3),matrix(25,400,3),NULL,matrix(0,400,3),'gamma_raw',meta)
truth <- pgamma(c(1e-8,420,1500),shape=25,scale=16,log.p=TRUE)
check(max(abs(gamma$log_p_lower-truth))<1e-10 && all(is.finite(gamma$z)), 'Gamma far lower tail avoids empirical zero')
check(!any(gamma$clipped),'no probability clipping')
# Nonconstant draws verify probability mixture, not average conditional z.
eta <- matrix(rep(c(-2,1),each=200),400,1)
r <- state_cdf_checked(.3,eta,matrix(.6,400,1),rep(5,400),matrix(0,400,1),'student_log1p',meta)
check(abs(r$z-qnorm(mean(pt((.3-eta)/.6,5))))<1e-10,'posterior mixture truth')
check(abs(r$z-mean(qnorm(pt((.3-eta)/.6,5))))>.1,'mixture differs from mean z')
# Correlated intercept/slope covariance including negative time.
cv <- array(0,c(400,2,2),dimnames=list(NULL,c('Intercept','Time_since_baseline_1y'),c('Intercept','Time_since_baseline_1y')))
cv[,1,1]<-.25;cv[,2,2]<-.09;cv[,1,2]<-cv[,2,1]<--.12
v <- state_re_variance(cv,data.frame(Time_since_baseline_1y=c(-1,0,1,3)))
check(max(abs(v[1,]-(.25+.09*c(-1,0,1,3)^2- .24*c(-1,0,1,3))))<1e-14,'correlated intercept and slope projection')
ri <- cv[,1,1,drop=FALSE]
check(all(state_re_variance(ri,data.frame(Time_since_baseline_1y=c(0,2)))==.25),'random intercept fallback')
for (family in c('gamma_raw','student_log1p')) {
 y <- if(family=='gamma_raw')c(.01,150,400,900,8000) else c(-10,4,5.9,7,20)
 eta <- if(family=='gamma_raw') log(400) else 5.9
 aux <- if(family=='gamma_raw') 25 else .35
 vv <- c(.01,.09,.25,.05,.2)
 out<-state_cdf_checked(y,matrix(eta,400,5),matrix(aux,400,5),if(family=='gamma_raw')NULL else rep(6,400),matrix(rep(vv,each=400),400),family,meta)
 for(j in seq_along(y)) {
  integrand<-if(family=='gamma_raw')function(b) pgamma(exp(log(y[j])+log(aux)-eta-b),shape=aux)*dnorm(b,0,sqrt(vv[j])) else function(b) pt((y[j]-eta-b)/aux,6)*dnorm(b,0,sqrt(vv[j]))
  # Scale integrand to prevent integrate() absolute tolerance hiding tiny tails.
  logtarget<-min(out$log_p_lower[j],out$log_p_upper[j])
  if(out$log_p_upper[j]<out$log_p_lower[j]) integrand<-if(family=='gamma_raw')function(b) pgamma(exp(log(y[j])+log(aux)-eta-b),shape=aux,lower.tail=FALSE)*dnorm(b,0,sqrt(vv[j])) else function(b) pt((y[j]-eta-b)/aux,6,lower.tail=FALSE)*dnorm(b,0,sqrt(vv[j]))
  actual<-integrate(function(b) integrand(b)/exp(logtarget),-Inf,Inf,rel.tol=1e-9,subdivisions=2000)$value
  check(abs(actual-1)<1e-4,paste(family,'independent integral row',j))
 }
 check(all(out$numerical_ok),paste(family,'adaptive integration converged'))
}
# Gaussian limiting Student marginal closed form; broad RE forces refinement.
out <- state_cdf_checked(1,matrix(0,400,1),matrix(1,400,1),rep(1e7,400),matrix(1,400,1),'student_log1p',meta)
check(abs(out$z-1/sqrt(2))<1e-6,'Gaussian convolution limit')
reject(state_cdf_checked(.01,matrix(0,400,1),matrix(.001,400,1),rep(5,400),matrix(25,400,1),'student_log1p',meta,quad_nodes=2,check_nodes=4,max_nodes=4),'converge')
reject(state_cdf_checked(0,matrix(NA_real_,400,1),matrix(1,400,1),rep(5,400),matrix(0,400,1),'student_log1p',meta),'finite')
reject(state_cdf_checked(1,matrix(1,400,1),matrix(0,400,1),NULL,matrix(0,400,1),'gamma_raw',meta),'positive')
reject(state_cdf_checked(1,matrix(1,400,1),matrix(1,400,1),NULL,matrix(-1,400,1),'gamma_raw',meta),'variance')
cv[,1,2]<-cv[,2,1]<-1
reject(state_re_variance(cv,data.frame(Time_since_baseline_1y=0)),'semidefinite')
# MCSE preserves autocorrelation and chain structure, including scrambled metadata.
set.seed(20)
x<-unlist(lapply(1:4,function(k)as.numeric(arima.sim(list(ar=.97),n=1000))))
metalong<-data.frame(.chain=rep(1:4,each=1000),.iteration=rep(1:1000,4),.draw=1:4000)
a<-state_posterior_precision(matrix(log(plogis(x)),4000,1),metalong,log(.5),0)
check(a$posterior_mcse_p > sd(plogis(x))/sqrt(4000)*2,'chain-aware MCSE retains autocorrelation')
ord<-sample(1:4000)
b<-state_posterior_precision(matrix(log(plogis(x[ord])),4000,1),metalong[ord,],log(.5),0)
check(abs(a$posterior_mcse_p-b$posterior_mcse_p)<1e-14,'MCSE uses chain and iteration metadata')
reject(state_posterior_precision(matrix(-1,400,1),meta[-1,],-1,0),'metadata')
# A finite converged quadrature estimate is insufficient when MCMC precision fails.
prefix<-file.path(tempdir(),'state_cdf_test')
state_write_diagnostics(student,prefix,'student_log1p','good')
check(file.exists(paste0(prefix,'_student_log1p_good_state_cdf_diagnostics.csv')),'diagnostic export')
bad<-student;bad$posterior_precision_ok[1]<-FALSE
err<-tryCatch(state_scoring_guard(state_write_diagnostics(bad,prefix,'student_log1p','bad')),error=identity)
check(inherits(err,'state_cdf_error')&&grepl('Scores withheld',conditionMessage(err)),'posterior precision failure is fatal and typed')
check(file.exists(paste0(prefix,'_student_log1p_bad_state_cdf_diagnostics.csv')),'failed precision diagnostics retained')
cat('STATE_CDF_NUMERICAL_SUITE_PASS\n')
