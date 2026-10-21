# Real numerical fixtures; no sampling and no text/source assertions.
a<-commandArgs(FALSE);root<-dirname(dirname(normalizePath(gsub('~+~',' ',sub('^--file=','',a[grepl('^--file=',a)][1]),fixed=TRUE))))
helper<-file.path(root,'R/state_precision.R')
if(file.exists(helper))source(helper)
if(!exists('swm_state_checked',mode='function'))stop('FAIL: primary state precision is unavailable')
near<-function(a,b,tol=1e-7)stopifnot(all(is.finite(a)),max(abs(a-b))<tol)
err<-function(x)inherits(tryCatch(force(x),error=identity),'error')
metadata<-data.frame(.chain=rep(1:4,each=100),.iteration=rep(1:100,4),.draw=1:400)
mat<-function(x)matrix(x,400,1)
# Fails if random effects are fixed to zero, covariance is omitted, or z is
# averaged before probabilities. Independent adaptive-integrate reference.
for(fam in c('hdg','beta')) {
  y<-if(fam=='hdg')5 else 8;grid<-if(fam=='hdg')NULL else 12
  eta<-if(fam=='hdg')log(4)else qlogis(.4);aux<-if(fam=='hdg')2 else 15
  ref<-function(u) {
    if(fam=='hdg') {
      mu<-exp(eta+.45*u);den<-pnorm(.5,mu,aux,lower.tail=FALSE)
      .12+.88*((pnorm(y-.5,mu,aux)-pnorm(.5,mu,aux))+.5*(pnorm(y+.5,mu,aux)-pnorm(y-.5,mu,aux)))/den
    } else {
      mu<-plogis(eta+.45*u)
      .12+.88*(pbeta((y-.5)/grid,mu*aux,(1-mu)*aux)+pbeta((y+.5)/grid,mu*aux,(1-mu)*aux))/2
    }
  }
  expected<-integrate(function(u)ref(u)*dnorm(u),-9,9,rel.tol=1e-10)$value
  out<-swm_state_checked(y,mat(eta),mat(aux),mat(qlogis(.12)),mat(.45^2),metadata,fam,grid)
  near(out$mid_p,expected);near(out$z,qnorm(expected),1e-6)
  stopifnot(out$numerical_ok,out$posterior_precision_ok,out$posterior_mcse_z==0,
    out$posterior_draws==400,abs(out$z-qnorm(ref(0)))>.005)
}
# Log upper tail stays finite when the ordinary CDF rounds to one.
tail<-swm_state_checked(12,mat(qlogis(.3)),mat(30),mat(-1000),mat(0),metadata,'beta',12)
truth<-qnorm(pbeta(11.5/12,9,21,lower.tail=FALSE,log.p=TRUE)-log(2),lower.tail=FALSE,log.p=TRUE)
near(tail$z,truth);stopifnot(tail$mid_p==1,is.finite(tail$log_mid_p_upper),tail$numerical_ok)
# A lower boundary includes exactly half of the floor/hurdle mass.
floor<-swm_state_checked(0,mat(log(4)),mat(2),mat(qlogis(.12)),mat(.2),metadata,'hdg')
near(floor$mid_p,.06);near(floor$z,qnorm(.06))
# Correlated intercept/slope contribution must retain its cross term.
cov<-array(rep(c(1,.3,.3,.25),each=400),c(400,2,2),dimnames=list(NULL,c('Intercept','Time_since_baseline_1y'),c('Intercept','Time_since_baseline_1y')))
near(swm_state_re_variance(cov,data.frame(Time_since_baseline_1y=c(0,2)))[1,],c(1,3.2))
bad<-cov;bad[,2,2]<--1;stopifnot(err(swm_state_re_variance(bad,data.frame(Time_since_baseline_1y=2))))
# True chain order matters: highly autocorrelated contributions have a larger
# error than pretending the retained draws were IID.
set.seed(28);n<-800
x<-unlist(lapply(1:4,function(i)as.numeric(arima.sim(list(ar=.98),n=n,sd=.02))))
p<-plogis(x);meta<-data.frame(.chain=rep(1:4,each=n),.iteration=rep(1:n,4),.draw=seq_along(x))
z<-qnorm(mean(p))
pr<-swm_state_posterior_precision(matrix(log(p),ncol=1),meta,log(mean(p)),z)
stopifnot(pr$posterior_mcse_z>3*sd(p)/sqrt(length(p))/dnorm(z))
short<-swm_state_posterior_precision(matrix(log(c(.4,.5,.6)),ncol=1),meta[1:3,],log(.5),0)
stopifnot(!short$posterior_precision_ok,is.na(short$posterior_mcse_z))
# Must persist failed row diagnostics before release is refused.
path<-tempfile(fileext='.csv');bad<-tail;bad$row_id_global<-99L;bad$posterior_precision_ok<-FALSE
stopifnot(err(swm_state_write_diagnostics(bad,path)),file.exists(path))
saved<-read.csv(path);stopifnot(saved$row_id_global==99L,!saved$posterior_precision_ok)
# Insufficient quadrature is separately rejected; it is not certified by a
# constant posterior contribution or by a finite z.
sharp<-swm_state_checked(5,mat(log(5)),mat(.02),mat(qlogis(.12)),mat(1),metadata,'hdg',quad_nodes=2,check_nodes=4,max_nodes=4)
stopifnot(!sharp$numerical_ok,err(swm_state_write_diagnostics(sharp,tempfile(fileext='.csv'))))
# An adapter/numerical exception must also leave row IDs and the failure cause.
path<-tempfile(fileext='.csv')
stopifnot(err(swm_state_predict_release(list(family=list(name='unsupported')),data.frame(row_id_global=7),'hdg',path)),file.exists(path))
failed<-read.csv(path);stopifnot(failed$row_id_global==7,!failed$numerical_ok,nzchar(failed$state_error))
cat('PASS: integrated discrete mid-P, direct log tails, covariance, chain-aware MCSE, diagnostic-first precision gates.\n')
