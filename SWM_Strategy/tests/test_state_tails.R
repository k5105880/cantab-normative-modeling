args<-commandArgs(FALSE)
root<-dirname(dirname(normalizePath(gsub('~+~',' ',sub('^--file=','',args[grepl('^--file=',args)][1]),fixed=TRUE))))
source(file.path(root,'R/load.R'))
near<-function(a,b,tol=1e-10)stopifnot(all(is.finite(a)),max(abs(a-b))<tol)

# An upper-tail category remains positive after the lower CDF rounds to one.
# The reference uses the beta survival function directly, not a summed CDF.
lp<-swm_logpmf(0:12,.3,30,12,0)
s<-swm_score_pmf(exp(lp),12,12,.5)
expected<-qnorm(pbeta(11.5/12,9,21,lower.tail=FALSE,log.p=TRUE)-log(2),
                lower.tail=FALSE,log.p=TRUE)
near(s$z_mid,expected)
stopifnot(s$mid_p==1,is.finite(s$log_mid_p_upper))

# Log probabilities preserve tails below floating-point probability range.
extreme<-c(0,-1000)
s<-swm_score_pmf(exp(extreme),1,1,.5,log_p=extreme)
near(s$z_mid,qnorm(-1000-log(2),lower.tail=FALSE,log.p=TRUE))
stopifnot(s$numeric_tail_status=='ok')
stopifnot(s$log_score==-1000,s$zero_probability==0L,s$probability_underflow==1L)
boundary<-swm_score_pmf(c(1,0),1,1,.5)
stopifnot(!is.finite(boundary$z_mid),boundary$numeric_tail_status=='numeric_tail_limit',
          boundary$zero_probability==1L,boundary$probability_underflow==0L)

# Posterior probabilities must be averaged before the inverse normal CDF.
logP<-log(rbind(c(.99,.01),c(.2,.8)))
averaged<-swm_average_logpmf(logP)
s<-swm_score_pmf(exp(averaged),1,1,.5,log_p=averaged)
near(s$z_mid,qnorm(.7975))
stopifnot(abs(s$z_mid-mean(qnorm(c(.995,.6))))>.1)
near(s$z_mid,-swm_score_pmf(rev(exp(averaged)),0,1,.5,log_p=rev(averaged))$z_mid)
cat('PASS: direct state log tails, unresolved exact boundaries, probability averaging and tail symmetry.\n')

# Analytic truncated-Normal references independently check the diagnostic
# randomized-PIT moments at ordinary and very narrow probability intervals.
moments<-swm_pit_moments(pnorm(-1),pnorm(1))
near(moments$ez,0)
near(moments$ez2,1-2*dnorm(1)/(pnorm(1)-pnorm(-1)))
floor<-swm_pit_moments(0,.2,-Inf,log(.2),log(.8))
near(floor$ez,-dnorm(qnorm(.2))/.2)
near(floor$ez2,1-qnorm(.2)*dnorm(qnorm(.2))/.2)
tiny<-swm_pit_moments(.4,.4,log(.4),-1000,log(.6))
near(tiny$ez,qnorm(.4));near(tiny$ez2,qnorm(.4)^2)
cat('PASS: stable Strategy randomized-PIT moments preserve ordinary and narrow intervals.\n')

# Randomizing within each category must recover a standard Normal across the
# complete PMF, including asymmetric beta tails and structural floor mass.
for(mu in c(.01,.3,.99))for(zi in c(0,.2,.99))for(grid in c(12,48)) {
  lp<-swm_logpmf(0:grid,mu,30,grid,zi)
  ll<-vapply(0:grid,function(k)swm_logsumexp(lp[seq_len(k)]),numeric(1))
  lu<-vapply(0:grid,function(k)swm_logsumexp(if(k<grid)lp[seq.int(k+2,grid+1)]else numeric()),numeric(1))
  p<-exp(lp);mom<-swm_pit_moments(exp(ll),exp(ll)+p,ll,lp,lu)
  near(sum(p*mom$ez),0,1e-10);near(sum(p*mom$ez2),1,1e-10)
}
# Independent 350-decimal-digit truncated-Normal references. These cases
# expose density cancellation near the transition to narrow-interval handling.
ll<-c(log(.4),log(1e-200),log1p(-1.00001e-23))
lm<-log(c(1e-5,1e-204,1e-28));lu<-c(log(.59999),log1p(-1.0001e-200),log(1e-23))
reference_mean<-c(-.2533341612807424133,-30.205592526120548,9.9730451233882693)
reference_second<-c(.064178197327647624,912.37781985403042,99.461629033138621)
reference_variance<-c(5.5830415409903785e-11,9.1127865106270473e-13,8.2155240723008551e-14)
mom<-swm_pit_moments(exp(ll),exp(ll)+exp(lm),ll,lm,lu)
near(mom$ez,reference_mean,1e-10);near(mom$ez2,reference_second,1e-10)
near(mom$ez2-mom$ez^2,reference_variance,4e-12)
cat('PASS: category-weighted PIT identities and high-precision tail moment references.\n')
