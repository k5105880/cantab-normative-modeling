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
