args <- gsub('~+~',' ',commandArgs(FALSE),fixed=TRUE)
root <- dirname(dirname(normalizePath(sub('^--file=','',args[grepl('^--file=',args)][1]))))
source(file.path(root,'R','precision_helpers.R'))
g <- normal_quadrature(64L)
stopifnot(abs(sum(g$weight)-1)<1e-12, abs(sum(g$weight*g$node))<1e-12,
          abs(sum(g$weight*g$node^2)-1)<1e-12)
# Exact marginal probability in a probit-normal model.
for(s in c(0,.2,1,2)) stopifnot(abs(sum(g$weight*pnorm(.7+s*g$node))-
                                        pnorm(.7/sqrt(1+s^2)))<1e-8)
# Uniform beta: interior CDF and boundary mid-P have known values.
u <- xbeta_bounds_logit(c(0,.4,1),eta=0,phi=2,kappa=0)
stopifnot(all.equal(u$lower,c(0,.4,1)),all.equal(u$upper,c(0,.4,1)))
b <- xbeta_bounds_logit(c(0,.4,1),eta=0,phi=2,kappa=.1)
stopifnot(max(abs(b$lower-c(0,.5/1.2,1.1/1.2)))<1e-12,
          max(abs(b$upper-c(.1/1.2,.5/1.2,1)))<1e-12)
# The subject-specific intercept/slope covariance is retained in marginal scoring.
stopifnot(abs(random_effect_variance(2,3,.4,1)-(4+9+4.8))<1e-12)
cat('PASS: normal quadrature, beta boundaries and random-effect covariance.\n')
