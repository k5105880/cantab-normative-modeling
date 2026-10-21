// Incomplete-beta continued fraction (DLMF 8.17.22), evaluated by modified
// Lentz iteration. Called on the convergence-favoured side of the distribution.
// This evaluates the same beta CDF; it is not a Gaussian approximation.
real swm_beta_logcdf_cf_left(real x, real a, real b) {
  real qab = a + b;
  real qap = a + 1;
  real qam = a - 1;
  real c = 1;
  real d = 1 - qab * x / qap;
  real h;
  real tiny = 1e-300;
  if (abs(d) < tiny) d = d < 0 ? -tiny : tiny;
  d = 1 / d;
  h = d;
  for (m in 1:10000) {
    real m2 = 2.0 * m;
    real aa = m * (b - m) * x / ((qam + m2) * (a + m2));
    real delta;
    d = 1 + aa * d;
    c = 1 + aa / c;
    if (abs(d) < tiny) d = d < 0 ? -tiny : tiny;
    if (abs(c) < tiny) c = c < 0 ? -tiny : tiny;
    d = 1 / d;
    h *= d * c;
    aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2));
    d = 1 + aa * d;
    c = 1 + aa / c;
    if (abs(d) < tiny) d = d < 0 ? -tiny : tiny;
    if (abs(c) < tiny) c = c < 0 ? -tiny : tiny;
    d = 1 / d;
    delta = d * c;
    h *= delta;
    // At integer b a coefficient can vanish before its parameter derivative
    // has converged. Do not terminate on an early primal-only equality.
    if (m >= 32 && abs(delta - 1) < 2e-14) {
      if (h <= 0 || is_inf(h) || is_nan(h))
        reject("Invalid incomplete-beta continued fraction");
      return a * log(x) + b * log1m(x) - lbeta(a, b) - log(a) + log(h);
    }
  }
  reject("Incomplete-beta continued fraction did not converge");
}

real swm_beta_logcdf(real x, real a, real b) {
  real log_lower_bound;
  if (a <= 0 || b <= 0 || is_inf(a) || is_inf(b) || is_nan(a) || is_nan(b))
    reject("Beta shape parameters must be positive and finite");
  if (x <= 0) return negative_infinity();
  if (x >= 1) return 0;
  // Integral bound: for b>=1, (1-t)^(b-1)>=(1-x)^(b-1) on [0,x];
  // for b<1 the same factor is >=1. Thus this is a LOWER bound on log CDF.
  log_lower_bound = a * log(x) - log(a) - lbeta(a, b);
  if (b >= 1) log_lower_bound += (b - 1) * log1m(x);
  // beta_cdf uses scaled derivatives, unlike the beta(a,b)-based beta_lcdf
  // derivatives in CmdStan 2.39. Never call it when its CDF may underflow:
  // branching after creating an invalid autodiff node would be too late.
  if (log_lower_bound > -600 && a + b <= 10000)
    return log(beta_cdf(x | a, b));
  if (x <= (a + 1) / (a + b + 2))
    return swm_beta_logcdf_cf_left(x, a, b);
  return log1m_exp(swm_beta_logcdf_cf_left(1 - x, b, a));
}

// Each integer is a category index on a 12- or 48-interval score grid.
real swm_disc_beta_lpmf(int k, real mu, real phi, int grid) {
  real a = mu * phi;
  real b = (1 - mu) * phi;
  real lo = fmax(0, (k - 0.5) / grid);
  real hi = fmin(1, (k + 0.5) / grid);
  real log_hi;
  if (k < 0 || k > grid) return negative_infinity();
  if (k == 0) return swm_beta_logcdf(hi, a, b);
  // Beta symmetry avoids cancellation in log(1 - inc_beta(...)).
  if (k == grid) return swm_beta_logcdf(1 - lo, b, a);
  log_hi = swm_beta_logcdf(hi, a, b);
  if (log_hi > log(0.5))
    return log_diff_exp(swm_beta_logcdf(1 - lo, b, a), swm_beta_logcdf(1 - hi, b, a));
  return log_diff_exp(log_hi, swm_beta_logcdf(lo, a, b));
}
real swm_floor_beta_lpmf(int k, real mu, real phi, real zi, int grid) {
  real lp = swm_disc_beta_lpmf(k | mu, phi, grid);
  if (k == 0) return log_sum_exp(log(zi), log1m(zi) + lp);
  return log1m(zi) + lp;
}
