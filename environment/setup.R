# Source this file, then explicitly call install_cantab_packages().
install_cantab_packages <- function(install_stan = FALSE, stan_version = "2.39.0") {
  if (!is.logical(install_stan) || length(install_stan) != 1L || is.na(install_stan))
    stop("install_stan must be TRUE or FALSE.")
  cran <- c("brms", "posterior", "mgcv", "dplyr", "loo", "digest", "jsonlite",
            "future", "callr", "zip")
  missing <- cran[!vapply(cran, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) install.packages(missing, repos = "https://cloud.r-project.org", Ncpus = 1L)
  if (!requireNamespace("cmdstanr", quietly = TRUE))
    install.packages("cmdstanr",
      repos = c("https://stan-dev.r-universe.dev", "https://cloud.r-project.org"), Ncpus = 1L)
  if (!"xbeta" %in% getNamespaceExports("brms"))
    stop("The installed brms is too old for the A-prime model. Install a version exporting xbeta.")
  if (install_stan) {
    cmdstanr::check_cmdstan_toolchain()
    cmdstanr::install_cmdstan(version = stan_version, cores = 1L)
  }
  invisible(TRUE)
}
