# Compare package versions without changing the R session configuration.
local({
expected <- c(brms="2.23.0", cmdstanr="0.9.0", posterior="1.7.0", mgcv="1.9.4",
              dplyr="1.2.1", loo="2.10.1", digest="0.6.39", jsonlite="2.0.0",
              future="1.75.0",
              callr="3.8.0", zip="2.3.3")
actual <- vapply(names(expected), function(p) if(requireNamespace(p,quietly=TRUE))
  as.character(packageVersion(p)) else "MISSING", character(1))
r <- data.frame(package=names(expected), expected=unname(expected), actual=actual,
                matches=actual==expected, row.names=NULL)
print(r, row.names=FALSE)
rv <- as.character(getRversion()); cat("R:",rv,"expected 4.6.1\n")
cv <- if(requireNamespace("cmdstanr",quietly=TRUE))
  tryCatch(as.character(cmdstanr::cmdstan_version()),error=function(e) "UNAVAILABLE") else "UNAVAILABLE"
cat("CmdStan:",cv,"expected 2.39.0\n")
environment_matches_recorded <- all(r$matches) && identical(rv,"4.6.1") && identical(cv,"2.39.0")
if (!environment_matches_recorded)
  message("Versions differ from the tested configuration. Review the installed packages before running.")
invisible(environment_matches_recorded)
})
