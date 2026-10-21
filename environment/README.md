# R and Stan setup

Model fitting requires R, CmdStan and a working C++ toolchain.

## Install packages

From the repository root:

```r
source("environment/setup.R")
install_cantab_packages()  # installs missing R dependencies; leaves existing versions in place
```

The script uses CRAN and the official Stan R-universe repository for cmdstanr. It does not automatically replace your existing package versions. A′ requires brms with the `xbeta` family. [Official CmdStanR setup](https://mc-stan.org/cmdstanr/articles/cmdstanr.html)

If brms lacks `xbeta`, update it in a fresh R session:

```r
install.packages("brms", repos = "https://cloud.r-project.org")
```

Restart R after updating, then rerun setup and the environment check. This installs the current CRAN version, not the recorded reference version.

## Install or select CmdStan

If CmdStan is already installed in a custom location, select the directory containing its `bin` folder and make it available to child R sessions:

```r
stan_dir <- normalizePath("path/to/cmdstan", winslash = "/", mustWork = TRUE)
Sys.setenv(CMDSTAN = stan_dir)
cmdstanr::set_cmdstan_path(stan_dir)
```

For a new installation, prepare the compiler toolchain and run:

```r
source("environment/setup.R")
install_cantab_packages(install_stan = TRUE)
```

Windows uses the Rtools version corresponding to the installed R version; macOS uses its C++ command-line tools. Follow the toolchain check's instructions for your system. [Official toolchain and installation reference](https://mc-stan.org/cmdstanr/reference/install_cmdstan.html)

## Tested environment

The workflows were checked with R 4.6.1, brms 2.23.0, mgcv 1.9-4, cmdstanr 0.9.0, posterior 1.7.0 and CmdStan 2.39.0. [recorded_R_packages.csv](recorded_R_packages.csv) lists package versions for reference. The list is not an automatically restored dependency lockfile. Execution also uses `callr` and `zip`.

```r
source("environment/check_environment.R")
```

The check reports differences from the tested configuration; matching versions alone does not establish model suitability. Use consistent package versions when resuming a run. Checks were executed on Linux; native Windows and macOS were not available for testing.

Prefer output paths without spaces where required by the CmdStan toolchain.

[日本語の環境準備](README_JA.md) · [Execution guide](../docs/USAGE.md) · [Model settings](../docs/SETTINGS.md) · [Logs and diagnostics](../docs/DIAGNOSTICS.md)
