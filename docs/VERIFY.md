# Software checks and limits

## Run the checks

From the repository root:

```r
source("scripts/run_tests.R")
source("scripts/verify_release.R")
```

`run_tests.R` runs the 43 scripts listed in [test_files.txt](../scripts/test_files.txt). They use synthetic fixtures to check input preparation and validation, participant identifiers, grouped cross-validation, model design, probability and score calculations, output assembly, execution isolation, and resume guards. They do not fit production models.

`verify_release.R` checks the complete file inventory, file sizes, SHA-256 hashes and R source syntax against the supplied manifests. Run it on an unmodified distribution; `.git/` is ignored. Missing files or local edits cause a failure. An integrity failure concerns the distributed files and does not by itself establish an error in the statistical calculations.

## Verified scope

The last full software validation passed all 43 tests and all six entry points in `mode = "check"`, including input, model-design and Stan-syntax checks without MCMC. These checks were performed for version `2026.10.21-r8`. Subsequent documentation and packaging revisions leave the executable analysis code and model settings unchanged. The current distribution is checked separately for file integrity, documentation links and R example syntax.

Model validation also included comparisons of formulas, priors, generated Stan code and Stan data for the full-control fit and all five training folds. Deterministic numerical calculations, available saved synthetic-fit predictions and completed score exports were compared. The [validation record](../provenance/validation_report.json) distinguishes these earlier checks from checks performed for the current distribution.

## Limits and checks for your data

Full production MCMC was not rerun for these distribution revisions, and the packaged tests do not establish posterior convergence or reproduce a complete fitted analysis. Saved synthetic fits test calculation consistency, not production convergence. Execution was checked on Linux; native Windows and macOS were not available for testing.

For each fitted run, inspect convergence, score precision and held-out calibration using the [diagnostic guide](DIAGNOSTICS.md). Assess whether the covariate ranges and model design match your data. Input validation and passing software tests do not establish model suitability. Resume only with unchanged input, code, software and settings.
