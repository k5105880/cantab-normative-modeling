# Software checks and limits

## Run checks from R

From the repository root:

```r
source("scripts/run_tests.R")
source("scripts/verify_release.R")
```

The test runner covers CSV creation, malformed-input rejection, participant-ID preservation, age derivation, numerical score calculations, grouped folds, model design, output assembly, execution and resume guards. Common result-export tests verify all six schemas, exact score text, participant IDs, row identity, reference-model routing, precision checks and protection of edited results. Input preparation ignores unrelated export fields, and the tests verify that these fields do not alter score eligibility. The tests use synthetic fixtures and do not fit production models.

The release verifier compares the file inventory and SHA-256 hashes and parses every R source file. Run it on a pristine release. Version control's `.git` directory is ignored. The complete list of additional tests is in `scripts/test_files.txt`. Model-design tests construct synthetic inputs by default; optional input arguments allow the same checks on a separately supplied dataset.

## Validation scope

All six R entry points have passed no-MCMC input, model-design and Stan-syntax checks. Tests also cover input without optional reference files, numeric-looking identifiers, the 14-column input schema, and source/output paths containing spaces. The [validation report](../provenance/validation_report.json) records these checks and their dates.

The available execution host used Linux. Native Windows and macOS were not available for execution tests. CmdStan and a compatible C++ toolchain are required for fitting on every operating system.

Model definitions were compared before and after source reorganization: formulas, prior contents, generated Stan code and Stan data were checked for the full-control model and all five training folds. Deterministic numerical checks and available saved synthetic-fit predictions were compared without refitting production models. Completed score tables were also used to verify exact score text, participant IDs and row order through the shared exporter.

Regression tests cover unused participant-factor levels in automatic fold assignment, cohorts containing only healthy controls, isolation of R session settings, and output directories accessed through symbolic links. These are software edge cases; corrected automatic folds must be generated in a fresh run. The validation report distinguishes unchanged model calculations from these execution fixes.

## Configuration and usability checks

The latest full validation passed all 43 public test scripts and all six no-MCMC entry-point checks, including restoration of the caller environment. This documentation revision leaves the executable analysis code and its defaults unchanged. Documentation links and R example syntax were checked again.

The active settings were traced to their use, and their remaining default values were compared. Removing ineffective draw-selection settings from movement/latency preserved full-reference z values and complete held-out prediction objects in saved synthetic-fit comparisons. These short fits test arithmetic equivalence; they do not establish production convergence. Model definitions, priors and sampling/recovery values remain unchanged. The settings and diagnostics guides are provided in English and Japanese; local links, code fences and R example syntax were checked.

Configuration checks reject unknown, duplicate, missing and incorrectly typed A′ settings before analysis. Reference tests verify unique fold-file selection, unchanged cohort labels and copied-fold hashes. Assembly tests retain ineligible rows whose platform has no eligible controls, while continuing to reject scored rows without that support. A′ prepared numerical inputs and all six full/held-out Stan datasets were compared after the cleanup.

CSV structure tests reject unclosed quotations, mismatched field counts and blank records before creating an analysis input. Valid quoted cells, Unicode, UTF-8 BOM, CRLF line endings and a missing final newline remain supported. The documented export-to-input workflow was executed in the `project/working_directory/` layout. Valid input tables and constructed CSVs were compared with the previous release without changing model calculations.

Diagnosis values are HC, MDD and BD throughout input validation, reference-cohort selection and results. Tests verify correct control/patient score routing and unchanged output labels for all six outcomes. Prepared numerical values, eligibility and row selections were compared with the preceding release after accounting for diagnosis spelling; A′ automatic folds were also identical. The IQ documentation distinguishes an estimated IQ value from the reading test used to obtain it; the FSIQ column and its numerical transformation are unchanged.

## Evaluating your fitted model

Input and software checks do not establish model suitability. For each fitted run, inspect convergence diagnostics, score precision and held-out calibration, and evaluate whether the covariate ranges and supported design match your data. Full production MCMC is not part of the packaged test suite. A new set of participants or folds requires a new run; resume only with unchanged input, code, software and settings.

The repository includes schema templates, synthetic examples and code.
