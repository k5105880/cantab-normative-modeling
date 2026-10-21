# Implementation and file integrity

The repository provides six outcome-specific normative modeling workflows. Each folder contains model definitions, sampling configuration, scoring routines and a sourceable R entry point. Shared R helpers prepare and validate inputs; `run_cantab()` dispatches to the selected outcome.

## Computation

Each workflow fits one specified reference model to eligible healthy-control observations. Participant-grouped five-fold cross-validation excludes the scored healthy participant from model training. Other participants are scored with the model fitted to all eligible healthy controls.

Source files are named by their role. Model definitions, probability calculations, fitting, prediction, calibration and execution checks are separated where appropriate. Generic numerical integration helpers support the distributions required by these models. The outcome READMEs identify the relevant files and explain model-specific settings.

R child sessions isolate model execution. Public entry points validate the output path before writing and require a directory outside the code repository. Resume checks compare inputs, source files, software versions and settings. Use a new output directory after changing code or the calculation environment.

Shared exports copy accepted score text without rounding and verify participant identities, row order, precision status and reference-model routing. `normative_scores.csv` and `normative_z` provide the same user-facing interface for all six outcomes.

## Integrity records

`provenance/source_manifest.csv` records the current source paths, roles and SHA-256 values. `SHA256SUMS` and `RELEASE_MANIFEST.json` describe the distributed files. These identify the exact implementation supplied with the release.

Run `source("scripts/verify_release.R")` from the repository root to inspect a pristine release. See [VERIFY.md](VERIFY.md) for numerical checks and their limits. Licensing and software attribution are provided in [LICENSE](../LICENSE) and [CITATION.cff](../CITATION.cff).
