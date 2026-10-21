# Publish this repository

For the documented project layout, the repository root is `project/working_directory/`. The prepared ZIP contains a `working_directory/` folder; place it inside your own `project/` folder for local use. Upload the contents of `working_directory/` so `README.md`, `LICENSE`, `run_analysis.R` and the six outcome folders appear at the top level. The adjacent `input/` and `output/` folders are not part of the upload. Run `source("scripts/verify_release.R")` from R before uploading a pristine release.

When creating the GitHub repository, do not generate another README, license or .gitignore: these are already supplied. Include the hidden `.gitignore` and `.gitattributes` files. The latter preserves file bytes across Git checkouts. For browser uploads, split the folder contents into batches of at most 100 files. Upload source files and folders so readers can browse them; a ZIP may additionally be attached to a tagged release.

After all files are present, create a release such as `v2026.10.21-r11`. Add the actual repository URL to `CITATION.cff` when available. A persistent software-release identifier can also be added if one is created. If files are edited after this prepared release, regenerate the manifest and hashes as part of a new version.

The repository contains normative-model R code, schema templates, synthetic input examples, model settings and tests. Code is MIT licensed. The implementation and integrity records are described in [PROVENANCE.md](PROVENANCE.md).

Official instructions: [create a repository](https://docs.github.com/en/repositories/creating-and-managing-repositories/creating-a-new-repository), [upload files](https://docs.github.com/en/repositories/working-with-files/managing-files/adding-a-file-to-a-repository), [create a release](https://docs.github.com/en/repositories/releasing-projects-on-github/managing-releases-in-a-repository).
