# CANTAB normative modeling

R code for preparing CANTAB data, fitting reference models in healthy controls, and expressing each assessment as a **normative z score**. Scores express a value's position in the predicted distribution; see the [interpretation guide](docs/USAGE.md). Positive scores indicate higher raw performance values; the interpretation depends on the outcome. 

## Start here

1. [Create the input CSV](docs/DATA.md): prepare participant and assessment tables, map columns, and check units.
2. [Set up R and Stan](environment/README.md): install the required packages and compiler tools.
3. [Run an analysis and read the results](docs/USAGE.md): arguments, output columns, terminology and common errors.
4. [Settings and source-file roles](docs/SETTINGS.md): what to edit, what the internal constants mean, and how settings are applied.
5. [Inspect the diagnostics](docs/DIAGNOSTICS.md): exact output paths, status fields and checks after fitting.

These workflows fit models to **your supplied healthy-control data**. The repository does not include a pre-trained reference model. Supported diagnosis codes are `HC` (healthy control), `MDD` (major depressive disorder) and `BD` (bipolar disorder); platforms are `RS6` and `Connect`.

| Outcome argument / folder | Measure | Direction of `normative_z` |
|---|---|---|
| `RTI_reaction` | Reaction time | Higher is worse |
| `RTI_movement` | Movement time | Higher is worse |
| `RVP_A_prime` | Attention sensitivity (A prime) | Higher is better |
| `RVP_latency` | Response latency | Higher is worse |
| `SWM_Between_errors` | Spatial working-memory errors | Higher is worse |
| `SWM_Strategy` | Spatial working-memory strategy | Higher is worse |

## Run one outcome

Run from `project/working_directory/` using the [documented folder layout](docs/USAGE.md) and a prepared input CSV.

```r
source("run_analysis.R")
input <- "../input/CANTAB_connect_RS6_rawscores.csv"

# Check input and model design without fitting.
run_cantab("RTI_reaction", input, "../output/reaction_check", mode = "check")

# Fit and score; use a different, new directory.
run_cantab("RTI_reaction", input, "../output/reaction_fit", mode = "fresh")

# Read the result while preserving participant identifiers.
scores <- read.csv("../output/reaction_fit/normative_scores.csv",
                   colClasses = c(Subject_ID = "character"))
```

All six outcomes use the same main result name, **`normative_scores.csv`**, inside their own output directories, and the same score column, **`normative_z`**. `score_dictionary.csv` explains every result column; `analysis_summary.txt` identifies the model and score direction. The other generated files support model diagnostics and computation.

Change the outcome argument to another folder name in the table. The editable [R example](examples/run_one_outcome.R) collects the user settings in one place. For interrupted runs, see the guide's `resume` instructions.

## Examples, checks and license

```r
source("data/example/create_example.R") # practice creating the CSV using synthetic records
source("scripts/run_tests.R")          # software tests; no production model fitting
source("scripts/verify_release.R")     # integrity check for an unchanged release
```

The small synthetic example demonstrates the input format; it is not suitable for fitting these models. See [validation scope](docs/VERIFY.md) and [implementation notes](docs/PROVENANCE.md). The code uses the [MIT License](LICENSE); software citation information is in [CITATION.cff](CITATION.cff).

