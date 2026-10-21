# All records in this example are synthetic. This demonstrates input creation only.
.example_sources <- unlist(lapply(sys.frames(), function(f) f$ofile), use.names = FALSE)
.example_self <- .example_sources[basename(.example_sources) == "create_example.R"]
if (!length(.example_self)) stop("Run source('data/example/create_example.R') from R.")
.example_dir <- dirname(normalizePath(tail(.example_self, 1L), winslash = "/"))
source(file.path(.example_dir, "..", "..", "R", "input_data.R"))
example_csv <- file.path(tempfile("cantab_synthetic_"), "CANTAB_connect_RS6_rawscores.csv")
create_cantab_input(file.path(.example_dir, "participants.csv"),
                    file.path(.example_dir, "assessments.csv"), example_csv)
print(read_cantab_input(example_csv)[, cantab_required_columns])
message("Synthetic CSV location: ", example_csv)
message("Input-format example only: this dataset is too small for model fitting.")
rm(.example_sources, .example_self, .example_dir)
