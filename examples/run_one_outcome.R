# Run from project/working_directory/.
# project/working_directory/ から実行。

# Edit these four settings / 変更する4項目:
outcome <- "RTI_reaction"
input_file <- "../input/CANTAB_connect_RS6_rawscores.csv"
output_dir <- "../output/RTI_reaction_check"
mode <- "check"   # check: validate; fresh: fit in a new directory; resume: continue fitting

# For fresh, choose a new output directory such as "../output/RTI_reaction_fit".
# For resume, use the same fit directory, input, code and software.
# freshは新しい保存先へ。resumeは中断した推定と同じ保存先・条件で実行します。
source("run_analysis.R")
run_cantab(outcome = outcome, input_file = input_file, output_dir = output_dir,
           mode = mode)
