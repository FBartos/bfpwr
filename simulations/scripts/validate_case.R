source(file.path("simulations", "scripts", "common.R"))

args <- bfpwr_sim_parse_args()
corpus_root <- args[["corpus-root"]]
analysis_case_id <- args[["analysis-case"]]
case_set <- if (is.null(args[["case-set"]])) "production" else args[["case-set"]]
if (is.null(corpus_root) || is.null(analysis_case_id)) {
    stop("usage: Rscript simulations/scripts/validate_case.R --corpus-root <path> --analysis-case <id>")
}

analysis <- bfpwr_sim_find_analysis_case(
    analysis_case_id, bfpwr_sim_analysis_case_set(case_set))
summary <- readRDS(file.path(bfpwr_sim_analysis_dir(corpus_root, analysis),
                             "summary.rds"))
validation <- bfpwr_sim_validate_summary(summary)
print(validation)
ok_cols <- grep("_ok$", names(validation), value = TRUE)
if (any(validation[ok_cols] == FALSE, na.rm = TRUE)) {
    stop("validation failed for ", analysis_case_id)
}
