source(file.path("simulations", "scripts", "common.R"))

args <- bfpwr_sim_parse_args()
corpus_root <- args[["corpus-root"]]
case_set <- if (is.null(args[["case-set"]])) "production" else args[["case-set"]]
if (is.null(corpus_root)) {
    stop("usage: Rscript simulations/scripts/validate_corpus.R --corpus-root <path>")
}
corpus_root <- bfpwr_sim_require_corpus_root(
    corpus_root,
    asset_set = "all",
    purpose = "corpus validation")

analyses <- bfpwr_sim_analysis_case_set(case_set)
if (length(analyses) == 0) {
    cat("no analysis cases in case set:", case_set, "\n")
    quit(status = 0)
}
results <- lapply(analyses, function(analysis) {
    summary_file <- file.path(bfpwr_sim_analysis_dir(corpus_root, analysis),
                              "summary.rds")
    if (!file.exists(summary_file)) {
        return(data.frame(analysis_case_id = analysis$analysis_case_id,
                          pH1_ok = FALSE, pH0_ok = FALSE, pInc_ok = FALSE))
    }
    bfpwr_sim_validate_summary(readRDS(summary_file))
})
results <- do.call(rbind, results)
print(results)
ok_cols <- grep("_ok$", names(results), value = TRUE)
if (any(results[ok_cols] == FALSE, na.rm = TRUE)) {
    stop("corpus validation failed")
}
