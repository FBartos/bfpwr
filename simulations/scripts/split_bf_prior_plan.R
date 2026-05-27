source(file.path("simulations", "scripts", "common.R"))

args <- bfpwr_sim_parse_args()
case_set <- if (is.null(args[["case-set"]])) "production" else args[["case-set"]]
plan_index_file <- args[["plan-index-file"]]
output_file <- args[["output-file"]]
merge_output_file <- args[["merge-output-file"]]
split_size <- if (is.null(args[["split-size"]])) 25L else as.integer(args[["split-size"]])
only_slow <- isTRUE(args[["only-slow"]])
skip_present <- isTRUE(args[["skip-present"]])
validate_present <- isTRUE(args[["validate-present"]])
corpus_root <- args[["corpus-root"]]

if (is.null(output_file)) {
    stop("usage: Rscript simulations/scripts/split_bf_prior_plan.R ",
         "--output-file <split.csv> [--merge-output-file <merge.csv>] ",
         "[--plan-index-file <csv>] [--split-size <n>] [--only-slow] ",
         "[--skip-present --corpus-root <path>] [--validate-present]")
}
if (!is.finite(split_size) || split_size < 1L) {
    stop("split-size must be a positive integer")
}
if (skip_present && is.null(corpus_root)) {
    stop("--corpus-root is required with --skip-present")
}
if (is.null(merge_output_file)) {
    merge_output_file <- sub("[.]csv$", "-merge.csv", output_file)
    if (identical(merge_output_file, output_file)) {
        merge_output_file <- paste0(output_file, "-merge.csv")
    }
}

plan <- bfpwr_sim_bf_prior_plan_rows(case_set = case_set,
                                     plan_index_file = plan_index_file)
split_plan <- bfpwr_sim_split_bf_prior_plan(
    plan,
    split_size = split_size,
    case_set = case_set,
    only_slow = only_slow,
    corpus_root = corpus_root,
    skip_present = skip_present,
    validate_present = validate_present
)
merge_plan <- bfpwr_sim_bf_prior_split_merge_plan(split_plan)

dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(merge_output_file), recursive = TRUE, showWarnings = FALSE)
utils::write.csv(split_plan, output_file, row.names = FALSE)
utils::write.csv(merge_plan, merge_output_file, row.names = FALSE)

cat("source plan rows:", nrow(plan), "\n")
cat("split plan rows:", nrow(split_plan), "\n")
cat("merge plan rows:", nrow(merge_plan), "\n")
cat("split size:", split_size, "\n")
cat("only slow:", only_slow, "\n")
cat("skip present:", skip_present, "\n")
cat("validate present:", validate_present, "\n")
if (nrow(split_plan) > 0L) {
    cat("\nby family/tier:\n")
    print(stats::aggregate(
        cbind(split_rows = rep(1, nrow(split_plan))) ~ family + tier,
        data = split_plan, FUN = sum
    ), row.names = FALSE)
    cat("\nslow-task reasons:\n")
    print(sort(table(split_plan$slow_reason), decreasing = TRUE))
}
cat("wrote split plan:", output_file, "\n")
cat("wrote merge plan:", merge_output_file, "\n")
