source(file.path("simulations", "scripts", "common.R"))

args <- bfpwr_sim_parse_args()
corpus_root <- args[["corpus-root"]]
plan_index <- args[["plan-index"]]
if (is.null(plan_index) || identical(plan_index, "")) {
    plan_index <- Sys.getenv("PBS_ARRAY_INDEX", unset = NA_character_)
}
plan_index <- as.integer(plan_index)
plan_index_file <- args[["plan-index-file"]]
case_set <- if (is.null(args[["case-set"]])) "production" else args[["case-set"]]
if (is.null(corpus_root) || is.na(plan_index)) {
    stop("usage: Rscript simulations/scripts/finalize_bf_prior_plan.R --corpus-root <path> --plan-index <index>")
}

if (!is.null(plan_index_file)) {
    finalize_plan <- utils::read.csv(plan_index_file, stringsAsFactors = FALSE)
    if (!"plan_index" %in% names(finalize_plan)) {
        stop("plan index file must contain a plan_index column: ",
             plan_index_file)
    }
    if (plan_index < 1L || plan_index > nrow(finalize_plan)) {
        stop("plan index file row out of range: ", plan_index,
             " (valid range: 1-", nrow(finalize_plan), ")")
    }
    array_index <- plan_index
    plan_index <- as.integer(finalize_plan$plan_index[[array_index]])
    if (is.na(plan_index)) {
        stop("plan_index is NA at row ", array_index, " in ",
             plan_index_file)
    }
    cat("mapped finalize plan index file row", array_index,
        "to BF prior index", plan_index, "\n")
}

designs <- bfpwr_sim_design_case_set(case_set)
bf_priors <- bfpwr_sim_bf_prior_case_set(case_set)
if (plan_index < 1L || plan_index > length(bf_priors)) {
    stop("plan index out of range: ", plan_index, " (valid range: 1-",
         length(bf_priors), ")")
}

bf_prior <- bf_priors[[plan_index]]
design <- bfpwr_sim_find_design_case(bf_prior$design_case_id, designs)
index <- bfpwr_sim_write_bf_prior_chunk_index(corpus_root, bf_prior, design)
cat("finalized", bf_prior$bf_prior_id, "with", nrow(index),
    "BF chunk files\n")
