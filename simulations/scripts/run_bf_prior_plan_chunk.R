source(file.path("simulations", "scripts", "common.R"))

args <- bfpwr_sim_parse_args()
corpus_root <- args[["corpus-root"]]
input_corpus_root <- if (is.null(args[["input-corpus-root"]])) {
    corpus_root
} else {
    args[["input-corpus-root"]]
}
plan_index <- args[["plan-index"]]
if (is.null(plan_index) || identical(plan_index, "")) {
    plan_index <- Sys.getenv("PBS_ARRAY_INDEX", unset = NA_character_)
}
plan_index <- as.integer(plan_index)
plan_index_file <- args[["plan-index-file"]]
case_set <- if (is.null(args[["case-set"]])) "production" else args[["case-set"]]
force <- isTRUE(args[["force"]])
allow_slow <- isTRUE(args[["allow-slow"]]) ||
    Sys.getenv("BFPWR_ALLOW_SLOW_BF_PRIOR") %in% c("1", "true", "TRUE", "yes", "YES")
if (is.null(corpus_root) || is.null(input_corpus_root) || is.na(plan_index)) {
    stop("usage: Rscript simulations/scripts/run_bf_prior_plan_chunk.R --corpus-root <path> --plan-index <index>")
}

if (!is.null(plan_index_file)) {
    missing_plan <- utils::read.csv(plan_index_file, stringsAsFactors = FALSE)
    if (!"plan_index" %in% names(missing_plan)) {
        stop("plan index file must contain a plan_index column: ",
             plan_index_file)
    }
    if (plan_index < 1L || plan_index > nrow(missing_plan)) {
        stop("plan index file row out of range: ", plan_index,
             " (valid range: 1-", nrow(missing_plan), ")")
    }
    array_index <- plan_index
    plan_index <- as.integer(missing_plan$plan_index[[array_index]])
    if (is.na(plan_index)) {
        stop("plan_index is NA at row ", array_index, " in ",
             plan_index_file)
    }
    cat("mapped plan index file row", array_index,
        "to BF prior plan index", plan_index, "\n")
}

bfpwr_sim_source_package()

designs <- bfpwr_sim_design_case_set(case_set)
bf_priors <- bfpwr_sim_bf_prior_case_set(case_set)
plan <- bfpwr_sim_bf_prior_run_plan(bf_priors, designs)
if (plan_index < 1L || plan_index > nrow(plan)) {
    stop("plan index out of range: ", plan_index, " (valid range: 1-",
         nrow(plan), ")")
}

row <- plan[plan_index, , drop = FALSE]
design <- bfpwr_sim_find_design_case(row$design_case_id, designs)
bf_prior <- bfpwr_sim_find_bf_prior_case(row$bf_prior_id, bf_priors)
chunk_id <- row$chunk_id[[1]]
bf_dir <- bfpwr_sim_bf_prior_dir(corpus_root, bf_prior)
bf_prior_file <- file.path(bf_dir, "bf-prior.rds")
bf_chunk_file <- file.path(bf_dir, "bf_chunks",
                           sprintf("chunk-%04d.rds", chunk_id))

if (file.exists(bf_chunk_file) && !force) {
    valid <- tryCatch({
        bf_rows <- readRDS(bf_chunk_file)
        bfpwr_sim_validate_bf_rows(bf_prior, design, bf_rows)
        expected <- seq(row$replicate_start[[1]], row$replicate_end[[1]])
        identical(sort(unique(bf_rows$replicate_id)), expected) &&
            all(bf_rows$chunk_id == chunk_id)
    }, error = function(e) {
        message("existing BF chunk failed validation and will be regenerated: ",
                conditionMessage(e))
        FALSE
    })
    if (isTRUE(valid)) {
        bfpwr_sim_write_rds(bf_prior, bf_prior_file)
        cat("skipped existing valid BF chunk", chunk_id, "for",
            bf_prior$bf_prior_id, "\n")
        quit(save = "no", status = 0)
    }
}

slow_reason <- bfpwr_sim_bf_prior_slow_reason(
    bf_prior, design, row$replicate_start[[1]], row$replicate_end[[1]])
if (!allow_slow && !is.na(slow_reason)) {
    stop("BF prior plan index ", plan_index, " is marked slow: ",
         slow_reason, ". Use split_bf_prior_plan.R and ",
         "run_bf_prior_split_plan_chunk.R, or pass --allow-slow to run ",
         "the original full chunk anyway.")
}

trajectories <- bfpwr_sim_read_design_chunk(input_corpus_root, design, chunk_id)
bfpwr_sim_validate_design_trajectories(design, trajectories)
bf_rows <- bfpwr_sim_bf_rows(bf_prior, trajectories)
bfpwr_sim_validate_bf_rows(bf_prior, design, bf_rows)

bfpwr_sim_write_rds(bf_prior, bf_prior_file)
bfpwr_sim_write_rds(bf_rows, bf_chunk_file)
cat("wrote BF prior plan index", plan_index, "chunk", chunk_id, "for",
    bf_prior$bf_prior_id, "\n")
