source(file.path("simulations", "scripts", "common.R"))

args <- bfpwr_sim_parse_args()
corpus_root <- args[["corpus-root"]]
case_set <- if (is.null(args[["case-set"]])) "production" else args[["case-set"]]
output_file <- args[["output-file"]]
validate_existing <- isTRUE(args[["validate-existing"]])
if (is.null(corpus_root)) {
    stop("usage: Rscript simulations/scripts/missing_bf_prior_plan_indices.R --corpus-root <path> [--output-file <csv>] [--validate-existing]")
}

if (validate_existing) {
    bfpwr_sim_source_package()
}

designs <- bfpwr_sim_design_case_set(case_set)
bf_priors <- bfpwr_sim_bf_prior_case_set(case_set)
plan <- bfpwr_sim_bf_prior_run_plan(bf_priors, designs)
plan$plan_index <- seq_len(nrow(plan))
plan$chunk_file <- NA_character_
plan$status <- "missing"
plan$error <- NA_character_

design_ids <- vapply(designs, function(x) x$design_case_id, character(1))
bf_prior_ids <- vapply(bf_priors, function(x) x$bf_prior_id, character(1))

for (i in seq_len(nrow(plan))) {
    design <- designs[[match(plan$design_case_id[[i]], design_ids)]]
    bf_prior <- bf_priors[[match(plan$bf_prior_id[[i]], bf_prior_ids)]]
    chunk_file <- file.path(
        bfpwr_sim_bf_prior_dir(corpus_root, bf_prior),
        "bf_chunks",
        sprintf("chunk-%04d.rds", plan$chunk_id[[i]])
    )
    plan$chunk_file[[i]] <- chunk_file
    if (!file.exists(chunk_file)) {
        next
    }
    if (!validate_existing) {
        plan$status[[i]] <- "present"
        next
    }
    ok <- tryCatch({
        bf_rows <- readRDS(chunk_file)
        bfpwr_sim_validate_bf_rows(bf_prior, design, bf_rows)
        expected <- seq(plan$replicate_start[[i]], plan$replicate_end[[i]])
        observed <- sort(unique(bf_rows$replicate_id))
        if (!identical(observed, expected)) {
            stop("unexpected replicate range")
        }
        if (!all(bf_rows$chunk_id == plan$chunk_id[[i]])) {
            stop("unexpected chunk_id values")
        }
        TRUE
    }, error = function(e) {
        plan$error[[i]] <<- conditionMessage(e)
        FALSE
    })
    plan$status[[i]] <- if (ok) "valid" else "invalid"
}

missing <- plan[plan$status %in% c("missing", "invalid"), ]
missing <- missing[, c(
    "plan_index", "status", "tier", "family", "design_case_id",
    "bf_prior_id", "chunk_id", "replicate_start", "replicate_end",
    "chunk_file", "error"
)]

if (!is.null(output_file)) {
    dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
    utils::write.csv(missing, output_file, row.names = FALSE)
}

cat("BF prior chunks expected:", nrow(plan), "\n")
cat("present:", sum(plan$status == "present"), "\n")
cat("valid:", sum(plan$status == "valid"), "\n")
cat("missing:", sum(plan$status == "missing"), "\n")
cat("invalid:", sum(plan$status == "invalid"), "\n")
if (nrow(missing) > 0) {
    if (nrow(missing) <= 200L) {
        cat("missing/invalid plan indices:",
            paste(missing$plan_index, collapse = ","), "\n")
    } else {
        cat("missing/invalid plan indices: see output CSV; first 200:",
            paste(head(missing$plan_index, 200L), collapse = ","), "\n")
    }
}
