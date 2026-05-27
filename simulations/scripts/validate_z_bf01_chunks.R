source(file.path("simulations", "scripts", "common.R"))

args <- bfpwr_sim_parse_args()
corpus_root <- args[["corpus-root"]]
case_set <- if (is.null(args[["case-set"]])) "production" else args[["case-set"]]
chunk_id <- if (is.null(args[["chunk-id"]])) 1L else as.integer(args[["chunk-id"]])

default_bf_prior_ids <- c(
    "z-bf01-point-null0-pm0p2-psd0-on-z-dpoint-0-usd1-short",
    "z-bf01-normal-null0-pm0-psd0p70710678-on-z-dpoint-0-usd1-short",
    "z-bf01-normal-null0-pm10-psd20-on-z-dpoint-10-usd50-short"
)
bf_prior_ids <- if (is.null(args[["bf-priors"]])) {
    default_bf_prior_ids
} else {
    trimws(strsplit(args[["bf-priors"]], ",", fixed = TRUE)[[1]])
}

if (is.null(corpus_root) || !nzchar(corpus_root) || is.na(chunk_id)) {
    stop("usage: Rscript simulations/scripts/validate_z_bf01_chunks.R ",
         "--corpus-root <path> [--case-set production] [--chunk-id 1] ",
         "[--bf-priors id1,id2]")
}
corpus_root <- bfpwr_sim_require_corpus_root(
    corpus_root,
    asset_set = "chunk-validation",
    purpose = "z BF01 chunk validation")

bfpwr_sim_source_package()

same_log_bf <- function(observed, expected, tolerance = 1e-12) {
    if (length(observed) != length(expected)) return(FALSE)
    both_nan <- is.nan(observed) & is.nan(expected)
    both_pos_inf <- is.infinite(observed) & is.infinite(expected) &
        observed > 0 & expected > 0
    both_neg_inf <- is.infinite(observed) & is.infinite(expected) &
        observed < 0 & expected < 0
    both_finite <- is.finite(observed) & is.finite(expected)

    ok <- both_nan | both_pos_inf | both_neg_inf
    ok[both_finite] <- abs(observed[both_finite] -
                               expected[both_finite]) <= tolerance
    all(ok)
}

designs <- bfpwr_sim_design_case_set(case_set)
bf_priors <- bfpwr_sim_bf_prior_case_set(case_set)
results <- vector("list", length(bf_prior_ids))

for (i in seq_along(bf_prior_ids)) {
    bf_prior_id <- bf_prior_ids[[i]]
    bf_prior <- bfpwr_sim_find_bf_prior_case(bf_prior_id, bf_priors)
    design <- bfpwr_sim_find_design_case(bf_prior$design_case_id, designs)
    if (!identical(bf_prior$test_family, "z") ||
        !identical(bf_prior$bf_type, "normal")) {
        stop("validate_z_bf01_chunks.R only supports z bf01 BF priors: ",
             bf_prior_id)
    }
    if (!chunk_id %in% design$chunks$chunk_id) {
        stop("chunk-id ", chunk_id, " is not part of design ",
             design$design_case_id)
    }

    bf_dir <- bfpwr_sim_bf_prior_dir(corpus_root, bf_prior)
    expected_chunk_files <- bfpwr_sim_bf_chunk_files(corpus_root, bf_prior,
                                                     design)
    missing <- expected_chunk_files[!file.exists(expected_chunk_files)]
    if (length(missing) > 0) {
        bfpwr_sim_stop_missing_corpus_files(
            corpus_root,
            missing = missing,
            asset_set = "chunk-validation",
            label = paste("BF chunk files for", bf_prior_id))
    }

    stored_rows <- bfpwr_sim_read_bf_prior_chunk(corpus_root, bf_prior,
                                                 design, chunk_id)
    bfpwr_sim_validate_bf_rows(bf_prior, design, stored_rows)

    trajectories <- bfpwr_sim_read_design_chunk(corpus_root, design, chunk_id)
    bfpwr_sim_validate_design_trajectories(design, trajectories)
    recomputed_rows <- bfpwr_sim_bf_rows(bf_prior, trajectories)

    key_cols <- c("design_case_id", "replicate_id", "chunk_id", "look", "n")
    if (!identical(stored_rows[key_cols], recomputed_rows[key_cols])) {
        stop("stored and recomputed BF row keys differ for ", bf_prior_id)
    }
    if (!same_log_bf(stored_rows$log_bf01, recomputed_rows$log_bf01)) {
        stop("stored and recomputed log_bf01 values differ for ",
             bf_prior_id)
    }

    results[[i]] <- data.frame(
        bf_prior_id = bf_prior_id,
        design_case_id = design$design_case_id,
        chunk_id = chunk_id,
        n_expected_chunks = length(expected_chunk_files),
        n_rows_checked = nrow(stored_rows),
        chunk_bytes = file.info(file.path(bf_dir, "bf_chunks",
                                          sprintf("chunk-%04d.rds",
                                                  chunk_id)))$size,
        stringsAsFactors = FALSE
    )
}

results <- do.call(rbind, results)
print(results)
cat("validated cached z-bf01 BF chunks against current bf01() output\n")
