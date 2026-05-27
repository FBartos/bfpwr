source(file.path("simulations", "scripts", "common.R"))

args <- bfpwr_sim_parse_args()
corpus_root <- args[["corpus-root"]]
case_set <- if (is.null(args[["case-set"]])) "production" else args[["case-set"]]
bf_prior_ids <- args[["bf-prior-id"]]
if (is.null(corpus_root) || is.null(bf_prior_ids)) {
    stop("usage: Rscript simulations/scripts/repair_bf_prior_nan_chunks.R ",
         "--corpus-root <path> --bf-prior-id <id1,id2,...> ",
         "[--case-set production]")
}

corpus_root <- bfpwr_sim_require_corpus_root(
    corpus_root,
    asset_set = "chunk-validation",
    purpose = "BF prior chunk repair")
bf_prior_ids <- trimws(strsplit(bf_prior_ids, ",", fixed = TRUE)[[1]])
bf_prior_ids <- bf_prior_ids[nzchar(bf_prior_ids)]

bfpwr_sim_source_package()
designs <- bfpwr_sim_design_case_set(case_set)
bf_priors <- bfpwr_sim_bf_prior_case_set(case_set)

repair_chunk <- function(bf_prior, design, chunk_id) {
    bf_file <- bfpwr_sim_expected_bf_chunk_file(corpus_root, bf_prior,
                                                chunk_id)
    if (!file.exists(bf_file)) {
        bfpwr_sim_stop_missing_corpus_files(
            corpus_root,
            missing = bf_file,
            asset_set = "chunk-validation",
            label = "BF chunk files")
    }
    bf_rows <- readRDS(bf_file)
    bad <- is.na(bf_rows$log_bf01) | is.nan(bf_rows$log_bf01)
    if (!any(bad)) {
        return(0L)
    }

    trajectories <- bfpwr_sim_read_design_chunk(corpus_root, design, chunk_id)
    key <- paste(trajectories$replicate_id, trajectories$look, sep = "\r")
    bad_key <- paste(bf_rows$replicate_id[bad], bf_rows$look[bad],
                     sep = "\r")
    idx <- match(bad_key, key)
    if (any(is.na(idx))) {
        stop("could not align BF rows to design trajectories for ",
             bf_prior$bf_prior_id, " chunk ", chunk_id)
    }

    repaired <- bfpwr_sim_log_bf(bf_prior, trajectories[idx, , drop = FALSE])
    if (any(is.na(repaired) | is.nan(repaired))) {
        stop("recomputed log_bf01 still contains NA/NaN for ",
             bf_prior$bf_prior_id, " chunk ", chunk_id)
    }

    bf_rows$log_bf01[bad] <- repaired
    bfpwr_sim_validate_bf_rows(bf_prior, design, bf_rows)
    bfpwr_sim_write_rds(bf_rows, bf_file)
    sum(bad)
}

summary <- list()
for (id in bf_prior_ids) {
    bf_prior <- bfpwr_sim_find_bf_prior_case(id, bf_priors)
    design <- bfpwr_sim_find_design_case(bf_prior$design_case_id, designs)
    chunk_ids <- design$chunks$chunk_id
    counts <- vapply(chunk_ids, function(chunk_id) {
        repair_chunk(bf_prior, design, chunk_id)
    }, integer(1))
    row <- data.frame(
        bf_prior_id = id,
        chunks_checked = length(chunk_ids),
        rows_repaired = sum(counts),
        stringsAsFactors = FALSE)
    summary[[length(summary) + 1L]] <- row
    cat(id, "repaired", row$rows_repaired, "rows across",
        row$chunks_checked, "chunks\n")
}

out <- do.call(rbind, summary)
print(out, row.names = FALSE)
