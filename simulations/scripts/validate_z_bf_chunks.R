source(file.path("simulations", "scripts", "common.R"))

args <- bfpwr_sim_parse_args()
corpus_root <- args[["corpus-root"]]
case_set <- if (is.null(args[["case-set"]])) "production" else args[["case-set"]]
chunk_id <- if (is.null(args[["chunk-id"]])) 1L else args[["chunk-id"]]
fixture_set <- args[["fixture-set"]]
bf_type <- if (is.null(args[["bf-type"]])) "all" else args[["bf-type"]]
tolerance <- if (is.null(args[["tolerance"]])) 1e-12 else
    as.numeric(args[["tolerance"]])
output_file <- args[["output-file"]]

if (is.null(corpus_root) || !nzchar(corpus_root) ||
    !bf_type %in% c("all", "normal", "moment", "directional") ||
    !is.finite(tolerance) || tolerance < 0) {
    stop("usage: Rscript simulations/scripts/validate_z_bf_chunks.R ",
         "--corpus-root <path> [--case-set production] ",
         "[--fixture-set fixture_id] [--bf-priors id1,id2] ",
         "[--bf-type all|normal|moment|directional] ",
         "[--chunk-id 1|all] [--tolerance 1e-12] ",
         "[--output-file <csv>]")
}

chunk_ids <- if (identical(chunk_id, "all")) {
    "all"
} else {
    as.integer(strsplit(as.character(chunk_id), ",", fixed = TRUE)[[1]])
}
if (!identical(chunk_ids, "all") && any(is.na(chunk_ids))) {
    stop("--chunk-id must be an integer, comma-separated integers, or all")
}
corpus_root <- bfpwr_sim_require_corpus_root(
    corpus_root,
    asset_set = "chunk-validation",
    purpose = "z BF chunk validation")

bfpwr_sim_source_package()

same_log_bf <- function(observed, expected, tolerance = 1e-12) {
    if (length(observed) != length(expected)) return(FALSE)
    both_missing <- is.na(observed) & is.na(expected)
    both_pos_inf <- is.infinite(observed) & is.infinite(expected) &
        observed > 0 & expected > 0
    both_neg_inf <- is.infinite(observed) & is.infinite(expected) &
        observed < 0 & expected < 0
    both_finite <- is.finite(observed) & is.finite(expected)

    ok <- both_missing | both_pos_inf | both_neg_inf
    ok[both_finite] <- abs(observed[both_finite] -
                               expected[both_finite]) <= tolerance
    all(ok)
}

designs <- bfpwr_sim_design_case_set(case_set)
bf_priors <- bfpwr_sim_bf_prior_case_set(case_set)

if (!is.null(args[["bf-priors"]])) {
    requested_ids <- trimws(strsplit(args[["bf-priors"]], ",",
                                     fixed = TRUE)[[1]])
    selected <- lapply(requested_ids, bfpwr_sim_find_bf_prior_case,
                       cases = bf_priors)
} else if (!is.null(fixture_set)) {
    fixture <- bfpwr_sim_find_fixture_spec(fixture_set)
    if (!identical(fixture$family, "z")) {
        stop("--fixture-set must refer to a z fixture: ", fixture_set)
    }
    selected <- bfpwr_sim_select_fixture_bf_priors(
        fixture, bf_priors, corpus_root = corpus_root,
        available_only = TRUE)
} else {
    default_bf_prior_ids <- c(
        "z-bf01-point-null0-pm0p2-psd0-on-z-dpoint-0-usd1-short",
        "z-nmbf01-moment-null0-psd0p35355339-on-z-dpoint-0-usd1-short",
        "z-dirbf01-directional-normal-null0-pm0-psd0p70710678-on-z-dpoint-0-usd1-short"
    )
    selected <- lapply(default_bf_prior_ids, bfpwr_sim_find_bf_prior_case,
                       cases = bf_priors)
}

if (!identical(bf_type, "all")) {
    selected <- selected[vapply(selected, function(x) {
        identical(x$bf_type, bf_type)
    }, logical(1))]
}
if (length(selected) == 0) {
    stop("no z BF priors selected")
}

results <- list()
idx <- 0L
for (bf_prior in selected) {
    if (!identical(bf_prior$test_family, "z") ||
        !bf_prior$bf_type %in% c("normal", "moment", "directional")) {
        stop("selected BF prior is not a supported z BF prior: ",
             bf_prior$bf_prior_id)
    }

    design <- bfpwr_sim_find_design_case(bf_prior$design_case_id, designs)
    ids <- if (identical(chunk_ids, "all")) {
        design$chunks$chunk_id
    } else {
        chunk_ids
    }
    missing_chunks <- setdiff(ids, design$chunks$chunk_id)
    if (length(missing_chunks) > 0) {
        stop("chunk-id not part of design ", design$design_case_id, ": ",
             paste(missing_chunks, collapse = ", "))
    }

    expected_chunk_files <- bfpwr_sim_bf_chunk_files(corpus_root, bf_prior,
                                                     design)
    missing <- expected_chunk_files[!file.exists(expected_chunk_files)]
    if (length(missing) > 0) {
        bfpwr_sim_stop_missing_corpus_files(
            corpus_root,
            missing = missing,
            asset_set = "chunk-validation",
            label = paste("BF chunk files for", bf_prior$bf_prior_id))
    }

    for (id in ids) {
        stored_rows <- bfpwr_sim_read_bf_prior_chunk(
            corpus_root, bf_prior, design, id)
        bfpwr_sim_validate_bf_rows(bf_prior, design, stored_rows)

        trajectories <- bfpwr_sim_read_design_chunk(corpus_root, design, id)
        bfpwr_sim_validate_design_trajectories(design, trajectories)
        recomputed_rows <- bfpwr_sim_bf_rows(bf_prior, trajectories)

        key_cols <- c("design_case_id", "replicate_id", "chunk_id",
                      "look", "n")
        if (!identical(stored_rows[key_cols], recomputed_rows[key_cols])) {
            stop("stored and recomputed BF row keys differ for ",
                 bf_prior$bf_prior_id, " chunk ", id)
        }
        if (!same_log_bf(stored_rows$log_bf01,
                         recomputed_rows$log_bf01,
                         tolerance = tolerance)) {
            stop("stored and recomputed log_bf01 values differ for ",
                 bf_prior$bf_prior_id, " chunk ", id)
        }

        idx <- idx + 1L
        results[[idx]] <- data.frame(
            bf_prior_id = bf_prior$bf_prior_id,
            design_case_id = design$design_case_id,
            bf_type = bf_prior$bf_type,
            chunk_id = id,
            n_expected_chunks = length(expected_chunk_files),
            n_rows_checked = nrow(stored_rows),
            chunk_bytes = file.info(file.path(
                bfpwr_sim_bf_prior_dir(corpus_root, bf_prior), "bf_chunks",
                sprintf("chunk-%04d.rds", id)))$size,
            stringsAsFactors = FALSE
        )
    }
}

results <- do.call(rbind, results)
print(results)
if (!is.null(output_file)) {
    dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
    utils::write.csv(results, output_file, row.names = FALSE)
}
cat("validated cached z BF chunks against current package output\n")
