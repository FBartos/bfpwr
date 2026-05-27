source(file.path("simulations", "scripts", "common.R"))

args <- bfpwr_sim_parse_args()
corpus_root <- args[["corpus-root"]]
case_set <- if (is.null(args[["case-set"]])) "production" else args[["case-set"]]
output_file <- args[["output-file"]]
if (is.null(corpus_root)) {
    stop("usage: Rscript simulations/scripts/summarize_bf_prior_corpus.R --corpus-root <path>")
}

designs <- bfpwr_sim_design_case_set(case_set)
bf_priors <- bfpwr_sim_bf_prior_case_set(case_set)
rows <- lapply(bf_priors, function(bf_prior) {
    design <- bfpwr_sim_find_design_case(bf_prior$design_case_id, designs)
    bf_dir <- bfpwr_sim_bf_prior_dir(corpus_root, bf_prior)
    index_file <- file.path(bf_dir, "bf-chunk-index.rds")
    chunk_files <- bfpwr_sim_bf_chunk_files(corpus_root, bf_prior, design)
    present_chunks <- chunk_files[file.exists(chunk_files)]
    missing_chunks <- chunk_files[!file.exists(chunk_files)]

    if (file.exists(index_file)) {
        index <- readRDS(index_file)
        n_index_chunks <- nrow(index)
        total_rows <- sum(index$n_rows)
        total_bytes <- sum(index$bytes)
        max_chunk_bytes <- max(index$bytes)
        n_finite_log_bf01 <- if ("n_finite_log_bf01" %in% names(index)) {
            sum(index$n_finite_log_bf01)
        } else {
            NA_real_
        }
        n_infinite_log_bf01 <- if ("n_infinite_log_bf01" %in% names(index)) {
            sum(index$n_infinite_log_bf01)
        } else {
            NA_real_
        }
        n_nan_log_bf01 <- if ("n_nan_log_bf01" %in% names(index)) {
            sum(index$n_nan_log_bf01)
        } else {
            NA_real_
        }
    } else {
        n_index_chunks <- 0L
        total_rows <- NA_real_
        total_bytes <- NA_real_
        max_chunk_bytes <- NA_real_
        n_finite_log_bf01 <- NA_real_
        n_infinite_log_bf01 <- NA_real_
        n_nan_log_bf01 <- NA_real_
    }

    data.frame(
        bf_prior_id = bf_prior$bf_prior_id,
        design_case_id = bf_prior$design_case_id,
        family = bf_prior$test_family,
        tier = bf_prior$tier,
        bf_type = bf_prior$bf_type,
        prior_form = bf_prior$prior_form,
        expected_chunks = nrow(design$chunks),
        present_chunks = length(present_chunks),
        indexed_chunks = n_index_chunks,
        missing_chunks = length(missing_chunks),
        total_rows = total_rows,
        n_finite_log_bf01 = n_finite_log_bf01,
        n_infinite_log_bf01 = n_infinite_log_bf01,
        n_nan_log_bf01 = n_nan_log_bf01,
        total_bytes = total_bytes,
        max_chunk_bytes = max_chunk_bytes,
        stringsAsFactors = FALSE
    )
})
summary <- do.call(rbind, rows)

if (!is.null(output_file)) {
    dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
    utils::write.csv(summary, output_file, row.names = FALSE)
}

cat("BF priors:", nrow(summary), "\n")
cat("finalized BF priors:",
    sum(summary$indexed_chunks == summary$expected_chunks), "\n")
cat("expected BF chunks:", sum(summary$expected_chunks), "\n")
cat("present BF chunk files:", sum(summary$present_chunks), "\n")
cat("indexed BF chunks:", sum(summary$indexed_chunks), "\n")
cat("missing BF chunk files:", sum(summary$missing_chunks), "\n")
cat("total indexed bytes:", sum(summary$total_bytes, na.rm = TRUE), "\n")
cat("total finite log BF values:",
    sum(summary$n_finite_log_bf01, na.rm = TRUE), "\n")
cat("total infinite log BF values:",
    sum(summary$n_infinite_log_bf01, na.rm = TRUE), "\n")
cat("total NaN log BF values:",
    sum(summary$n_nan_log_bf01, na.rm = TRUE), "\n")
max_chunk_bytes <- suppressWarnings(max(summary$max_chunk_bytes, na.rm = TRUE))
if (!is.finite(max_chunk_bytes)) max_chunk_bytes <- NA_real_
cat("max BF chunk bytes:", max_chunk_bytes, "\n")
cat("max BF chunk MB:", round(max_chunk_bytes / 1024^2, 2), "\n")
cat("\nby tier/family:\n")
print(stats::aggregate(
    cbind(bf_priors = rep(1, nrow(summary)),
          expected_chunks = summary$expected_chunks,
          present_chunks = summary$present_chunks,
          indexed_chunks = summary$indexed_chunks,
          total_bytes = ifelse(is.na(summary$total_bytes), 0,
                               summary$total_bytes)) ~ tier + family,
    data = summary,
    FUN = sum
), row.names = FALSE)
cat("\nby family/prior form:\n")
print(stats::aggregate(
    cbind(bf_priors = rep(1, nrow(summary)),
          expected_chunks = summary$expected_chunks,
          present_chunks = summary$present_chunks,
          indexed_chunks = summary$indexed_chunks) ~ family + prior_form,
    data = summary,
    FUN = sum
), row.names = FALSE)

bad <- summary$indexed_chunks != summary$expected_chunks |
    summary$missing_chunks > 0
if (any(bad)) {
    cat("\nproblem BF priors:\n")
    bad_rows <- summary[bad, c("bf_prior_id", "expected_chunks",
                               "indexed_chunks", "present_chunks",
                               "missing_chunks")]
    print(head(bad_rows, 50L), row.names = FALSE)
    if (nrow(bad_rows) > 50L) {
        cat("... omitted", nrow(bad_rows) - 50L,
            "additional problem BF priors; see output CSV for details\n")
    }
    quit(save = "no", status = 1)
}

cat("\nBF prior corpus ok\n")
