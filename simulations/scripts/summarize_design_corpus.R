source(file.path("simulations", "scripts", "common.R"))

args <- bfpwr_sim_parse_args()
corpus_root <- args[["corpus-root"]]
case_set <- if (is.null(args[["case-set"]])) "production" else args[["case-set"]]
output_file <- args[["output-file"]]
if (is.null(corpus_root)) {
    stop("usage: Rscript simulations/scripts/summarize_design_corpus.R --corpus-root <path>")
}

designs <- bfpwr_sim_design_case_set(case_set)
rows <- lapply(designs, function(design) {
    design_dir <- bfpwr_sim_design_dir(corpus_root, design)
    index_file <- file.path(design_dir, "chunk-index.rds")
    chunk_files <- bfpwr_sim_design_chunk_files(corpus_root, design)
    missing_chunks <- chunk_files[!file.exists(chunk_files)]

    if (file.exists(index_file)) {
        index <- readRDS(index_file)
        n_index_chunks <- nrow(index)
        total_rows <- sum(index$n_rows)
        total_bytes <- sum(index$bytes)
        max_chunk_bytes <- max(index$bytes)
    } else {
        n_index_chunks <- 0L
        total_rows <- NA_real_
        total_bytes <- NA_real_
        max_chunk_bytes <- NA_real_
    }

    data.frame(
        design_case_id = design$design_case_id,
        family = design$test_family,
        tier = design$tier,
        look_grid = design$look_grid_name,
        expected_chunks = nrow(design$chunks),
        indexed_chunks = n_index_chunks,
        missing_chunks = length(missing_chunks),
        total_rows = total_rows,
        total_bytes = total_bytes,
        max_chunk_bytes = max_chunk_bytes,
        has_trajectories = file.exists(file.path(design_dir, "trajectories.rds")),
        stringsAsFactors = FALSE
    )
})
summary <- do.call(rbind, rows)

if (!is.null(output_file)) {
    dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
    utils::write.csv(summary, output_file, row.names = FALSE)
}

cat("designs:", nrow(summary), "\n")
cat("finalized designs:", sum(summary$indexed_chunks == summary$expected_chunks), "\n")
cat("expected chunks:", sum(summary$expected_chunks), "\n")
cat("indexed chunks:", sum(summary$indexed_chunks), "\n")
cat("missing chunk files:", sum(summary$missing_chunks), "\n")
cat("monolithic trajectories:", sum(summary$has_trajectories), "\n")
cat("total indexed bytes:", sum(summary$total_bytes, na.rm = TRUE), "\n")
cat("max chunk bytes:", max(summary$max_chunk_bytes, na.rm = TRUE), "\n")
cat("max chunk MB:", round(max(summary$max_chunk_bytes, na.rm = TRUE) / 1024^2, 2), "\n")
cat("\nby tier/family:\n")
print(stats::aggregate(
    cbind(designs = rep(1, nrow(summary)),
          expected_chunks = summary$expected_chunks,
          indexed_chunks = summary$indexed_chunks,
          total_bytes = summary$total_bytes) ~ tier + family,
    data = summary,
    FUN = sum
), row.names = FALSE)

bad <- summary$indexed_chunks != summary$expected_chunks |
    summary$missing_chunks > 0 |
    summary$has_trajectories
if (any(bad)) {
    cat("\nproblem designs:\n")
    print(summary[bad, c("design_case_id", "expected_chunks", "indexed_chunks",
                         "missing_chunks", "has_trajectories")],
          row.names = FALSE)
    quit(save = "no", status = 1)
}

cat("\ndesign corpus ok\n")
