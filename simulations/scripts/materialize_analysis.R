source(file.path("simulations", "scripts", "common.R"))

args <- bfpwr_sim_parse_args()
corpus_root <- args[["corpus-root"]]
analysis_case_id <- args[["analysis-case"]]
case_set <- if (is.null(args[["case-set"]])) "production" else args[["case-set"]]
include_bf_trajectories <- isTRUE(args[["include-bf-trajectories"]])
if (is.null(corpus_root) || is.null(analysis_case_id)) {
    stop("usage: Rscript simulations/scripts/materialize_analysis.R --corpus-root <path> --analysis-case <id>")
}

bfpwr_sim_source_package()
analysis <- bfpwr_sim_find_analysis_case(
    analysis_case_id, bfpwr_sim_analysis_case_set(case_set))
design <- bfpwr_sim_find_design_case(
    analysis$design_case_id, bfpwr_sim_design_case_set(case_set))
analysis_dir <- bfpwr_sim_analysis_dir(corpus_root, analysis)
bf_chunk_dir <- if (include_bf_trajectories) {
    file.path(analysis_dir, "bf_chunks")
} else {
    NULL
}
materialized <- bfpwr_sim_materialize_analysis_from_chunks(
    analysis, design, corpus_root,
    include_bf_trajectories = include_bf_trajectories,
    bf_chunk_dir = bf_chunk_dir)
bfpwr_sim_write_rds(analysis, file.path(analysis_dir, "analysis.rds"))
bfpwr_sim_write_rds(materialized$outcomes, file.path(analysis_dir, "outcomes.rds"))
bfpwr_sim_write_rds(materialized$summary, file.path(analysis_dir, "summary.rds"))
if (include_bf_trajectories) {
    bf_index <- data.frame(
        analysis_case_id = analysis$analysis_case_id,
        design_case_id = design$design_case_id,
        chunk_id = design$chunks$chunk_id,
        file = file.path("bf_chunks", basename(materialized$bf_chunk_files)),
        bytes = file.info(materialized$bf_chunk_files)$size,
        sha256 = unname(tools::sha256sum(materialized$bf_chunk_files)),
        stringsAsFactors = FALSE
    )
    bfpwr_sim_write_rds(bf_index, file.path(analysis_dir, "bf-chunk-index.rds"))
    utils::write.csv(bf_index, file.path(analysis_dir, "bf-chunk-index.csv"),
                     row.names = FALSE)
}
utils::write.csv(materialized$summary, file.path(analysis_dir, "summary.csv"),
                 row.names = FALSE)
hash_files <- c(file.path(analysis_dir, "analysis.rds"),
                file.path(analysis_dir, "outcomes.rds"),
                file.path(analysis_dir, "summary.rds"),
                file.path(analysis_dir, "summary.csv"))
if (include_bf_trajectories) {
    hash_files <- c(hash_files,
                    file.path(analysis_dir, "bf-chunk-index.rds"),
                    file.path(analysis_dir, "bf-chunk-index.csv"),
                    materialized$bf_chunk_files)
}
bfpwr_sim_write_hashes(hash_files, file.path(analysis_dir, "sha256.txt"))
cat("materialized", analysis_case_id, "\n")
