source(file.path("simulations", "scripts", "common.R"))

args <- bfpwr_sim_parse_args()
corpus_root <- args[["corpus-root"]]
design_case_id <- args[["design-case"]]
case_set <- if (is.null(args[["case-set"]])) "production" else args[["case-set"]]
write_trajectories <- isTRUE(args[["write-trajectories"]])
if (is.null(corpus_root) || is.null(design_case_id)) {
    stop("usage: Rscript simulations/scripts/merge_chunks.R --corpus-root <path> --design-case <id>")
}

design <- bfpwr_sim_find_design_case(design_case_id,
                                     bfpwr_sim_design_case_set(case_set))
design_dir <- bfpwr_sim_design_dir(corpus_root, design)
index <- bfpwr_sim_write_design_chunk_index(corpus_root, design)

if (write_trajectories) {
    trajectories <- bfpwr_sim_read_design_trajectories(corpus_root, design)
    bfpwr_sim_validate_design_trajectories(design, trajectories)
    bfpwr_sim_write_rds(trajectories, file.path(design_dir, "trajectories.rds"))
    bfpwr_sim_write_hashes(c(file.path(design_dir, "design.rds"),
                             file.path(design_dir, "chunk-index.rds"),
                             file.path(design_dir, "chunk-index.csv"),
                             file.path(design_dir, "trajectories.rds"),
                             bfpwr_sim_design_chunk_files(corpus_root, design)),
                           file.path(design_dir, "sha256.txt"))
    cat("merged", nrow(index), "chunks and wrote trajectories.rds for",
        design_case_id, "\n")
} else {
    cat("finalized", nrow(index), "chunk files for", design_case_id,
        "without writing trajectories.rds\n")
}
