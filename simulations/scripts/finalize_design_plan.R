source(file.path("simulations", "scripts", "common.R"))

args <- bfpwr_sim_parse_args()
corpus_root <- args[["corpus-root"]]
plan_index <- args[["plan-index"]]
if (is.null(plan_index) || identical(plan_index, "")) {
    plan_index <- Sys.getenv("PBS_ARRAY_INDEX", unset = NA_character_)
}
plan_index <- as.integer(plan_index)
case_set <- if (is.null(args[["case-set"]])) "production" else args[["case-set"]]
write_trajectories <- isTRUE(args[["write-trajectories"]])
if (is.null(corpus_root) || is.na(plan_index)) {
    stop("usage: Rscript simulations/scripts/finalize_design_plan.R --corpus-root <path> --plan-index <index>")
}

designs <- bfpwr_sim_design_case_set(case_set)
if (plan_index < 1L || plan_index > length(designs)) {
    stop("plan index out of range: ", plan_index, " (valid range: 1-",
         length(designs), ")")
}

design <- designs[[plan_index]]
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
    cat("finalized", design$design_case_id, "with", nrow(index),
        "chunks and wrote trajectories.rds\n")
} else {
    cat("finalized", design$design_case_id, "with", nrow(index),
        "chunk files and no trajectories.rds\n")
}
