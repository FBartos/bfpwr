source(file.path("simulations", "scripts", "common.R"))

args <- bfpwr_sim_parse_args()
corpus_root <- args[["corpus-root"]]
design_case_id <- args[["design-case"]]
chunk_id <- as.integer(args[["chunk"]])
case_set <- if (is.null(args[["case-set"]])) "production" else args[["case-set"]]
if (is.null(corpus_root) || is.null(design_case_id) || is.na(chunk_id)) {
    stop("usage: Rscript simulations/scripts/run_chunk.R --corpus-root <path> --design-case <id> --chunk <id>")
}

bfpwr_sim_source_package()
design <- bfpwr_sim_find_design_case(design_case_id,
                                     bfpwr_sim_design_case_set(case_set))
chunk <- bfpwr_sim_run_design_chunk(design, chunk_id)
design_dir <- bfpwr_sim_design_dir(corpus_root, design)
bfpwr_sim_write_rds(design, file.path(design_dir, "design.rds"))
bfpwr_sim_write_rds(chunk, file.path(design_dir, "chunks",
                                     sprintf("chunk-%04d.rds", chunk_id)))
cat("wrote chunk", chunk_id, "for", design_case_id, "\n")
