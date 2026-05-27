source(file.path("simulations", "scripts", "common.R"))

args <- bfpwr_sim_parse_args()
corpus_root <- args[["corpus-root"]]
plan_index <- args[["plan-index"]]
if (is.null(plan_index) || identical(plan_index, "")) {
    plan_index <- Sys.getenv("PBS_ARRAY_INDEX", unset = NA_character_)
}
plan_index <- as.integer(plan_index)
case_set <- if (is.null(args[["case-set"]])) "production" else args[["case-set"]]
force <- isTRUE(args[["force"]])
if (is.null(corpus_root) || is.na(plan_index)) {
    stop("usage: Rscript simulations/scripts/run_design_plan_chunk.R --corpus-root <path> --plan-index <index>")
}

bfpwr_sim_source_package()

designs <- bfpwr_sim_design_case_set(case_set)
plan <- bfpwr_sim_run_plan(designs, analyses = list())
plan <- plan[plan$phase == "generate_design", , drop = FALSE]
if (plan_index < 1L || plan_index > nrow(plan)) {
    stop("plan index out of range: ", plan_index, " (valid range: 1-", nrow(plan), ")")
}

row <- plan[plan_index, , drop = FALSE]
design <- bfpwr_sim_find_design_case(row$design_case_id, designs)
chunk_id <- row$chunk_id[[1]]
design_dir <- bfpwr_sim_design_dir(corpus_root, design)
design_file <- file.path(design_dir, "design.rds")
chunk_file <- file.path(design_dir, "chunks",
                        sprintf("chunk-%04d.rds", chunk_id))

if (file.exists(chunk_file) && !force) {
    valid <- tryCatch({
        chunk <- readRDS(chunk_file)
        bfpwr_sim_validate_design_trajectories(design, chunk)
        expected <- seq(row$replicate_start[[1]], row$replicate_end[[1]])
        identical(sort(unique(chunk$replicate_id)), expected) &&
            all(chunk$chunk_id == chunk_id)
    }, error = function(e) {
        message("existing chunk failed validation and will be regenerated: ",
                conditionMessage(e))
        FALSE
    })
    if (isTRUE(valid)) {
        bfpwr_sim_write_rds(design, design_file)
        cat("skipped existing valid chunk", chunk_id, "for",
            design$design_case_id, "\n")
        quit(save = "no", status = 0)
    }
}

chunk <- bfpwr_sim_run_design_chunk(design, chunk_id)
bfpwr_sim_write_rds(design, design_file)
bfpwr_sim_write_rds(chunk, chunk_file)
cat("wrote plan index", plan_index, "chunk", chunk_id, "for",
    design$design_case_id, "\n")
