source(file.path("simulations", "scripts", "common.R"))

args <- bfpwr_sim_parse_args()
corpus_root <- args[["corpus-root"]]
fixture_set_id <- args[["fixture-set"]]
case_set <- if (is.null(args[["case-set"]])) "production" else args[["case-set"]]
force <- isTRUE(args[["force"]])
validate_existing <- isTRUE(args[["validate-existing"]])
skip_index <- isTRUE(args[["skip-index"]])

if (is.null(corpus_root) || is.null(fixture_set_id)) {
    stop("usage: Rscript simulations/scripts/materialize_fixture_bf_priors.R ",
         "--corpus-root <path> --fixture-set <id> [--force] ",
         "[--validate-existing] [--skip-index]")
}

bfpwr_sim_source_package()

fixture <- bfpwr_sim_find_fixture_spec(fixture_set_id)
designs <- bfpwr_sim_design_case_set(case_set)
bf_priors <- bfpwr_sim_bf_prior_case_set(case_set)
selected <- bfpwr_sim_select_fixture_bf_priors(
    fixture = fixture,
    bf_priors = bf_priors,
    corpus_root = corpus_root,
    available_only = FALSE)

if (length(selected) == 0) {
    stop("no BF prior cases selected for fixture: ", fixture_set_id)
}

cat("fixture:", fixture_set_id, "\n")
cat("selected BF priors:", length(selected), "\n")

is_valid_bf_chunk <- function(file, bf_prior, design, chunk) {
    if (!file.exists(file)) return(FALSE)
    if (!validate_existing) return(TRUE)
    tryCatch({
        bf_rows <- readRDS(file)
        bfpwr_sim_validate_bf_rows(bf_prior, design, bf_rows)
        expected <- seq(chunk$replicate_start, chunk$replicate_end)
        identical(sort(unique(bf_rows$replicate_id)), expected) &&
            all(bf_rows$chunk_id == chunk$chunk_id)
    }, error = function(e) {
        message("invalid existing chunk ", file, ": ", conditionMessage(e))
        FALSE
    })
}

selected_ids <- vapply(selected, function(x) x$bf_prior_id, character(1))
design_ids <- unique(vapply(selected, function(x) x$design_case_id,
                            character(1)))
wrote_chunks <- 0L
skipped_chunks <- 0L

for (design_id in design_ids) {
    design <- bfpwr_sim_find_design_case(design_id, designs)
    priors_for_design <- selected[vapply(selected, function(x) {
        identical(x$design_case_id, design_id)
    }, logical(1))]
    cat("design:", design_id, "BF priors:", length(priors_for_design), "\n")

    for (i in seq_len(nrow(design$chunks))) {
        chunk <- design$chunks[i, , drop = FALSE]
        chunk_id <- chunk$chunk_id[[1]]
        trajectories <- NULL

        for (bf_prior in priors_for_design) {
            bf_dir <- bfpwr_sim_bf_prior_dir(corpus_root, bf_prior)
            bf_prior_file <- file.path(bf_dir, "bf-prior.rds")
            bf_chunk_file <- file.path(
                bf_dir, "bf_chunks", sprintf("chunk-%04d.rds", chunk_id))

            if (!force &&
                is_valid_bf_chunk(bf_chunk_file, bf_prior, design, chunk)) {
                bfpwr_sim_write_rds(bf_prior, bf_prior_file)
                skipped_chunks <- skipped_chunks + 1L
                next
            }

            if (is.null(trajectories)) {
                trajectories <- bfpwr_sim_read_design_chunk(
                    corpus_root, design, chunk_id)
                bfpwr_sim_validate_design_trajectories(design, trajectories)
            }

            bf_rows <- bfpwr_sim_bf_rows(bf_prior, trajectories)
            bfpwr_sim_validate_bf_rows(bf_prior, design, bf_rows)
            bfpwr_sim_write_rds(bf_prior, bf_prior_file)
            bfpwr_sim_write_rds(bf_rows, bf_chunk_file)
            wrote_chunks <- wrote_chunks + 1L
        }

        cat("  chunk", chunk_id, "done; wrote", wrote_chunks,
            "skipped", skipped_chunks, "\n")
    }
}

if (!skip_index) {
    for (bf_prior_id in selected_ids) {
        bf_prior <- bfpwr_sim_find_bf_prior_case(bf_prior_id, selected)
        design <- bfpwr_sim_find_design_case(bf_prior$design_case_id, designs)
        bfpwr_sim_write_bf_prior_chunk_index(corpus_root, bf_prior, design)
    }
} else {
    cat("skipped BF-prior chunk index finalization\n")
}

cat("completed fixture BF-prior materialization:", fixture_set_id, "\n")
cat("wrote BF chunks:", wrote_chunks, "\n")
cat("skipped BF chunks:", skipped_chunks, "\n")
