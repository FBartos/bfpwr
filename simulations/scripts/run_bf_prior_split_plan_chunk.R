source(file.path("simulations", "scripts", "common.R"))

args <- bfpwr_sim_parse_args()
corpus_root <- args[["corpus-root"]]
input_corpus_root <- if (is.null(args[["input-corpus-root"]])) {
    corpus_root
} else {
    args[["input-corpus-root"]]
}
split_plan_file <- args[["split-plan-file"]]
split_plan_index <- args[["split-plan-index"]]
if (is.null(split_plan_index) || identical(split_plan_index, "")) {
    split_plan_index <- Sys.getenv("PBS_ARRAY_INDEX", unset = NA_character_)
}
split_plan_index <- as.integer(split_plan_index)
case_set <- if (is.null(args[["case-set"]])) "production" else args[["case-set"]]
force <- isTRUE(args[["force"]])

if (is.null(corpus_root) || is.null(input_corpus_root) ||
    is.null(split_plan_file) || is.na(split_plan_index)) {
    stop("usage: Rscript simulations/scripts/run_bf_prior_split_plan_chunk.R ",
         "--corpus-root <path> --split-plan-file <csv> ",
         "--split-plan-index <index>")
}

bfpwr_sim_source_package()

split_plan <- utils::read.csv(split_plan_file, stringsAsFactors = FALSE)
if (split_plan_index < 1L || split_plan_index > nrow(split_plan)) {
    stop("split plan index out of range: ", split_plan_index,
         " (valid range: 1-", nrow(split_plan), ")")
}
row <- split_plan[split_plan_index, , drop = FALSE]

designs <- bfpwr_sim_design_case_set(case_set)
bf_priors <- bfpwr_sim_bf_prior_case_set(case_set)
design <- bfpwr_sim_find_design_case(row$design_case_id[[1]], designs)
bf_prior <- bfpwr_sim_find_bf_prior_case(row$bf_prior_id[[1]], bf_priors)
chunk_id <- as.integer(row$chunk_id[[1]])
replicate_start <- as.integer(row$replicate_start[[1]])
replicate_end <- as.integer(row$replicate_end[[1]])
subchunk_id <- as.integer(row$subchunk_id[[1]])

valid_bf_file <- function(paths, replicate_start, replicate_end) {
    any(vapply(unique(paths), bfpwr_sim_bf_file_valid, logical(1),
               bf_prior = bf_prior, design = design, chunk_id = chunk_id,
               replicate_start = replicate_start,
               replicate_end = replicate_end))
}

final_files <- c(
    bfpwr_sim_expected_bf_chunk_file(input_corpus_root, bf_prior, chunk_id),
    bfpwr_sim_expected_bf_chunk_file(corpus_root, bf_prior, chunk_id)
)
if (!force && valid_bf_file(
    final_files, row$chunk_replicate_start[[1]],
    row$chunk_replicate_end[[1]])) {
    cat("skipped split plan index", split_plan_index,
        "because final BF chunk already exists for", bf_prior$bf_prior_id,
        "chunk", chunk_id, "\n")
    quit(save = "no", status = 0)
}

partial_file <- bfpwr_sim_bf_prior_split_chunk_file(
    corpus_root, bf_prior, chunk_id, subchunk_id)
partial_files <- c(
    bfpwr_sim_bf_prior_split_chunk_file(
        input_corpus_root, bf_prior, chunk_id, subchunk_id),
    partial_file
)
if (!force && valid_bf_file(partial_files, replicate_start, replicate_end)) {
    cat("skipped existing valid BF split chunk", subchunk_id, "for",
        bf_prior$bf_prior_id, "chunk", chunk_id, "\n")
    quit(save = "no", status = 0)
}

trajectories <- bfpwr_sim_read_design_chunk(input_corpus_root, design, chunk_id)
keep <- trajectories$replicate_id >= replicate_start &
    trajectories$replicate_id <= replicate_end
trajectories <- trajectories[keep, , drop = FALSE]
bfpwr_sim_validate_design_trajectories(design, trajectories)
observed <- sort(unique(trajectories$replicate_id))
expected <- seq.int(replicate_start, replicate_end)
if (!identical(observed, expected)) {
    stop("filtered design trajectories do not contain expected replicate range")
}

bf_rows <- bfpwr_sim_bf_rows(bf_prior, trajectories)
bfpwr_sim_validate_bf_rows(bf_prior, design, bf_rows)
bfpwr_sim_write_rds(bf_rows, partial_file)

cat("wrote BF split plan index", split_plan_index,
    "plan index", row$plan_index[[1]],
    "chunk", chunk_id,
    "subchunk", subchunk_id,
    "replicates", replicate_start, "-", replicate_end,
    "for", bf_prior$bf_prior_id, "\n")
