source(file.path("simulations", "scripts", "common.R"))

args <- bfpwr_sim_parse_args()
corpus_root <- args[["corpus-root"]]
split_plan_file <- args[["split-plan-file"]]
merge_plan_file <- args[["merge-plan-file"]]
merge_plan_index <- args[["merge-plan-index"]]
if (is.null(merge_plan_index) || identical(merge_plan_index, "")) {
    merge_plan_index <- Sys.getenv("PBS_ARRAY_INDEX", unset = NA_character_)
}
merge_plan_index <- as.integer(merge_plan_index)
case_set <- if (is.null(args[["case-set"]])) "production" else args[["case-set"]]
force <- isTRUE(args[["force"]])

if (is.null(corpus_root) || is.null(split_plan_file) ||
    is.null(merge_plan_file) || is.na(merge_plan_index)) {
    stop("usage: Rscript simulations/scripts/merge_bf_prior_split_plan_chunk.R ",
         "--corpus-root <path> --split-plan-file <csv> ",
         "--merge-plan-file <csv> --merge-plan-index <index>")
}

bfpwr_sim_source_package()

merge_plan <- utils::read.csv(merge_plan_file, stringsAsFactors = FALSE)
split_plan <- utils::read.csv(split_plan_file, stringsAsFactors = FALSE)
if (merge_plan_index < 1L || merge_plan_index > nrow(merge_plan)) {
    stop("merge plan index out of range: ", merge_plan_index,
         " (valid range: 1-", nrow(merge_plan), ")")
}
row <- merge_plan[merge_plan_index, , drop = FALSE]

designs <- bfpwr_sim_design_case_set(case_set)
bf_priors <- bfpwr_sim_bf_prior_case_set(case_set)
design <- bfpwr_sim_find_design_case(row$design_case_id[[1]], designs)
bf_prior <- bfpwr_sim_find_bf_prior_case(row$bf_prior_id[[1]], bf_priors)
chunk_id <- as.integer(row$chunk_id[[1]])
replicate_start <- as.integer(row$replicate_start[[1]])
replicate_end <- as.integer(row$replicate_end[[1]])

final_file <- bfpwr_sim_expected_bf_chunk_file(corpus_root, bf_prior, chunk_id)
if (!force && bfpwr_sim_bf_file_valid(final_file, bf_prior, design, chunk_id,
                                      replicate_start, replicate_end)) {
    cat("skipped existing valid merged BF chunk", chunk_id, "for",
        bf_prior$bf_prior_id, "\n")
    quit(save = "no", status = 0)
}

parts <- split_plan[split_plan$plan_index == row$plan_index[[1]], ,
                    drop = FALSE]
parts <- parts[order(parts$subchunk_id), , drop = FALSE]
if (nrow(parts) != row$n_subchunks[[1]]) {
    stop("merge plan expects ", row$n_subchunks[[1]], " subchunks but found ",
         nrow(parts), " split-plan rows")
}

part_files <- vapply(seq_len(nrow(parts)), function(i) {
    bfpwr_sim_bf_prior_split_chunk_file(
        corpus_root, bf_prior, chunk_id, parts$subchunk_id[[i]])
}, character(1))
missing <- part_files[!file.exists(part_files) | file.info(part_files)$size <= 0]
if (length(missing) > 0L) {
    stop("missing BF split chunk files: ", paste(missing, collapse = ", "))
}

bf_parts <- vector("list", length(part_files))
for (i in seq_along(part_files)) {
    bfpwr_sim_validate_bf_file_range(
        part_files[[i]], bf_prior, design, chunk_id,
        parts$replicate_start[[i]], parts$replicate_end[[i]])
    bf_parts[[i]] <- readRDS(part_files[[i]])
}
bf_rows <- do.call(rbind, bf_parts)
bf_rows <- bf_rows[order(bf_rows$replicate_id, bf_rows$look), , drop = FALSE]
bfpwr_sim_validate_bf_rows(bf_prior, design, bf_rows)

expected <- seq.int(replicate_start, replicate_end)
observed <- sort(unique(bf_rows$replicate_id))
if (!identical(observed, expected)) {
    stop("merged BF chunk does not contain expected replicate range")
}
if (!all(bf_rows$chunk_id == chunk_id)) {
    stop("merged BF chunk contains unexpected chunk_id values")
}

bf_dir <- bfpwr_sim_bf_prior_dir(corpus_root, bf_prior)
bfpwr_sim_write_rds(bf_prior, file.path(bf_dir, "bf-prior.rds"))
tmp_final_file <- paste0(final_file, ".tmp.", Sys.getpid())
on.exit(unlink(tmp_final_file), add = TRUE)
bfpwr_sim_write_rds(bf_rows, tmp_final_file)
if (!file.rename(tmp_final_file, final_file)) {
    stop("failed to move merged BF chunk into place: ", final_file)
}

cat("merged", nrow(parts), "BF split chunks into chunk", chunk_id,
    "for", bf_prior$bf_prior_id, "\n")
