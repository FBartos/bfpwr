source(file.path("simulations", "scripts", "common.R"))

args <- bfpwr_sim_parse_args()
corpus_root <- args[["corpus-root"]]
fixture_set_id <- args[["fixture-set"]]
family <- args[["family"]]
case_set <- if (is.null(args[["case-set"]])) "production" else args[["case-set"]]
output_dir <- args[["output-dir"]]
n_fractions <- if (is.null(args[["n-fractions"]])) {
    c(0.05, 0.25, 0.50, 0.75, 1.00)
} else {
    as.numeric(strsplit(args[["n-fractions"]], ",", fixed = TRUE)[[1]])
}

if (is.null(corpus_root) || is.null(fixture_set_id) || is.null(family) ||
    is.null(output_dir)) {
    stop("usage: Rscript simulations/scripts/materialize_fixed_reference_plot_diagnostics.R ",
         "--corpus-root <path> --fixture-set <fixture-id> --family <family> ",
         "--output-dir <path> [--case-set production] ",
         "[--n-fractions 0.05,0.25,0.50,0.75,1.00]")
}

corpus_root <- bfpwr_sim_require_corpus_root(
    corpus_root,
    asset_set = "fixture-tests",
    purpose = "fixed-reference plot diagnostics")
bfpwr_sim_source_package()

fixture <- bfpwr_sim_read_fixture_summary(
    corpus_root = corpus_root,
    fixture_set_id = fixture_set_id,
    family = family,
    mode = "fixed")
fixture$tail_summary <- bfpwr_sim_select_fixed_reference_plot_rows(
    fixture$tail_summary,
    n_fractions = n_fractions)

cat("materializing fixed reference plot diagnostics for",
    fixture_set_id, "\n")
cat("rows selected:", nrow(fixture$tail_summary), "\n")

t0 <- proc.time()[["elapsed"]]
reference <- bfpwr_sim_validate_fixed_references(
    fixture = fixture,
    designs = bfpwr_sim_design_case_set(case_set),
    bf_priors = bfpwr_sim_bf_prior_case_set(case_set))
elapsed <- proc.time()[["elapsed"]] - t0

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
prefix <- file.path(output_dir, fixture_set_id)
utils::write.csv(reference$summary,
                 paste0(prefix, "-reference-plot-summary.csv"),
                 row.names = FALSE)
utils::write.csv(reference$group_summary,
                 paste0(prefix, "-reference-plot-groups.csv"),
                 row.names = FALSE)
utils::write.csv(reference$outliers,
                 paste0(prefix, "-reference-plot-outliers.csv"),
                 row.names = FALSE)
utils::write.csv(reference$diagnostic_rows,
                 paste0(prefix, "-reference-plot-diagnostic-rows.csv"),
                 row.names = FALSE)

cat("elapsed seconds:", elapsed, "\n")
print(reference$summary)
