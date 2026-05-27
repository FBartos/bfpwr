source("simulations/scripts/common.R")

args <- bfpwr_sim_parse_args()
corpus_root <- args[["corpus-root"]]
fixture_set_id <- args[["fixture-set"]]
case_set <- if (is.null(args[["case-set"]])) "production" else args[["case-set"]]
available_only <- !isTRUE(args[["all-registry-priors"]])
validate <- isTRUE(args[["validate"]])

if (is.null(corpus_root) || is.null(fixture_set_id)) {
    stop("usage: Rscript simulations/scripts/materialize_fixed_fixture.R ",
         "--corpus-root <path> --fixture-set <id>")
}
corpus_root <- bfpwr_sim_require_corpus_root(
    corpus_root,
    asset_set = "chunk-validation",
    purpose = "fixed fixture materialization")

fixture <- bfpwr_sim_find_fixture_spec(fixture_set_id)
if (!identical(fixture$mode, "fixed")) {
    stop("fixture-set is not a fixed fixture: ", fixture_set_id)
}

res <- bfpwr_sim_materialize_fixed_fixture(
    corpus_root = corpus_root,
    fixture = fixture,
    bf_priors = bfpwr_sim_bf_prior_case_set(case_set),
    designs = bfpwr_sim_design_case_set(case_set),
    available_only = available_only,
    validate = validate)

cat("wrote fixed fixture", fixture_set_id, "to", res$fixture_dir, "\n")
print(res$overview)
