source(file.path("simulations", "scripts", "common.R"))

args <- bfpwr_sim_parse_args()
corpus_root <- args[["corpus-root"]]
case_set <- if (is.null(args[["case-set"]])) "production" else args[["case-set"]]
fixture_set <- if (is.null(args[["fixture-set"]])) "all" else args[["fixture-set"]]
families <- if (is.null(args[["family"]])) NULL else {
    strsplit(args[["family"]], ",", fixed = TRUE)[[1]]
}
modes <- if (is.null(args[["mode"]])) NULL else {
    strsplit(args[["mode"]], ",", fixed = TRUE)[[1]]
}
level <- if (is.null(args[["level"]])) "exhaustive" else args[["level"]]
output_dir <- args[["output-dir"]]
require_all <- isTRUE(args[["require-all"]])
sequential_reference_profile <- if (is.null(args[["sequential-reference-profile"]])) {
    "curated"
} else {
    args[["sequential-reference-profile"]]
}

if (is.null(corpus_root) || !nzchar(corpus_root)) {
    stop("usage: Rscript simulations/scripts/validate_fixture_suite.R ",
         "--corpus-root <path> [--fixture-set all|id1,id2] ",
         "[--family z,t,binomial] [--mode fixed,sequential] ",
         "[--level integrity|exhaustive] [--require-all] ",
         "[--sequential-reference-profile curated|none] ",
         "[--output-dir <path>]")
}
if (!level %in% c("integrity", "exhaustive")) {
    stop("--level must be either integrity or exhaustive")
}
if (!sequential_reference_profile %in% c("curated", "none")) {
    stop("--sequential-reference-profile must be curated or none")
}
corpus_root <- bfpwr_sim_require_corpus_root(
    corpus_root,
    asset_set = "fixture-tests",
    purpose = "fixture-suite validation")

bfpwr_sim_source_package()

specs <- bfpwr_sim_fixture_specs(families = families, modes = modes)
if (!identical(fixture_set, "all")) {
    requested <- trimws(strsplit(fixture_set, ",", fixed = TRUE)[[1]])
    specs <- lapply(requested, bfpwr_sim_find_fixture_spec,
                    specs = bfpwr_sim_fixture_specs())
}
if (length(specs) == 0) {
    stop("no fixture specs selected")
}

fixture_exists <- vapply(specs, function(spec) {
    dir.exists(bfpwr_sim_fixture_dir(corpus_root, spec))
}, logical(1))
missing_specs <- specs[!fixture_exists]
if (length(missing_specs) > 0) {
    missing_ids <- vapply(missing_specs, function(x) x$fixture_set_id,
                          character(1))
    msg <- paste("missing materialized fixture directories:",
                 paste(missing_ids, collapse = ", "))
    if (require_all || !identical(fixture_set, "all")) {
        stop(bfpwr_sim_corpus_instruction_text(
            corpus_root = corpus_root,
            reason = msg,
            asset_set = "fixture-tests"),
            call. = FALSE)
    }
    warning(bfpwr_sim_corpus_instruction_text(
        corpus_root = corpus_root,
        reason = msg,
        asset_set = "fixture-tests"),
        call. = FALSE)
}
specs <- specs[fixture_exists]
if (length(specs) == 0) {
    stop(bfpwr_sim_corpus_instruction_text(
        corpus_root = corpus_root,
        reason = "No materialized fixture directories were selected.",
        asset_set = "fixture-tests"),
        call. = FALSE)
}

results <- vector("list", length(specs))
for (i in seq_along(specs)) {
    spec <- specs[[i]]
    cat("validating", spec$fixture_set_id, "(", spec$family, spec$mode, ")\n")
    results[[i]] <- bfpwr_sim_validate_fixture_summary(
        corpus_root = corpus_root,
        fixture_set_id = spec$fixture_set_id,
        family = spec$family,
        mode = spec$mode,
        deterministic = TRUE,
        references = identical(level, "exhaustive"),
        case_set = case_set,
        sequential_reference_profile = sequential_reference_profile,
        output_dir = output_dir)
    print(results[[i]]$summary)
    if (!is.null(results[[i]]$reference)) {
        print(results[[i]]$reference$summary)
    }
    if (nrow(results[[i]]$failures) > 0) {
        print(results[[i]]$failures)
    }
}

summary <- do.call(rbind, lapply(results, function(x) x$summary))
failures <- do.call(rbind, lapply(results, function(x) {
    if (nrow(x$failures) == 0) {
        return(data.frame())
    }
    cbind(fixture_set_id = x$summary$fixture_set_id[[1]],
          x$failures,
          stringsAsFactors = FALSE)
}))

cat("\nfixture validation summary\n")
print(summary)
if (!is.null(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    utils::write.csv(summary, file.path(output_dir, "fixture-suite-summary.csv"),
                     row.names = FALSE)
    utils::write.csv(failures,
                     file.path(output_dir, "fixture-suite-failures.csv"),
                     row.names = FALSE)
}

if (nrow(failures) > 0 || any(!summary$passed)) {
    stop("fixture validation failed")
}
cat("fixture validation passed\n")
