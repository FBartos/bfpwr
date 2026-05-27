bfpwr_sim_this_file <- function() {
    if (!is.null(sys.frames()[[1]]$ofile)) {
        return(normalizePath(sys.frames()[[1]]$ofile, winslash = "/", mustWork = TRUE))
    }
    normalizePath("simulations/scripts/common.R", winslash = "/", mustWork = TRUE)
}

bfpwr_sim_repo_root <- function() {
    common <- bfpwr_sim_this_file()
    normalizePath(file.path(dirname(common), "..", ".."), winslash = "/", mustWork = TRUE)
}

bfpwr_sim_source_all <- function(repo_root = bfpwr_sim_repo_root()) {
    files <- c(
        "simulations/R/corpus-release.R",
        "simulations/R/grids.R",
        "simulations/R/cases.R",
        "simulations/registry/factories.R",
        "simulations/R/registry.R",
        "simulations/R/validate.R",
        "simulations/R/generators.R",
        "simulations/R/bayes-factors.R",
        "simulations/R/bf-prior-splitting.R",
        "simulations/R/summarize.R",
        "simulations/R/query.R",
        "simulations/R/io.R",
        "simulations/R/fixture-summaries.R",
        "simulations/R/fixture-validation.R",
        "simulations/registry/designs/z.R",
        "simulations/registry/designs/t.R",
        "simulations/registry/designs/binomial.R",
        "simulations/registry/analyses/z.R",
        "simulations/registry/analyses/t.R",
        "simulations/registry/analyses/binomial.R",
        "simulations/registry/bf-priors/z.R",
        "simulations/registry/bf-priors/t.R",
        "simulations/registry/bf-priors/binomial.R",
        "simulations/registry/design-cases.R",
        "simulations/registry/analysis-cases.R",
        "simulations/registry/bf-prior-cases.R",
        "simulations/registry/fixtures/thresholds.R",
        "simulations/registry/fixtures/z.R",
        "simulations/registry/fixtures/t.R",
        "simulations/registry/fixtures/binomial.R",
        "simulations/registry/fixture-cases.R",
        "simulations/R/smoke.R"
    )
    for (file in file.path(repo_root, files)) {
        source(file, local = .GlobalEnv)
    }
    invisible(files)
}

bfpwr_sim_source_all()
