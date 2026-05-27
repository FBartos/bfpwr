## Shared helper for gated fixture-corpus tests. These helpers load the
## simulation registry and compare cached summaries to package references; they
## are not themselves manuscript calculations.

bfpwr_fixture_guess_repo_root <- function(start = getwd()) {
    cur <- normalizePath(start, winslash = "/", mustWork = FALSE)
    for (i in 0:8) {
        if (file.exists(file.path(cur, "simulations", "scripts",
                                  "common.R"))) {
            return(cur)
        }
        parent <- dirname(cur)
        if (identical(parent, cur)) {
            break
        }
        cur <- parent
    }
    normalizePath(".", winslash = "/", mustWork = FALSE)
}

bfpwr_fixture_corpus_download_instructions <- function(
        corpus_root = Sys.getenv("BFPWR_SIM_CORPUS", unset = ""),
        repo_root = Sys.getenv("BFPWR_SIM_REPO", unset = ""),
        asset_set = "fixture-tests",
        reason = NULL) {
    if (!nzchar(repo_root)) {
        repo_root <- bfpwr_fixture_guess_repo_root()
    }
    if (!nzchar(corpus_root)) {
        corpus_root <- file.path(repo_root, "simulations", "corpus", "v1")
    }
    corpus_root <- normalizePath(corpus_root, winslash = "/", mustWork = FALSE)
    repo_root <- normalizePath(repo_root, winslash = "/", mustWork = FALSE)
    lines <- c()
    if (!is.null(reason) && nzchar(reason)) {
        lines <- c(lines, reason, "")
    }
    paste(c(
        lines,
        paste0("The bfpwr simulation corpus is required at: ", corpus_root),
        "",
        "Download and extract the needed release assets with:",
        paste0("  Rscript simulations/scripts/download_corpus_release.R ",
               "--corpus-root ", shQuote(corpus_root),
               " --asset-set ", asset_set),
        "",
        "Release: https://github.com/FBartos/bfpwr/releases/tag/sim-corpus-v1",
        "",
        "Then run fixture tests with:",
        paste0("  BFPWR_SIM_CORPUS=", corpus_root),
        paste0("  BFPWR_SIM_REPO=", repo_root)
    ), collapse = "\n")
}

bfpwr_fixture_source_simulations <- function(repo_root) {
    simulation_files <- c(
        "simulations/R/corpus-release.R",
        "simulations/R/grids.R",
        "simulations/R/cases.R",
        "simulations/registry/factories.R",
        "simulations/R/registry.R",
        "simulations/R/validate.R",
        "simulations/R/generators.R",
        "simulations/R/bayes-factors.R",
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
    for (file in file.path(repo_root, simulation_files)) {
        if (!file.exists(file)) {
            stop("required simulation source file does not exist: ", file)
        }
        source(file, local = .GlobalEnv)
    }
    invisible(simulation_files)
}

bfpwr_fixture_context <- function(required_fixture_dir = NULL) {
    corpus_root <- Sys.getenv("BFPWR_SIM_CORPUS", unset = "")
    if (!nzchar(corpus_root) || !dir.exists(corpus_root)) {
        return(list(skip = bfpwr_fixture_corpus_download_instructions(
            corpus_root = corpus_root,
            reason = paste0("Set BFPWR_SIM_CORPUS to run large ",
                            "fixture-corpus tests."))))
    }
    corpus_root <- normalizePath(corpus_root, winslash = "/", mustWork = TRUE)

    repo_root <- Sys.getenv("BFPWR_SIM_REPO", unset = "")
    if (!nzchar(repo_root)) {
        repo_root <- normalizePath(file.path(corpus_root, "..", "..", ".."),
                                   winslash = "/", mustWork = FALSE)
    } else {
        repo_root <- normalizePath(repo_root, winslash = "/", mustWork = TRUE)
    }

    common_file <- file.path(repo_root, "simulations", "scripts", "common.R")
    if (!file.exists(common_file)) {
        return(list(skip = bfpwr_fixture_corpus_download_instructions(
            corpus_root = corpus_root,
            repo_root = repo_root,
            reason = paste0("Set BFPWR_SIM_REPO to a checkout containing ",
                            "simulations/scripts/common.R."))))
    }

    if (!is.null(required_fixture_dir)) {
        fixture_dir <- file.path(corpus_root, required_fixture_dir)
        if (!dir.exists(fixture_dir)) {
            return(list(skip = bfpwr_fixture_corpus_download_instructions(
                corpus_root = corpus_root,
                repo_root = repo_root,
                reason = paste0("The fixture directory is missing: ",
                                fixture_dir))))
        }
    }

    old_wd <- setwd(repo_root)
    bfpwr_fixture_source_simulations(repo_root)

    list(
        corpus_root = corpus_root,
        repo_root = repo_root,
        old_wd = old_wd
    )
}

bfpwr_fixture_manifest_status <- function(fixture_dir) {
    manifest_file <- file.path(fixture_dir, "manifest.csv")
    if (!file.exists(manifest_file)) {
        stop("fixture manifest does not exist: ", manifest_file)
    }
    manifest <- utils::read.csv(manifest_file, stringsAsFactors = FALSE)
    required <- c("file", "bytes", "sha256")
    missing <- setdiff(required, names(manifest))
    if (length(missing) > 0) {
        stop("fixture manifest is missing columns: ",
             paste(missing, collapse = ", "))
    }

    paths <- file.path(fixture_dir, manifest$file)
    exists <- file.exists(paths)
    actual_bytes <- rep(NA_real_, length(paths))
    actual_sha256 <- rep(NA_character_, length(paths))
    if (any(exists)) {
        actual_bytes[exists] <- file.info(paths[exists])$size
        actual_sha256[exists] <- unname(tools::sha256sum(paths[exists]))
    }

    data.frame(
        file = manifest$file,
        exists = exists,
        expected_bytes = manifest$bytes,
        actual_bytes = actual_bytes,
        bytes_ok = exists & actual_bytes == manifest$bytes,
        expected_sha256 = manifest$sha256,
        actual_sha256 = actual_sha256,
        sha256_ok = exists & actual_sha256 == manifest$sha256,
        stringsAsFactors = FALSE
    )
}

bfpwr_fixture_expect_fixed_summary_integrity <- function(fixture,
                                                         label = "fixture") {
    expect_true(nrow(fixture$tail_summary) > 0,
                info = paste(label, "tail summary is available"))
    expect_true(nrow(fixture$decision_summary) > 0,
                info = paste(label, "decision summary is available"))
    expect_true(nrow(fixture$search_summary) > 0,
                info = paste(label, "search summary is available"))

    manifest_status <- bfpwr_fixture_manifest_status(fixture$fixture_dir)
    expect_true(all(manifest_status$exists),
                info = paste(label, "manifest files exist"))
    expect_true(all(manifest_status$bytes_ok),
                info = paste(label, "manifest byte sizes match"))
    expect_true(all(manifest_status$sha256_ok),
                info = paste(label, "manifest SHA256 hashes match"))

    failures <- bfpwr_sim_validate_fixed_fixture_deterministic(fixture)
    expect_equal(nrow(failures), 0L,
                 info = paste(label, "fixed-summary deterministic invariants pass"))
    invisible(failures)
}

bfpwr_fixture_expect_sequential_summary_integrity <- function(fixture,
                                                              label = "fixture",
                                                              require_search = FALSE) {
    expect_true(nrow(fixture$cumulative_summary) > 0,
                info = paste(label, "cumulative summary is available"))
    expect_true(nrow(fixture$final_summary) > 0,
                info = paste(label, "final summary is available"))
    if (isTRUE(require_search)) {
        expect_true(nrow(fixture$search_summary) > 0,
                    info = paste(label, "search summary is available"))
    } else {
        expect_equal(nrow(fixture$search_summary), 0L,
                     info = paste(label, "search summary is intentionally empty"))
    }

    manifest_status <- bfpwr_fixture_manifest_status(fixture$fixture_dir)
    expect_true(all(manifest_status$exists),
                info = paste(label, "manifest files exist"))
    expect_true(all(manifest_status$bytes_ok),
                info = paste(label, "manifest byte sizes match"))
    expect_true(all(manifest_status$sha256_ok),
                info = paste(label, "manifest SHA256 hashes match"))

    failures <- bfpwr_sim_validate_sequential_fixture_deterministic(fixture)
    expect_equal(nrow(failures), 0L,
                 info = paste(label,
                              "sequential-summary deterministic invariants pass"))
    invisible(failures)
}

bfpwr_fixture_expect_fixed_search_first_crossing <- function(fixture,
                                                             bf_prior_id,
                                                             tail,
                                                             evidence_threshold,
                                                             target_prob,
                                                             label = "fixture") {
    search_row <- fixture$search_summary[
        fixture$search_summary$bf_prior_id == bf_prior_id &
            fixture$search_summary$tail == tail &
            fixture$search_summary$evidence_threshold ==
                evidence_threshold &
            fixture$search_summary$target_prob == target_prob,
        , drop = FALSE]
    expect_equal(nrow(search_row), 1L,
                 info = paste(label, "selected fixed-summary search row exists"))
    if (nrow(search_row) != 1L) {
        return(invisible(NULL))
    }

    curve <- fixture$tail_summary[
        fixture$tail_summary$bf_prior_id == search_row$bf_prior_id &
            fixture$tail_summary$tail == search_row$tail &
            fixture$tail_summary$evidence_threshold ==
                search_row$evidence_threshold,
        , drop = FALSE]
    curve <- curve[order(curve$n), , drop = FALSE]
    hit <- which(curve$prob >= search_row$target_prob)
    if (length(hit) == 0) {
        expect_false(isTRUE(search_row$achieved),
                     info = paste(label, "search row records no crossing"))
        return(invisible(search_row))
    }

    first_hit <- hit[[1]]
    expect_equal(search_row$n_found, curve$n[[first_hit]],
                 info = paste(label,
                              "fixed-summary search n matches first threshold crossing"))
    expect_equal(search_row$prob_found, curve$prob[[first_hit]],
                 info = paste(label,
                              "fixed-summary search probability matches first threshold crossing"))
    invisible(search_row)
}

bfpwr_fixture_sequential_schedule <- function(fixture, schedule_id) {
    hits <- vapply(fixture$spec$schedules, function(x) {
        identical(x$schedule_id, schedule_id)
    }, logical(1))
    if (sum(hits) != 1L) {
        stop("schedule not found or not unique: ", schedule_id)
    }
    fixture$spec$schedules[[which(hits)]]
}

bfpwr_fixture_z_sequential_reference <- function(fixture,
                                                 row,
                                                 bf_prior,
                                                 design,
                                                 strict = TRUE) {
    schedule <- bfpwr_fixture_sequential_schedule(fixture, row$schedule_id[[1]])
    n <- schedule$n
    se <- design$generation$usd / sqrt(n)
    dp <- bfpwr_sim_design_prior_mean_sd(design)
    prior <- bf_prior$analysis_prior
    dpm <- dp$dpm - prior$null

    if (identical(bf_prior$bf_type, "normal")) {
        pm <- prior$pm - prior$null
        psd <- prior$psd
        type <- "normal"
    } else if (identical(bf_prior$bf_type, "moment")) {
        pm <- NULL
        psd <- prior$psd
        type <- "moment"
    } else if (identical(bf_prior$bf_type, "directional")) {
        pm <- prior$pm - prior$null
        psd <- prior$psd
        type <- "directional"
    } else {
        stop("unsupported z sequential BF type: ", bf_prior$bf_type)
    }

    pbf01seq(k1 = row$k1[[1]],
             k0 = row$k0[[1]],
             se = se,
             n = n,
             pm = pm,
             psd = psd,
             dpm = dpm,
             dpsd = dp$dpsd,
             type = type,
             strict = strict)
}

bfpwr_fixture_z_sequential_reference_case <- function(fixture,
                                                      bf_prior_id,
                                                      schedule_id,
                                                      evidence_threshold,
                                                      designs,
                                                      bf_priors,
                                                      strict = TRUE) {
    final_row <- fixture$final_summary[
        fixture$final_summary$bf_prior_id == bf_prior_id &
            fixture$final_summary$schedule_id == schedule_id &
            fixture$final_summary$evidence_threshold == evidence_threshold,
        , drop = FALSE]
    if (nrow(final_row) != 1L) {
        stop("selected sequential final row is not unique: ",
             paste(bf_prior_id, schedule_id, evidence_threshold,
                   sep = ", "))
    }

    curve <- fixture$cumulative_summary[
        fixture$cumulative_summary$bf_prior_id == bf_prior_id &
            fixture$cumulative_summary$schedule_id == schedule_id &
            fixture$cumulative_summary$evidence_threshold ==
                evidence_threshold,
        , drop = FALSE]
    curve <- curve[order(curve$look), , drop = FALSE]
    if (nrow(curve) != final_row$n_looks[[1]]) {
        stop("selected sequential cumulative curve has wrong length: ",
             paste(bf_prior_id, schedule_id, evidence_threshold,
                   sep = ", "))
    }

    bf_prior <- bfpwr_sim_find_bf_prior_case(bf_prior_id, bf_priors)
    design <- bfpwr_sim_find_design_case(final_row$design_case_id[[1]],
                                         designs)
    reference <- bfpwr_fixture_z_sequential_reference(
        fixture = fixture,
        row = final_row,
        bf_prior = bf_prior,
        design = design,
        strict = strict)

    list(final_row = final_row, curve = curve, reference = reference)
}

bfpwr_fixture_expect_z_sequential_reference <- function(fixture,
                                                        bf_prior_id,
                                                        schedule_id,
                                                        evidence_threshold,
                                                        designs,
                                                        bf_priors,
                                                        strict = TRUE,
                                                        label = "fixture") {
    final_row <- fixture$final_summary[
        fixture$final_summary$bf_prior_id == bf_prior_id &
            fixture$final_summary$schedule_id == schedule_id &
            fixture$final_summary$evidence_threshold == evidence_threshold,
        , drop = FALSE]
    expect_equal(nrow(final_row), 1L,
                 info = paste(label, "selected sequential final row exists",
                              bf_prior_id, schedule_id, evidence_threshold))
    if (nrow(final_row) != 1L) {
        return(invisible(NULL))
    }

    curve <- fixture$cumulative_summary[
        fixture$cumulative_summary$bf_prior_id == bf_prior_id &
            fixture$cumulative_summary$schedule_id == schedule_id &
            fixture$cumulative_summary$evidence_threshold ==
                evidence_threshold,
        , drop = FALSE]
    curve <- curve[order(curve$look), , drop = FALSE]
    expect_equal(nrow(curve), final_row$n_looks[[1]],
                 info = paste(label, "selected sequential cumulative curve exists",
                              bf_prior_id, schedule_id, evidence_threshold))
    if (nrow(curve) != final_row$n_looks[[1]]) {
        return(invisible(final_row))
    }

    bf_prior <- bfpwr_sim_find_bf_prior_case(bf_prior_id, bf_priors)
    design <- bfpwr_sim_find_design_case(final_row$design_case_id[[1]],
                                         designs)
    reference <- bfpwr_fixture_z_sequential_reference(
        fixture = fixture,
        row = final_row,
        bf_prior = bf_prior,
        design = design,
        strict = strict)

    expect_true(
        all(bfpwr_sim_mc_close(curve$cum_pH1, reference$cumpH1, curve$nsim,
                               observed_mcse = curve$mcse_cum_pH1)),
        info = paste(label, "cumulative H1 curve agrees with pbf01seq",
                     bf_prior_id, schedule_id, evidence_threshold))
    expect_true(
        all(bfpwr_sim_mc_close(curve$cum_pH0, reference$cumpH0, curve$nsim,
                               observed_mcse = curve$mcse_cum_pH0)),
        info = paste(label, "cumulative H0 curve agrees with pbf01seq",
                     bf_prior_id, schedule_id, evidence_threshold))
    expect_true(
        all(bfpwr_sim_mc_close(curve$cum_pInc, reference$cumpInc, curve$nsim,
                               observed_mcse = curve$mcse_cum_pInc)),
        info = paste(label, "cumulative inconclusive curve agrees with pbf01seq",
                     bf_prior_id, schedule_id, evidence_threshold))
    expect_true(
        bfpwr_sim_mc_close(final_row$pH1, reference$cumpH1[[length(reference$cumpH1)]],
                           final_row$nsim, observed_mcse = final_row$mcse_pH1),
        info = paste(label, "final H1 probability agrees with pbf01seq",
                     bf_prior_id, schedule_id, evidence_threshold))
    expect_true(
        bfpwr_sim_mc_close(final_row$pH0, reference$cumpH0[[length(reference$cumpH0)]],
                           final_row$nsim, observed_mcse = final_row$mcse_pH0),
        info = paste(label, "final H0 probability agrees with pbf01seq",
                     bf_prior_id, schedule_id, evidence_threshold))
    expect_true(
        bfpwr_sim_mc_close(final_row$pInc, reference$cumpInc[[length(reference$cumpInc)]],
                           final_row$nsim, observed_mcse = final_row$mcse_pInc),
        info = paste(label, "final inconclusive probability agrees with pbf01seq",
                     bf_prior_id, schedule_id, evidence_threshold))
    expect_true(abs(final_row$EN[[1]] - reference$EN) <=
                    max(1, 4 * sqrt(max(final_row$VarN[[1]], 0) /
                                    final_row$nsim[[1]])),
                info = paste(label, "expected sample size agrees with pbf01seq",
                             bf_prior_id, schedule_id, evidence_threshold))
    invisible(reference)
}

bfpwr_fixture_fixed_check_grid <- function(representative_cases,
                                           n_by_grid = list(
                                               short = c(100L, 500L),
                                               long = c(100L, 1000L, 10000L)
                                           )) {
    do.call(rbind, lapply(seq_len(nrow(representative_cases)), function(i) {
        case <- representative_cases[i, , drop = FALSE]
        expand.grid(
            bf_prior_id = case$bf_prior_id,
            design_role = case$design_role,
            prior_form = case$prior_form,
            look_grid_name = case$look_grid_name,
            tail = c("H1", "H0"),
            evidence_threshold = c(3, 10, 30),
            n = n_by_grid[[case$look_grid_name]],
            KEEP.OUT.ATTRS = FALSE,
            stringsAsFactors = FALSE
        )
    }))
}

bfpwr_fixture_t_n2 <- function(n1, design) {
    if (identical(design$generation$type, "two.sample")) {
        pmax(2, round(n1 * design$generation$n2_multiplier))
    } else {
        n1
    }
}

bfpwr_fixture_t_fixed_reference <- function(row, designs, bf_priors) {
    if (!is.data.frame(row) || nrow(row) != 1L) {
        stop("t fixed reference requires exactly one summary row")
    }
    bf_prior <- bfpwr_sim_find_bf_prior_case(row$bf_prior_id[[1]], bf_priors)
    design <- bfpwr_sim_find_design_case(row$design_case_id[[1]], designs)
    prior <- bf_prior$analysis_prior
    dp <- bfpwr_sim_design_prior_mean_sd(design)
    n1 <- row$n[[1]]
    n2 <- bfpwr_fixture_t_n2(n1, design)

    ptbf01(
        k = row$threshold[[1]],
        n = n1,
        n1 = n1,
        n2 = n2,
        null = prior$null,
        plocation = prior$plocation,
        pscale = prior$pscale,
        pdf = prior$pdf,
        dpm = dp$dpm,
        dpsd = dp$dpsd,
        type = prior$type,
        alternative = prior$alternative,
        lower.tail = row$tail[[1]] == "H1")
}

bfpwr_fixture_t_sequential_reference <- function(fixture,
                                                 row,
                                                 bf_prior,
                                                 design,
                                                 strict = TRUE) {
    schedule <- bfpwr_fixture_sequential_schedule(fixture, row$schedule_id[[1]])
    n1 <- schedule$n
    n2 <- bfpwr_fixture_t_n2(n1, design)
    prior <- bf_prior$analysis_prior
    dp <- bfpwr_sim_design_prior_mean_sd(design)

    ptbf01seq(
        k1 = row$k1[[1]],
        k0 = row$k0[[1]],
        n = n1,
        n1 = n1,
        n2 = n2,
        plocation = prior$plocation - prior$null,
        pscale = prior$pscale,
        pdf = prior$pdf,
        dpm = dp$dpm - prior$null,
        dpsd = dp$dpsd,
        type = prior$type,
        alternative = prior$alternative,
        strict = strict)
}

bfpwr_fixture_t_sequential_reference_case <- function(fixture,
                                                      bf_prior_id,
                                                      schedule_id,
                                                      evidence_threshold,
                                                      designs,
                                                      bf_priors,
                                                      strict = TRUE) {
    final_row <- fixture$final_summary[
        fixture$final_summary$bf_prior_id == bf_prior_id &
            fixture$final_summary$schedule_id == schedule_id &
            fixture$final_summary$evidence_threshold == evidence_threshold,
        , drop = FALSE]
    if (nrow(final_row) != 1L) {
        stop("selected sequential final row is not unique: ",
             paste(bf_prior_id, schedule_id, evidence_threshold,
                   sep = ", "))
    }

    curve <- fixture$cumulative_summary[
        fixture$cumulative_summary$bf_prior_id == bf_prior_id &
            fixture$cumulative_summary$schedule_id == schedule_id &
            fixture$cumulative_summary$evidence_threshold ==
                evidence_threshold,
        , drop = FALSE]
    curve <- curve[order(curve$look), , drop = FALSE]
    if (nrow(curve) != final_row$n_looks[[1]]) {
        stop("selected sequential cumulative curve has wrong length: ",
             paste(bf_prior_id, schedule_id, evidence_threshold,
                   sep = ", "))
    }

    bf_prior <- bfpwr_sim_find_bf_prior_case(bf_prior_id, bf_priors)
    design <- bfpwr_sim_find_design_case(final_row$design_case_id[[1]],
                                         designs)
    reference <- bfpwr_fixture_t_sequential_reference(
        fixture = fixture,
        row = final_row,
        bf_prior = bf_prior,
        design = design,
        strict = strict)

    list(final_row = final_row, curve = curve, reference = reference)
}

bfpwr_fixture_binomial_fixed_reference <- function(row, designs, bf_priors) {
    if (!is.data.frame(row) || nrow(row) != 1L) {
        stop("binomial fixed reference requires exactly one summary row")
    }
    bf_prior <- bfpwr_sim_find_bf_prior_case(row$bf_prior_id[[1]], bf_priors)
    design <- bfpwr_sim_find_design_case(row$design_case_id[[1]], designs)
    prior <- bf_prior$analysis_prior
    design_args <- bfpwr_sim_binomial_design_args(design)
    pbinbf01(
        k = row$threshold[[1]],
        n = row$n[[1]],
        p0 = prior$p0,
        type = bf_prior$bf_type,
        a = prior$a,
        b = prior$b,
        dp = design_args$dp,
        da = design_args$da,
        db = design_args$db,
        dl = design_args$dl,
        du = design_args$du,
        lower.tail = row$tail[[1]] == "H1")
}
