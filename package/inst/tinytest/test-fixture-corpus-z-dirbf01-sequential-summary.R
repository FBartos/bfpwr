library(tinytest)
library(bfpwr)

## Gated fixture-corpus test for sequential directional z BF simulations. No
## manuscript row derives these directional-normal sequential cases; they are
## package-only checks of pbf01seq directional behavior.

helper_file <- system.file("tinytest", "helper-fixture-corpus.R",
                           package = "bfpwr")
helper_candidates <- c(
    helper_file,
    "helper-fixture-corpus.R",
    file.path("package", "inst", "tinytest", "helper-fixture-corpus.R")
)
helper_candidates <- helper_candidates[nzchar(helper_candidates) &
                                           file.exists(helper_candidates)]
if (length(helper_candidates) == 0) {
    stop("could not locate helper-fixture-corpus.R")
}
source(helper_candidates[[1]])

required_fixture_dir <- file.path("fixtures", "z", "sequential",
                                  "z-dirbf01-sequential-core-v1")
fixture_context <- bfpwr_fixture_context(required_fixture_dir)
if (!is.null(fixture_context$skip)) {
    exit_file(fixture_context$skip)
}
on.exit(setwd(fixture_context$old_wd), add = TRUE)
corpus_root <- fixture_context$corpus_root

fixture <- bfpwr_sim_read_fixture_summary(
    corpus_root = corpus_root,
    fixture_set_id = "z-dirbf01-sequential-core-v1",
    family = "z",
    mode = "sequential")

bfpwr_fixture_expect_sequential_summary_integrity(
    fixture,
    label = "sequential z-dirbf01 fixture")

designs <- bfpwr_sim_design_case_set("production")
bf_priors <- bfpwr_sim_bf_prior_case_set("production")

expect_equal(nrow(fixture$search_summary), 0L,
             info = "sequential z-dirbf01 search targets are disabled")
expect_true(all(c("start20-by10-looks02", "start20-by10-looks20",
                  "start10-by10-looks50", "start10-by10-looks100",
                  "start10-by10-looks200") %in%
                    unique(fixture$final_summary$schedule_id)),
            info = "sequential z-dirbf01 fixture includes all registered schedules")

representative_cases <- data.frame(
    bf_prior_id = c(
        "z-dirbf01-directional-normal-null0-pm0-psd0p70710678-on-z-dpoint-0-usd1-short",
        "z-dirbf01-directional-normal-null0-pm0-psd0p70710678-on-z-dpoint-0p2-usd1-short",
        "z-dirbf01-directional-normal-null0p2-pm0p2-psd0p70710678-on-z-dpoint-0p2-usd1-short",
        "z-dirbf01-directional-normal-null0-pm0-psd0p70710678-on-z-dpoint-m0p3-usd1-short",
        "z-dirbf01-directional-normal-null0-pm0-psd20-on-z-dpoint-10-usd50-short",
        "z-dirbf01-directional-normal-null0p2-pm0p2-psd0p70710678-on-z-dnorm-0p2-s0p5-usd1-long"
    ),
    design_role = c(
        "null", "alternative", "alternative", "adversarial",
        "adversarial", "alternative"
    ),
    prior_form = c(
        "directional-normal", "directional-normal",
        "shifted-directional-normal", "directional-normal",
        "extreme-directional-normal", "shifted-directional-normal"
    ),
    look_grid_name = c(
        "short", "short", "short", "short", "short", "long"
    ),
    schedule_id = rep("start20-by10-looks05", 6),
    stringsAsFactors = FALSE
)

expect_true(all(c("short", "long") %in% representative_cases$look_grid_name),
            info = "stratified checks include short and long grids")
expect_true(all(c("null", "alternative", "adversarial") %in%
                    representative_cases$design_role),
            info = "stratified checks include null, alternative, and adversarial designs")
expect_true(any(grepl("^shifted-", representative_cases$prior_form)),
            info = "stratified checks include shifted-null directional priors")

for (i in seq_len(nrow(representative_cases))) {
    for (threshold in c(3, 10, 30)) {
        bf_prior_id <- representative_cases$bf_prior_id[[i]]
        schedule_id <- representative_cases$schedule_id[[i]]
        final_row <- fixture$final_summary[
            fixture$final_summary$bf_prior_id == bf_prior_id &
                fixture$final_summary$schedule_id == schedule_id &
                fixture$final_summary$evidence_threshold == threshold,
            , drop = FALSE]
        expect_equal(nrow(final_row), 1L,
                     info = paste("selected sequential final row exists",
                                  bf_prior_id, schedule_id, threshold))
        if (nrow(final_row) != 1L) next

        curve <- fixture$cumulative_summary[
            fixture$cumulative_summary$bf_prior_id == bf_prior_id &
                fixture$cumulative_summary$schedule_id == schedule_id &
                fixture$cumulative_summary$evidence_threshold == threshold,
            , drop = FALSE]
        curve <- curve[order(curve$look), , drop = FALSE]
        expect_equal(nrow(curve), final_row$n_looks[[1]],
                     info = paste("selected sequential cumulative curve exists",
                                  bf_prior_id, schedule_id, threshold))
        if (nrow(curve) != final_row$n_looks[[1]]) next

        bf_prior <- bfpwr_sim_find_bf_prior_case(bf_prior_id, bf_priors)
        design <- bfpwr_sim_find_design_case(final_row$design_case_id[[1]],
                                             designs)
        reference <- bfpwr_fixture_z_sequential_reference(
            fixture = fixture,
            row = final_row,
            bf_prior = bf_prior,
            design = design)

        expect_true(
            all(bfpwr_sim_mc_close(curve$cum_pH1, reference$cumpH1,
                                   curve$nsim,
                                   observed_mcse = curve$mcse_cum_pH1)),
            info = paste("cumulative H1 curve agrees with pbf01seq",
                         bf_prior_id, schedule_id, threshold))
        expect_true(
            all(bfpwr_sim_mc_close(curve$cum_pH0, reference$cumpH0,
                                   curve$nsim,
                                   observed_mcse = curve$mcse_cum_pH0)),
            info = paste("cumulative H0 curve agrees with pbf01seq",
                         bf_prior_id, schedule_id, threshold))
        expect_true(
            all(bfpwr_sim_mc_close(curve$cum_pInc, reference$cumpInc,
                                   curve$nsim,
                                   observed_mcse = curve$mcse_cum_pInc)),
            info = paste("cumulative inconclusive curve agrees with pbf01seq",
                         bf_prior_id, schedule_id, threshold))
        expect_true(abs(final_row$EN[[1]] - reference$EN) <=
                        max(1, 4 * sqrt(max(final_row$VarN[[1]], 0) /
                                        final_row$nsim[[1]])),
                    info = paste("expected sample size agrees with pbf01seq",
                                 bf_prior_id, schedule_id, threshold))
    }
}

schedule_shape_cases <- data.frame(
    schedule_id = c("start20-by10-looks02",
                    "start20-by10-looks20",
                    "start10-by10-looks50"),
    evidence_threshold = c(10, 10, 10),
    strict = c(TRUE, TRUE, TRUE),
    stringsAsFactors = FALSE
)
for (i in seq_len(nrow(schedule_shape_cases))) {
    bf_prior_id <- "z-dirbf01-directional-normal-null0-pm0-psd0p70710678-on-z-dpoint-0-usd1-short"
    schedule_id <- schedule_shape_cases$schedule_id[[i]]
    threshold <- schedule_shape_cases$evidence_threshold[[i]]
    final_row <- fixture$final_summary[
        fixture$final_summary$bf_prior_id == bf_prior_id &
            fixture$final_summary$schedule_id == schedule_id &
            fixture$final_summary$evidence_threshold == threshold,
        , drop = FALSE]
    expect_equal(nrow(final_row), 1L,
                 info = paste("schedule-shape final row exists",
                              bf_prior_id, schedule_id, threshold))
    if (nrow(final_row) != 1L) next

    curve <- fixture$cumulative_summary[
        fixture$cumulative_summary$bf_prior_id == bf_prior_id &
            fixture$cumulative_summary$schedule_id == schedule_id &
            fixture$cumulative_summary$evidence_threshold == threshold,
        , drop = FALSE]
    curve <- curve[order(curve$look), , drop = FALSE]
    expect_equal(nrow(curve), final_row$n_looks[[1]],
                 info = paste("schedule-shape cumulative curve exists",
                              bf_prior_id, schedule_id, threshold))
    if (nrow(curve) != final_row$n_looks[[1]]) next

    bf_prior <- bfpwr_sim_find_bf_prior_case(bf_prior_id, bf_priors)
    design <- bfpwr_sim_find_design_case(final_row$design_case_id[[1]],
                                         designs)
    reference <- bfpwr_fixture_z_sequential_reference(
        fixture = fixture,
        row = final_row,
        bf_prior = bf_prior,
        design = design,
        strict = schedule_shape_cases$strict[[i]])

    expect_true(
        all(bfpwr_sim_mc_close(curve$cum_pH1, reference$cumpH1,
                               curve$nsim,
                               observed_mcse = curve$mcse_cum_pH1)),
        info = paste("schedule-shape cumulative H1 agrees with pbf01seq",
                     bf_prior_id, schedule_id, threshold))
    expect_true(
        all(bfpwr_sim_mc_close(curve$cum_pH0, reference$cumpH0,
                               curve$nsim,
                               observed_mcse = curve$mcse_cum_pH0)),
        info = paste("schedule-shape cumulative H0 agrees with pbf01seq",
                     bf_prior_id, schedule_id, threshold))
    expect_true(
        all(bfpwr_sim_mc_close(curve$cum_pInc, reference$cumpInc,
                               curve$nsim,
                               observed_mcse = curve$mcse_cum_pInc)),
        info = paste("schedule-shape cumulative inconclusive agrees with pbf01seq",
                     bf_prior_id, schedule_id, threshold))
    expect_true(abs(final_row$EN[[1]] - reference$EN) <=
                    max(1, 4 * sqrt(max(final_row$VarN[[1]], 0) /
                                    final_row$nsim[[1]])),
                info = paste("schedule-shape expected sample size agrees with pbf01seq",
                             bf_prior_id, schedule_id, threshold))
}

long_run_cases <- data.frame(
    bf_prior_id = c(
        "z-dirbf01-directional-normal-null0-pm0-psd0p70710678-on-z-dpoint-0p2-usd1-long",
        "z-dirbf01-directional-normal-null0p2-pm0p2-psd0p70710678-on-z-dnorm-0p2-s0p5-usd1-long"
    ),
    schedule_id = c("start10-by10-looks100",
                    "start10-by10-looks200"),
    evidence_threshold = c(10, 10),
    stringsAsFactors = FALSE
)
for (i in seq_len(nrow(long_run_cases))) {
    case <- bfpwr_fixture_z_sequential_reference_case(
        fixture = fixture,
        bf_prior_id = long_run_cases$bf_prior_id[[i]],
        schedule_id = long_run_cases$schedule_id[[i]],
        evidence_threshold = long_run_cases$evidence_threshold[[i]],
        designs = designs,
        bf_priors = bf_priors,
        strict = TRUE)
    expect_true(
        all(bfpwr_sim_mc_close(case$curve$cum_pH1, case$reference$cumpH1,
                               case$curve$nsim,
                               observed_mcse = case$curve$mcse_cum_pH1)),
        info = paste("long-run cumulative H1 agrees with pbf01seq",
                     long_run_cases$bf_prior_id[[i]],
                     long_run_cases$schedule_id[[i]]))
    expect_true(
        all(bfpwr_sim_mc_close(case$curve$cum_pH0, case$reference$cumpH0,
                               case$curve$nsim,
                               observed_mcse = case$curve$mcse_cum_pH0)),
        info = paste("long-run cumulative H0 agrees with pbf01seq",
                     long_run_cases$bf_prior_id[[i]],
                     long_run_cases$schedule_id[[i]]))
    expect_true(
        all(bfpwr_sim_mc_close(case$curve$cum_pInc, case$reference$cumpInc,
                               case$curve$nsim,
                               observed_mcse = case$curve$mcse_cum_pInc)),
        info = paste("long-run cumulative inconclusive agrees with pbf01seq",
                     long_run_cases$bf_prior_id[[i]],
                     long_run_cases$schedule_id[[i]]))
    expect_true(abs(case$final_row$EN[[1]] - case$reference$EN) <=
                    max(1, 4 * sqrt(max(case$final_row$VarN[[1]], 0) /
                                    case$final_row$nsim[[1]])),
                info = paste("long-run expected sample size agrees with pbf01seq",
                             long_run_cases$bf_prior_id[[i]],
                             long_run_cases$schedule_id[[i]]))
}
