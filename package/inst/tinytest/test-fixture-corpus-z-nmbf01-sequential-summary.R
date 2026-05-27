library(tinytest)
library(bfpwr)

## Gated fixture-corpus test for sequential normal-moment z BF simulations.
## Fixed-sample manuscript source is paper/bfssd.Rnw 1672-1709; the sequential
## fixture rows are package-level checks, not direct manuscript examples.

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
                                  "z-nmbf01-sequential-core-v1")
fixture_context <- bfpwr_fixture_context(required_fixture_dir)
if (!is.null(fixture_context$skip)) {
    exit_file(fixture_context$skip)
}
on.exit(setwd(fixture_context$old_wd), add = TRUE)
corpus_root <- fixture_context$corpus_root

fixture <- bfpwr_sim_read_fixture_summary(
    corpus_root = corpus_root,
    fixture_set_id = "z-nmbf01-sequential-core-v1",
    family = "z",
    mode = "sequential")

bfpwr_fixture_expect_sequential_summary_integrity(
    fixture,
    label = "sequential z-nmbf01 fixture")

designs <- bfpwr_sim_design_case_set("production")
bf_priors <- bfpwr_sim_bf_prior_case_set("production")

expect_equal(nrow(fixture$search_summary), 0L,
             info = "sequential z-nmbf01 search targets are disabled")
expect_true(all(c("start20-by10-looks02", "start20-by10-looks20",
                  "start10-by10-looks50", "start10-by10-looks100",
                  "start10-by10-looks200") %in%
                    unique(fixture$final_summary$schedule_id)),
            info = "sequential z-nmbf01 fixture includes all registered schedules")

representative_cases <- data.frame(
    bf_prior_id = c(
        "z-nmbf01-moment-null0-psd0p35355339-on-z-dpoint-0-usd1-short",
        "z-nmbf01-moment-null0-psd0p70710678-on-z-dpoint-0-usd1-short",
        "z-nmbf01-moment-null0-psd0p35355339-on-z-dpoint-0p5-usd1-short",
        "z-nmbf01-moment-null0p2-psd0p35355339-on-z-dpoint-0p2-usd1-short",
        "z-nmbf01-moment-null0-psd0p35355339-on-z-dpoint-m0p3-usd1-short",
        "z-nmbf01-moment-null0-psd0p70710678-on-z-dnorm-0p2-s0p5-usd1-long"
    ),
    design_role = c(
        "null", "null", "alternative", "alternative", "adversarial",
        "alternative"
    ),
    prior_form = c(
        "moment-narrow", "moment-wide", "moment-narrow", "moment-shifted",
        "moment-narrow", "moment-wide"
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
expect_true(any(representative_cases$prior_form == "moment-narrow") &&
                any(representative_cases$prior_form == "moment-wide"),
            info = "stratified checks include narrow and wide moment priors")

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
