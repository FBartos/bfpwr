library(tinytest)
library(bfpwr)

## Gated fixture-corpus test for sequential t BF simulations. It compares cached
## MC summaries with ptbf01seq references. Related source: BFGSD appendix
## one-sided JZS design and bfssd.Rnw 1481-1530; fixture rows are regressions.

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

required_fixture_dir <- file.path("fixtures", "t", "sequential",
                                  "t-tbf01-sequential-core-v1")
fixture_context <- bfpwr_fixture_context(required_fixture_dir)
if (!is.null(fixture_context$skip)) {
    exit_file(fixture_context$skip)
}
on.exit(setwd(fixture_context$old_wd), add = TRUE)
corpus_root <- fixture_context$corpus_root

fixture <- bfpwr_sim_read_fixture_summary(
    corpus_root = corpus_root,
    fixture_set_id = "t-tbf01-sequential-core-v1",
    family = "t",
    mode = "sequential")

expect_true(nrow(fixture$cumulative_summary) > 0,
            info = "sequential t-tbf01 cumulative summary is available")
expect_true(nrow(fixture$final_summary) > 0,
            info = "sequential t-tbf01 final summary is available")
expect_equal(nrow(fixture$search_summary), 0L,
             info = "sequential t-tbf01 search summary is intentionally empty")

manifest_status <- bfpwr_fixture_manifest_status(fixture$fixture_dir)
expect_true(all(manifest_status$exists),
            info = "sequential t-tbf01 manifest files exist")
expect_true(all(manifest_status$bytes_ok),
            info = "sequential t-tbf01 manifest byte sizes match")
expect_true(all(manifest_status$sha256_ok),
            info = "sequential t-tbf01 manifest SHA256 hashes match")

fixture_failures <- bfpwr_sim_validate_sequential_fixture_deterministic(fixture)
expect_equal(nrow(fixture_failures), 0L,
             info = "sequential t-tbf01 deterministic invariants pass")

designs <- bfpwr_sim_design_case_set("production")
bf_priors <- bfpwr_sim_bf_prior_case_set("production")

expect_equal(nrow(fixture$search_summary), 0L,
             info = "sequential t-tbf01 search targets are disabled")
registered_schedule_ids <- vapply(fixture$spec$schedules,
                                   function(x) x$schedule_id,
                                   character(1))
expect_true(all(registered_schedule_ids %in%
                    unique(fixture$final_summary$schedule_id)),
            info = "sequential t-tbf01 fixture includes all registered schedules")

representative_cases <- data.frame(
    bf_prior_id = c(
        "t-tbf01-cauchy-null0-loc0-scale0p70710678-df1-two-sided-on-t-two-dpoint-0-short",
        "t-tbf01-cauchy-null0-loc0-scale0p70710678-df1-greater-on-t-two-dpoint-0p5-short",
        "t-tbf01-cauchy-null0-loc0-scale0p70710678-df1-less-on-t-two-dpoint-0p5-short",
        "t-tbf01-student-t-null0-loc0p5-scale0p1-df3-greater-on-t-one-dpoint-0p5-short",
        "t-tbf01-cauchy-null0-loc0-scale0p70710678-df1-two-sided-on-t-paired-dpoint-0p5-short",
        "t-tbf01-student-t-null0p2-loc0p5-scale0p35-df3-greater-on-t-two-dpoint-0p35-n2x1p1-short",
        "t-tbf01-cauchy-null0-loc0-scale0p70710678-df1-greater-on-t-two-dpoint-m0p4-short",
        "t-tbf01-cauchy-null0-loc0-scale0p70710678-df1-two-sided-on-t-two-dpoint-0-long",
        "t-tbf01-student-t-null0p2-loc0p5-scale0p35-df3-greater-on-t-two-dpoint-0p2-long"
    ),
    design_role = c(
        "null", "alternative", "alternative", "alternative",
        "alternative", "alternative", "adversarial", "null",
        "alternative"
    ),
    sampling = c(
        "two.sample", "two.sample", "two.sample", "one.sample",
        "paired", "two.sample-unequal", "two.sample", "two.sample",
        "two.sample"
    ),
    prior_form = c(
        "cauchy", "cauchy", "cauchy", "student-t-local",
        "cauchy", "shifted-null-df3", "wrong-direction",
        "cauchy", "shifted-null-df3"
    ),
    alternative = c(
        "two.sided", "greater", "less", "greater", "two.sided",
        "greater", "greater", "two.sided", "greater"
    ),
    look_grid_name = c(
        rep("short", 7), "long", "long"
    ),
    schedule_id = rep("start20-by10-looks05", 9),
    stringsAsFactors = FALSE
)

expect_true(all(c("short", "long") %in% representative_cases$look_grid_name),
            info = "stratified sequential checks include short and long grids")
expect_true(all(c("null", "alternative", "adversarial") %in%
                    representative_cases$design_role),
            info = "stratified sequential checks include null, alternative, and adversarial designs")
expect_true(all(c("two.sample", "one.sample", "paired",
                  "two.sample-unequal") %in% representative_cases$sampling),
            info = "stratified sequential checks include all t-test sampling paths")
expect_true(all(c("two.sided", "greater", "less") %in%
                    representative_cases$alternative),
            info = "stratified sequential checks include two-sided, greater, and less alternatives")

for (i in seq_len(nrow(representative_cases))) {
    for (threshold in c(3, 10, 30)) {
        final_row <- fixture$final_summary[
            fixture$final_summary$bf_prior_id ==
                representative_cases$bf_prior_id[[i]] &
                fixture$final_summary$schedule_id ==
                    representative_cases$schedule_id[[i]] &
                fixture$final_summary$evidence_threshold == threshold,
            , drop = FALSE]
        expect_equal(nrow(final_row), 1L,
                     info = paste("selected sequential t final row exists",
                                  representative_cases$bf_prior_id[[i]],
                                  representative_cases$schedule_id[[i]],
                                  threshold))
        if (nrow(final_row) != 1L) next

        curve <- fixture$cumulative_summary[
            fixture$cumulative_summary$bf_prior_id ==
                representative_cases$bf_prior_id[[i]] &
                fixture$cumulative_summary$schedule_id ==
                    representative_cases$schedule_id[[i]] &
                fixture$cumulative_summary$evidence_threshold == threshold,
            , drop = FALSE]
        curve <- curve[order(curve$look), , drop = FALSE]
        expect_equal(nrow(curve), final_row$n_looks[[1]],
                     info = paste("selected sequential t cumulative curve exists",
                                  representative_cases$bf_prior_id[[i]],
                                  representative_cases$schedule_id[[i]],
                                  threshold))
    }
}

strict_reference_cases <- data.frame(
    bf_prior_id = c(
        "t-tbf01-cauchy-null0-loc0-scale0p70710678-df1-greater-on-t-two-dpoint-0p5-short",
        "t-tbf01-cauchy-null0-loc0-scale0p70710678-df1-less-on-t-two-dpoint-0p5-short",
        "t-tbf01-student-t-null0-loc0p5-scale0p1-df3-greater-on-t-one-dpoint-0p5-short",
        "t-tbf01-cauchy-null0-loc0-scale0p70710678-df1-greater-on-t-two-dpoint-m0p4-short"
    ),
    schedule_id = rep("start20-by10-looks20", 4),
    evidence_threshold = c(10, 10, 10, 3),
    stringsAsFactors = FALSE
)

for (i in seq_len(nrow(strict_reference_cases))) {
    case <- bfpwr_fixture_t_sequential_reference_case(
        fixture = fixture,
        bf_prior_id = strict_reference_cases$bf_prior_id[[i]],
        schedule_id = strict_reference_cases$schedule_id[[i]],
        evidence_threshold = strict_reference_cases$evidence_threshold[[i]],
        designs = designs,
        bf_priors = bf_priors,
        strict = TRUE)

    expect_equal(case$final_row$n_looks[[1]], 20L,
                 info = paste("strict t reference case uses 20 looks",
                              strict_reference_cases$bf_prior_id[[i]]))
    expect_true(
        all(bfpwr_sim_mc_close(case$curve$cum_pH1,
                               case$reference$cumpH1,
                               case$curve$nsim,
                               observed_mcse = case$curve$mcse_cum_pH1,
                               z = 5)),
        info = paste("20-look strict cumulative H1 agrees with ptbf01seq",
                     strict_reference_cases$bf_prior_id[[i]],
                     strict_reference_cases$evidence_threshold[[i]]))
    expect_true(
        all(bfpwr_sim_mc_close(case$curve$cum_pH0,
                               case$reference$cumpH0,
                               case$curve$nsim,
                               observed_mcse = case$curve$mcse_cum_pH0,
                               z = 5)),
        info = paste("20-look strict cumulative H0 agrees with ptbf01seq",
                     strict_reference_cases$bf_prior_id[[i]],
                     strict_reference_cases$evidence_threshold[[i]]))
    expect_true(
        all(bfpwr_sim_mc_close(case$curve$cum_pInc,
                               case$reference$cumpInc,
                               case$curve$nsim,
                               observed_mcse = case$curve$mcse_cum_pInc,
                               z = 5)),
        info = paste("20-look strict cumulative inconclusive agrees with ptbf01seq",
                     strict_reference_cases$bf_prior_id[[i]],
                     strict_reference_cases$evidence_threshold[[i]]))
    expect_true(abs(case$final_row$EN[[1]] - case$reference$EN1) <=
                    max(1, 4 * sqrt(max(case$final_row$VarN[[1]], 0) /
                                    case$final_row$nsim[[1]])),
                info = paste("20-look strict expected sample size agrees with ptbf01seq",
                             strict_reference_cases$bf_prior_id[[i]],
                             strict_reference_cases$evidence_threshold[[i]]))
}

stress_cases <- data.frame(
    bf_prior_id = c(
        "t-tbf01-cauchy-null0-loc0-scale0p70710678-df1-greater-on-t-two-dpoint-m0p4-short",
        "t-tbf01-cauchy-null0-loc0-scale0p70710678-df1-two-sided-on-t-two-dpoint-0-long",
        "t-tbf01-student-t-null0p2-loc0p5-scale0p35-df3-greater-on-t-two-dpoint-0p2-long"
    ),
    schedule_id = c("start10-by10-looks50",
                    "start10-by10-looks100",
                    "start10-by10-looks200"),
    evidence_threshold = c(10, 10, 10),
    stringsAsFactors = FALSE
)

for (i in seq_len(nrow(stress_cases))) {
    final_row <- fixture$final_summary[
        fixture$final_summary$bf_prior_id == stress_cases$bf_prior_id[[i]] &
            fixture$final_summary$schedule_id == stress_cases$schedule_id[[i]] &
            fixture$final_summary$evidence_threshold ==
                stress_cases$evidence_threshold[[i]],
        , drop = FALSE]
    expect_equal(nrow(final_row), 1L,
                 info = paste("selected t stress final row exists",
                              stress_cases$bf_prior_id[[i]],
                              stress_cases$schedule_id[[i]],
                              stress_cases$evidence_threshold[[i]]))
    if (nrow(final_row) != 1L) next

    curve <- fixture$cumulative_summary[
        fixture$cumulative_summary$bf_prior_id == stress_cases$bf_prior_id[[i]] &
            fixture$cumulative_summary$schedule_id == stress_cases$schedule_id[[i]] &
            fixture$cumulative_summary$evidence_threshold ==
                stress_cases$evidence_threshold[[i]],
        , drop = FALSE]
    expect_equal(nrow(curve), final_row$n_looks[[1]],
                 info = paste("selected t stress cumulative curve exists",
                              stress_cases$bf_prior_id[[i]],
                              stress_cases$schedule_id[[i]],
                              stress_cases$evidence_threshold[[i]]))
}
