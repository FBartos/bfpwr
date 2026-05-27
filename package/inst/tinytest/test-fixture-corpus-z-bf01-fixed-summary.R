library(tinytest)
library(bfpwr)

## Gated fixture-corpus test for fixed-sample z BF01 simulations. It compares
## cached Monte Carlo summaries with pbf01 references. Manuscript source:
## paper/bfssd.Rnw 354-369 and 449-606; fixture rows are not manuscript examples.

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

required_fixture_dir <- file.path("fixtures", "z", "fixed",
                                  "z-bf01-fixed-core-v1")
fixture_context <- bfpwr_fixture_context(required_fixture_dir)
if (!is.null(fixture_context$skip)) {
    exit_file(fixture_context$skip)
}
on.exit(setwd(fixture_context$old_wd), add = TRUE)
corpus_root <- fixture_context$corpus_root

fixture <- bfpwr_sim_read_fixture_summary(
    corpus_root = corpus_root,
    fixture_set_id = "z-bf01-fixed-core-v1",
    family = "z",
    mode = "fixed")

bfpwr_fixture_expect_fixed_summary_integrity(
    fixture,
    label = "fixed z-bf01 fixture")

designs <- bfpwr_sim_design_case_set("production")
bf_priors <- bfpwr_sim_bf_prior_case_set("production")

representative_cases <- data.frame(
    bf_prior_id = c(
        "z-bf01-point-null0-pm0p2-psd0-on-z-dpoint-0-usd1-short",
        "z-bf01-normal-null0-pm0-psd0p70710678-on-z-dpoint-0-usd1-short",
        "z-bf01-point-null0-pm0p5-psd0-on-z-dpoint-0p5-usd1-short",
        "z-bf01-normal-null0-pm0p2-psd0p3-on-z-dpoint-0p5-usd1-short",
        "z-bf01-point-null0-pm0p2-psd0-on-z-dpoint-m0p3-usd1-short",
        "z-bf01-normal-null0-pm0-psd0p70710678-on-z-dnorm-0-s0p5-usd1-short",
        "z-bf01-point-null0-pm10-psd0-on-z-dpoint-10-usd50-short",
        "z-bf01-normal-null0-pm10-psd20-on-z-dnorm-10-s20-usd50-short",
        "z-bf01-point-null0-pm0p2-psd0-on-z-dpoint-0-usd1-long",
        "z-bf01-normal-null0-pm0-psd0p70710678-on-z-dpoint-0-usd1-long",
        "z-bf01-point-null0-pm0p2-psd0-on-z-dpoint-0p2-usd1-long",
        "z-bf01-normal-null0-pm0p2-psd0p3-on-z-dnorm-0p2-s0p5-usd1-long"
    ),
    design_role = c(
        "null", "null", "alternative", "alternative",
        "adversarial", "adversarial", "adversarial", "adversarial",
        "null", "null", "alternative", "alternative"
    ),
    prior_form = c(
        "point", "normal", "point", "normal",
        "point", "normal", "point", "normal",
        "point", "normal", "point", "normal"
    ),
    look_grid_name = c(
        "short", "short", "short", "short",
        "short", "short", "short", "short",
        "long", "long", "long", "long"
    ),
    stringsAsFactors = FALSE
)

expect_true(all(c("short", "long") %in% representative_cases$look_grid_name),
            info = "stratified checks include short and long grids")
expect_true(all(c("null", "alternative", "adversarial") %in%
                    representative_cases$design_role),
            info = "stratified checks include null, alternative, and adversarial designs")
expect_true(all(c("point", "normal") %in% representative_cases$prior_form),
            info = "stratified checks include point and normal analysis priors")

check_rows <- bfpwr_fixture_fixed_check_grid(representative_cases)

expect_true(all(c(3, 10, 30) %in% check_rows$evidence_threshold),
            info = "stratified checks include BF thresholds 3, 10, and 30")
expect_true(all(c("H1", "H0") %in% check_rows$tail),
            info = "stratified checks include H1 and H0 evidence tails")

for (i in seq_len(nrow(check_rows))) {
    row <- fixture$tail_summary[
        fixture$tail_summary$bf_prior_id == check_rows$bf_prior_id[[i]] &
            fixture$tail_summary$tail == check_rows$tail[[i]] &
            fixture$tail_summary$evidence_threshold ==
                check_rows$evidence_threshold[[i]] &
            fixture$tail_summary$n == check_rows$n[[i]],
        , drop = FALSE]
    expect_equal(nrow(row), 1L,
                 info = paste("selected fixed-summary row exists",
                              check_rows$bf_prior_id[[i]],
                              check_rows$tail[[i]],
                              check_rows$evidence_threshold[[i]],
                              check_rows$n[[i]]))

    bf_prior <- bfpwr_sim_find_bf_prior_case(row$bf_prior_id, bf_priors)
    design <- bfpwr_sim_find_design_case(row$design_case_id, designs)
    dp <- bfpwr_sim_design_prior_mean_sd(design)
    prior <- bf_prior$analysis_prior
    reference <- pbf01(
        k = row$threshold,
        n = row$n,
        usd = design$generation$usd,
        null = prior$null,
        pm = prior$pm,
        psd = prior$psd,
        dpm = dp$dpm,
        dpsd = dp$dpsd,
        lower.tail = row$tail == "H1")

    expect_true(
        bfpwr_sim_mc_close(row$prob, reference, row$nsim,
                           observed_mcse = row$mcse),
        info = paste("fixed-summary probability agrees with pbf01:",
                     row$bf_prior_id, row$tail,
                     row$evidence_threshold, row$n))
}

bfpwr_fixture_expect_fixed_search_first_crossing(
    fixture,
    bf_prior_id = "z-bf01-point-null0-pm0p2-psd0-on-z-dpoint-0-usd1-short",
    tail = "H0",
    evidence_threshold = 10,
    target_prob = 0.8,
    label = "fixed z-bf01 fixture")
