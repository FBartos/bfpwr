library(tinytest)
library(bfpwr)

## Gated fixture-corpus test for fixed-sample t BF simulations. It compares
## cached MC summaries with tbf01/ptbf01 references. Manuscript source:
## paper/bfssd.Rnw 1481-1530 and 1609-1627; fixture rows are regressions.

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

required_fixture_dir <- file.path("fixtures", "t", "fixed",
                                  "t-tbf01-fixed-core-v1")
fixture_context <- bfpwr_fixture_context(required_fixture_dir)
if (!is.null(fixture_context$skip)) {
    exit_file(fixture_context$skip)
}
on.exit(setwd(fixture_context$old_wd), add = TRUE)
corpus_root <- fixture_context$corpus_root

fixture <- bfpwr_sim_read_fixture_summary(
    corpus_root = corpus_root,
    fixture_set_id = "t-tbf01-fixed-core-v1",
    family = "t",
    mode = "fixed")

expect_true(nrow(fixture$tail_summary) > 0,
            info = "fixed t-tbf01 tail summary is available")
expect_true(nrow(fixture$decision_summary) > 0,
            info = "fixed t-tbf01 decision summary is available")
expect_true(nrow(fixture$search_summary) > 0,
            info = "fixed t-tbf01 search summary is available")

manifest_status <- bfpwr_fixture_manifest_status(fixture$fixture_dir)
expect_true(all(manifest_status$exists),
            info = "fixed t-tbf01 manifest files exist")
expect_true(all(manifest_status$bytes_ok),
            info = "fixed t-tbf01 manifest byte sizes match")
expect_true(all(manifest_status$sha256_ok),
            info = "fixed t-tbf01 manifest SHA256 hashes match")

fixture_failures <- bfpwr_sim_validate_fixed_fixture_deterministic(fixture)
expect_equal(nrow(fixture_failures), 0L,
             info = "fixed t-tbf01 deterministic invariants pass")

designs <- bfpwr_sim_design_case_set("production")
bf_priors <- bfpwr_sim_bf_prior_case_set("production")

representative_cases <- data.frame(
    bf_prior_id = c(
        "t-tbf01-cauchy-null0-loc0-scale0p70710678-df1-two-sided-on-t-two-dpoint-0-short",
        "t-tbf01-cauchy-null0-loc0-scale0p70710678-df1-greater-on-t-two-dpoint-0p5-short",
        "t-tbf01-cauchy-null0-loc0-scale0p70710678-df1-less-on-t-two-dpoint-0p5-short",
        "t-tbf01-student-t-null0-loc0-scale0p70710678-df3-two-sided-on-t-two-dpoint-0p5-short",
        "t-tbf01-student-t-null0-loc0-scale0p70710678-df30-two-sided-on-t-two-dpoint-0-short",
        "t-tbf01-cauchy-null0-loc0-scale0p70710678-df1-two-sided-on-t-one-dpoint-0-short",
        "t-tbf01-student-t-null0-loc0p5-scale0p1-df3-greater-on-t-one-dpoint-0p5-short",
        "t-tbf01-cauchy-null0-loc0-scale0p70710678-df1-two-sided-on-t-paired-dpoint-0p5-short",
        "t-tbf01-student-t-null0p2-loc0p5-scale0p35-df30-greater-on-t-paired-dpoint-0-short",
        "t-tbf01-cauchy-null0-loc0-scale0p70710678-df1-greater-on-t-two-dpoint-m0p4-short",
        "t-tbf01-student-t-null0p2-loc0p5-scale0p35-df30-greater-on-t-two-dpoint-0-n2x1p1-short",
        "t-tbf01-student-t-null0p2-loc0p5-scale0p35-df3-greater-on-t-two-dpoint-0p35-n2x1p1-short",
        "t-tbf01-cauchy-null0-loc0-scale0p70710678-df1-two-sided-on-t-two-dpoint-0-long",
        "t-tbf01-student-t-null0-loc0p5-scale0p1-df3-greater-on-t-two-dpoint-0p2-long",
        "t-tbf01-student-t-null0p2-loc0p5-scale0p35-df3-greater-on-t-two-dnorm-0p2-s0p5-long",
        "t-tbf01-student-t-null0-loc0-scale0p70710678-df30-two-sided-on-t-two-dnorm-0p2-s0p5-long"
    ),
    design_role = c(
        "null", "alternative", "alternative", "alternative", "null",
        "null", "alternative", "alternative", "null", "adversarial",
        "null", "alternative", "null", "alternative", "adversarial",
        "adversarial"
    ),
    sampling = c(
        "two.sample", "two.sample", "two.sample", "two.sample",
        "two.sample", "one.sample", "one.sample", "paired", "paired",
        "two.sample", "two.sample-unequal", "two.sample-unequal",
        "two.sample", "two.sample", "two.sample", "two.sample"
    ),
    prior_form = c(
        "cauchy", "cauchy", "cauchy", "student-t-df3",
        "student-t-df30", "cauchy", "student-t-local",
        "cauchy", "shifted-null-df30", "wrong-direction",
        "shifted-null-df30", "shifted-null-df3", "cauchy",
        "student-t-local", "shifted-null-df3", "student-t-df30"
    ),
    alternative = c(
        "two.sided", "greater", "less", "two.sided", "two.sided",
        "two.sided", "greater", "two.sided", "greater", "greater",
        "greater", "greater", "two.sided", "greater", "greater",
        "two.sided"
    ),
    look_grid_name = c(
        rep("short", 12), rep("long", 4)
    ),
    stringsAsFactors = FALSE
)

expect_true(all(c("short", "long") %in% representative_cases$look_grid_name),
            info = "stratified checks include short and long grids")
expect_true(all(c("null", "alternative", "adversarial") %in%
                    representative_cases$design_role),
            info = "stratified checks include null, alternative, and adversarial designs")
expect_true(all(c("two.sample", "one.sample", "paired",
                  "two.sample-unequal") %in% representative_cases$sampling),
            info = "stratified checks include all t-test sampling paths")
expect_true(all(c("cauchy", "student-t-df3", "student-t-df30",
                  "student-t-local", "shifted-null-df3",
                  "shifted-null-df30") %in%
                    representative_cases$prior_form),
            info = "stratified checks include Cauchy, Student-t, localized, and shifted-null priors")
expect_true(all(c("two.sided", "greater", "less") %in%
                    representative_cases$alternative),
            info = "stratified checks include two-sided, greater, and less alternatives")

check_rows <- bfpwr_fixture_fixed_check_grid(
    representative_cases,
    n_by_grid = list(short = c(100L, 500L),
                     long = c(100L, 1000L, 10000L)))

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
                 info = paste("selected fixed t row exists",
                              check_rows$bf_prior_id[[i]],
                              check_rows$tail[[i]],
                              check_rows$evidence_threshold[[i]],
                              check_rows$n[[i]]))
    if (nrow(row) != 1L) next

    reference <- bfpwr_fixture_t_fixed_reference(
        row, designs = designs, bf_priors = bf_priors)
    expect_true(
        bfpwr_sim_mc_close(row$prob, reference, row$nsim,
                           observed_mcse = row$mcse),
        info = paste("fixed t probability agrees with ptbf01:",
                     row$bf_prior_id, row$tail,
                     row$evidence_threshold, row$n))
}

bfpwr_fixture_expect_fixed_search_first_crossing(
    fixture,
    bf_prior_id = "t-tbf01-cauchy-null0-loc0-scale0p70710678-df1-greater-on-t-two-dpoint-0p5-short",
    tail = "H1",
    evidence_threshold = 10,
    target_prob = 0.8,
    label = "fixed t-tbf01 fixture")

bfpwr_fixture_expect_fixed_search_first_crossing(
    fixture,
    bf_prior_id = "t-tbf01-cauchy-null0-loc0-scale0p70710678-df1-two-sided-on-t-two-dpoint-0-short",
    tail = "H0",
    evidence_threshold = 3,
    target_prob = 0.3,
    label = "fixed t-tbf01 fixture")
