library(tinytest)
library(bfpwr)

## Gated fixture-corpus test for fixed-sample binomial directional BF
## simulations. No bfssd manuscript formula covers this family; binary outcomes
## are future work in paper/bfssd.Rnw 1854-1856.

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

required_fixture_dir <- file.path("fixtures", "binomial", "fixed",
                                  "binom-binbf01-direction-fixed-core-v1")
fixture_context <- bfpwr_fixture_context(required_fixture_dir)
if (!is.null(fixture_context$skip)) {
    exit_file(fixture_context$skip)
}
on.exit(setwd(fixture_context$old_wd), add = TRUE)
corpus_root <- fixture_context$corpus_root

fixture <- bfpwr_sim_read_fixture_summary(
    corpus_root = corpus_root,
    fixture_set_id = "binom-binbf01-direction-fixed-core-v1",
    family = "binomial",
    mode = "fixed")

expect_true(nrow(fixture$tail_summary) > 0,
            info = "fixed binomial direction tail summary is available")
expect_true(nrow(fixture$decision_summary) > 0,
            info = "fixed binomial direction decision summary is available")
expect_true(nrow(fixture$search_summary) > 0,
            info = "fixed binomial direction search summary is available")

manifest_status <- bfpwr_fixture_manifest_status(fixture$fixture_dir)
expect_true(all(manifest_status$exists),
            info = "fixed binomial direction manifest files exist")
expect_true(all(manifest_status$bytes_ok),
            info = "fixed binomial direction manifest byte sizes match")
expect_true(all(manifest_status$sha256_ok),
            info = "fixed binomial direction manifest SHA256 hashes match")

fixture_failures <- bfpwr_sim_validate_fixed_fixture_deterministic(fixture)
expect_equal(nrow(fixture_failures), 0L,
             info = "fixed binomial direction deterministic invariants pass")

designs <- bfpwr_sim_design_case_set("production")
bf_priors <- bfpwr_sim_bf_prior_case_set("production")

representative_cases <- data.frame(
    bf_prior_id = c(
        "binom-binbf01-direction-p00p5-a1-b1-on-binom-dpoint-0p5-short",
        "binom-binbf01-direction-p00p5-a1-b1-on-binom-dpoint-0p6-short",
        "binom-binbf01-direction-p00p5-a1-b1-on-binom-dpoint-0p2-short",
        "binom-binbf01-direction-p00p2-a1-b1-on-binom-dpoint-0p2-short",
        "binom-binbf01-direction-p00p75-a1-b1-on-binom-dpoint-0p75-short",
        "binom-binbf01-direction-p00p5-a1-b1-on-binom-dbeta-1-1-h1gt0p5-short",
        "binom-binbf01-direction-p00p5-a1-b1-on-binom-dbeta-1-1-h0le0p5-short",
        "binom-binbf01-direction-p00p6-a1-b1-on-binom-dbeta-60-40-short",
        "binom-binbf01-direction-p00p5-a1-b1-on-binom-dpoint-0p5-long",
        "binom-binbf01-direction-p00p55-a1-b1-on-binom-dpoint-0p55-long"
    ),
    design_role = c(
        "null", "alternative", "adversarial", "null", "null",
        "alternative", "adversarial", "alternative", "null",
        "alternative"
    ),
    prior_form = c(
        "default", "default", "default", "matched-low", "matched-high",
        "directional-h1", "directional-h0", "matched-beta",
        "default", "matched-long"
    ),
    look_grid_name = c(
        "short", "short", "short", "short", "short",
        "short", "short", "short", "long", "long"
    ),
    stringsAsFactors = FALSE
)

expect_true(all(c("short", "long") %in% representative_cases$look_grid_name),
            info = "stratified checks include short and long grids")
expect_true(all(c("null", "alternative", "adversarial") %in%
                    representative_cases$design_role),
            info = "stratified checks include null, alternative, and adversarial designs")
expect_true(all(c("default", "matched-low", "matched-high",
                  "directional-h1", "directional-h0", "matched-beta",
                  "matched-long") %in% representative_cases$prior_form),
            info = "stratified checks include directional prior variants")

check_rows <- bfpwr_fixture_fixed_check_grid(
    representative_cases,
    n_by_grid = list(short = c(100L, 300L, 500L),
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
                 info = paste("selected fixed-summary row exists",
                              check_rows$bf_prior_id[[i]],
                              check_rows$tail[[i]],
                              check_rows$evidence_threshold[[i]],
                              check_rows$n[[i]]))
    if (nrow(row) != 1L) next

    reference <- bfpwr_fixture_binomial_fixed_reference(
        row, designs = designs, bf_priors = bf_priors)
    expect_true(
        bfpwr_sim_mc_close(row$prob, reference, row$nsim,
                           observed_mcse = row$mcse),
        info = paste("fixed-summary probability agrees with pbinbf01:",
                     row$bf_prior_id, row$tail,
                     row$evidence_threshold, row$n))
}

bfpwr_fixture_expect_fixed_search_first_crossing(
    fixture,
    bf_prior_id = "binom-binbf01-direction-p00p5-a1-b1-on-binom-dpoint-0p6-short",
    tail = "H1",
    evidence_threshold = 10,
    target_prob = 0.8,
    label = "fixed binomial direction fixture")

bfpwr_fixture_expect_fixed_search_first_crossing(
    fixture,
    bf_prior_id = "binom-binbf01-direction-p00p5-a1-b1-on-binom-dpoint-0p2-short",
    tail = "H0",
    evidence_threshold = 10,
    target_prob = 0.8,
    label = "fixed binomial direction fixture")
