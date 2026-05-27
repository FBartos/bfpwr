library(tinytest)
library(bfpwr)

## Gated fixture-corpus test for sequential binomial point-null BF simulations.
## There is no manuscript source for binomial sequential designs; this checks
## package-only cached summaries and schedule coverage.

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

required_fixture_dir <- file.path("fixtures", "binomial", "sequential",
                                  "binom-binbf01-point-sequential-core-v1")
fixture_context <- bfpwr_fixture_context(required_fixture_dir)
if (!is.null(fixture_context$skip)) {
    exit_file(fixture_context$skip)
}
on.exit(setwd(fixture_context$old_wd), add = TRUE)
corpus_root <- fixture_context$corpus_root

fixture <- bfpwr_sim_read_fixture_summary(
    corpus_root = corpus_root,
    fixture_set_id = "binom-binbf01-point-sequential-core-v1",
    family = "binomial",
    mode = "sequential")

expect_true(nrow(fixture$cumulative_summary) > 0,
            info = "sequential binomial point cumulative summary is available")
expect_true(nrow(fixture$final_summary) > 0,
            info = "sequential binomial point final summary is available")
expect_equal(nrow(fixture$search_summary), 0L,
             info = "sequential binomial point search summary is intentionally empty")

manifest_status <- bfpwr_fixture_manifest_status(fixture$fixture_dir)
expect_true(all(manifest_status$exists),
            info = "sequential binomial point manifest files exist")
expect_true(all(manifest_status$bytes_ok),
            info = "sequential binomial point manifest byte sizes match")
expect_true(all(manifest_status$sha256_ok),
            info = "sequential binomial point manifest SHA256 hashes match")

fixture_failures <- bfpwr_sim_validate_sequential_fixture_deterministic(fixture)
expect_equal(nrow(fixture_failures), 0L,
             info = "sequential binomial point deterministic invariants pass")

expect_equal(nrow(fixture$search_summary), 0L,
             info = "sequential binomial point search targets are disabled")
expect_true(all(c("start20-by10-looks02", "start20-by10-looks20",
                  "start10-by10-looks50", "start10-by10-looks100",
                  "start10-by10-looks200") %in%
                    unique(fixture$final_summary$schedule_id)),
            info = "sequential binomial point fixture includes all registered schedules")

representative_cases <- data.frame(
    bf_prior_id = c(
        "binom-binbf01-point-p00p5-a1-b1-on-binom-dpoint-0p5-short",
        "binom-binbf01-point-p00p5-a1-b1-on-binom-dpoint-0p6-short",
        "binom-binbf01-point-p00p5-a60-b40-on-binom-dbeta-60-40-short",
        "binom-binbf01-point-p00p5-a40-b60-on-binom-dbeta-1-1-h0le0p5-short"
    ),
    design_role = c("null", "alternative", "alternative", "adversarial"),
    prior_variant = c("default", "default", "concentrated-h1",
                      "concentrated-h0"),
    look_grid_name = rep("short", 4),
    schedule_id = rep("start20-by10-looks05", 4),
    thresholds = I(rep(list(c(3, 10, 30)), 4)),
    stringsAsFactors = FALSE
)

for (i in seq_len(nrow(representative_cases))) {
    rows <- fixture$final_summary[
        fixture$final_summary$bf_prior_id == representative_cases$bf_prior_id[[i]] &
            fixture$final_summary$schedule_id == representative_cases$schedule_id[[i]],
        , drop = FALSE]
    expect_true(nrow(rows) > 0,
                info = paste("selected sequential point summary rows exist",
                             representative_cases$bf_prior_id[[i]],
                             representative_cases$schedule_id[[i]]))
    if (nrow(rows) == 0) next
    expect_true(all(representative_cases$thresholds[[i]] %in%
                        rows$evidence_threshold),
                info = paste("selected sequential point thresholds exist",
                             representative_cases$bf_prior_id[[i]],
                             representative_cases$schedule_id[[i]]))
}

expect_true(all(c("null", "alternative", "adversarial") %in%
                    representative_cases$design_role),
            info = "stratified checks include null, alternative, and adversarial designs")
expect_true(all(c("default", "concentrated-h1", "concentrated-h0") %in%
                    representative_cases$prior_variant),
            info = "stratified checks include default and concentrated priors")

schedule_shape_cases <- data.frame(
    bf_prior_id = "binom-binbf01-point-p00p5-a1-b1-on-binom-dpoint-0p5-short",
    schedule_id = c("start20-by10-looks02",
                    "start20-by10-looks20",
                    "start10-by10-looks50"),
    evidence_threshold = c(10, 10, 10),
    stringsAsFactors = FALSE
)

long_run_cases <- data.frame(
    bf_prior_id = c(
        "binom-binbf01-point-p00p5-a1-b1-on-binom-dpoint-0p5-long",
        "binom-binbf01-point-p00p5-a5100-b4900-on-binom-dpoint-0p55-long"
    ),
    schedule_id = c("start10-by10-looks100",
                    "start10-by10-looks200"),
    evidence_threshold = c(10, 10),
    stringsAsFactors = FALSE
)

coverage_rows <- rbind(schedule_shape_cases, long_run_cases)
for (i in seq_len(nrow(coverage_rows))) {
    final_row <- fixture$final_summary[
        fixture$final_summary$bf_prior_id == coverage_rows$bf_prior_id[[i]] &
            fixture$final_summary$schedule_id == coverage_rows$schedule_id[[i]] &
            fixture$final_summary$evidence_threshold ==
                coverage_rows$evidence_threshold[[i]],
        , drop = FALSE]
    expect_equal(nrow(final_row), 1L,
                 info = paste("selected sequential point coverage row exists",
                              coverage_rows$bf_prior_id[[i]],
                              coverage_rows$schedule_id[[i]],
                              coverage_rows$evidence_threshold[[i]]))
    if (nrow(final_row) != 1L) next

    curve <- fixture$cumulative_summary[
        fixture$cumulative_summary$bf_prior_id == coverage_rows$bf_prior_id[[i]] &
            fixture$cumulative_summary$schedule_id == coverage_rows$schedule_id[[i]] &
            fixture$cumulative_summary$evidence_threshold ==
                coverage_rows$evidence_threshold[[i]],
        , drop = FALSE]
    expect_equal(nrow(curve), final_row$n_looks[[1]],
                 info = paste("selected sequential point cumulative curve exists",
                              coverage_rows$bf_prior_id[[i]],
                              coverage_rows$schedule_id[[i]],
                              coverage_rows$evidence_threshold[[i]]))
}
