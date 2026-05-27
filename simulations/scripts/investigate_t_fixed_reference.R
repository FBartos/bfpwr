source(file.path("simulations", "scripts", "common.R"))
bfpwr_sim_source_package()

corpus_root <- bfpwr_sim_require_corpus_root(
    "simulations/corpus/v1",
    asset_set = "fixture-tests",
    purpose = "t fixed-reference investigation")
fixture <- bfpwr_sim_read_fixture_summary(
    corpus_root = corpus_root,
    fixture_set_id = "t-tbf01-fixed-core-v1",
    family = "t",
    mode = "fixed")
designs <- bfpwr_sim_design_case_set("production")
bf_priors <- bfpwr_sim_bf_prior_case_set("production")

find_row <- function(design_case_id, bf_prior_id, tail,
                     evidence_threshold, n) {
    row <- fixture$tail_summary[
        fixture$tail_summary$design_case_id == design_case_id &
            fixture$tail_summary$bf_prior_id == bf_prior_id &
            fixture$tail_summary$tail == tail &
            fixture$tail_summary$evidence_threshold == evidence_threshold &
            fixture$tail_summary$n == n,
        , drop = FALSE]
    if (nrow(row) != 1L) {
        stop("expected exactly one target row, found ", nrow(row))
    }
    row
}

exact_interval_probability <- function(roots, k, lower.tail, n1, n2, type,
                                       alternative, plocation, pscale, pdf,
                                       df, neff, dpm, dpsd) {
    roots <- roots[is.finite(roots)]
    breaks <- c(-Inf, sort(unique(roots)), Inf)
    probability <- 0
    intervals <- vector("list", length(breaks) - 1L)
    for (i in seq_len(length(breaks) - 1L)) {
        lo <- breaks[[i]]
        hi <- breaks[[i + 1L]]
        mid <- if (is.infinite(lo)) {
            hi - 8
        } else if (is.infinite(hi)) {
            lo + 8
        } else {
            mean(c(lo, hi))
        }
        bf_le_k <- tbf01(
            t = mid, n1 = n1, n2 = n2,
            plocation = plocation, pscale = pscale, pdf = pdf,
            type = type, alternative = alternative, log = TRUE) <= log(k)
        include <- if (lower.tail) bf_le_k else !bf_le_k
        interval_prob_given_delta <- function(delta) {
            ncp <- sqrt(neff) * delta
            stats::pt(hi, df = df, ncp = ncp) -
                stats::pt(lo, df = df, ncp = ncp)
        }
        p <- if (dpsd == 0) {
            interval_prob_given_delta(dpm)
        } else {
            f <- function(delta) {
                interval_prob_given_delta(delta) *
                    stats::dnorm(delta, mean = dpm, sd = dpsd)
            }
            stats::integrate(f, lower = dpm - 8 * dpsd,
                             upper = dpm + 8 * dpsd,
                             subdivisions = 200L,
                             rel.tol = 1e-8)$value
        }
        intervals[[i]] <- data.frame(
            lower = lo, upper = hi, midpoint = mid,
            bf01_le_k = bf_le_k, probability = p,
            included = include,
            stringsAsFactors = FALSE)
        if (include) probability <- probability + p
    }
    list(probability = probability, intervals = do.call(rbind, intervals))
}

exact_point_or_mixture_reference <- function(row) {
    design <- bfpwr_sim_find_design_case(row$design_case_id[[1]], designs)
    bf_prior <- bfpwr_sim_find_bf_prior_case(row$bf_prior_id[[1]], bf_priors)
    prior <- bf_prior$analysis_prior
    dp <- bfpwr_sim_design_prior_mean_sd(design)

    n1 <- row$n[[1]]
    n2 <- if (identical(design$generation$type, "two.sample")) {
        pmax(2, round(n1 * design$generation$n2_multiplier))
    } else {
        n1
    }
    pars <- .tbf01_pars(n1 = n1, n2 = n2, type = prior$type)
    roots <- tcrit(k = row$threshold[[1]], n1 = n1, n2 = n2,
                   plocation = prior$plocation - prior$null,
                   pscale = prior$pscale, pdf = prior$pdf,
                   type = prior$type, alternative = prior$alternative,
                   tol = 1e-9)
    exact <- exact_interval_probability(
        roots = roots,
        k = row$threshold[[1]],
        lower.tail = row$tail[[1]] == "H1",
        n1 = n1,
        n2 = n2,
        type = prior$type,
        alternative = prior$alternative,
        plocation = prior$plocation - prior$null,
        pscale = prior$pscale,
        pdf = prior$pdf,
        df = pars$df,
        neff = pars$neff,
        dpm = dp$dpm - prior$null,
        dpsd = dp$dpsd)
    approx <- ptbf01(k = row$threshold[[1]], n = n1, n1 = n1, n2 = n2,
                     null = prior$null,
                     plocation = prior$plocation,
                     pscale = prior$pscale,
                     pdf = prior$pdf,
                     dpm = dp$dpm,
                     dpsd = dp$dpsd,
                     type = prior$type,
                     alternative = prior$alternative,
                     lower.tail = row$tail[[1]] == "H1")
    data.frame(
        design_case_id = row$design_case_id,
        bf_prior_id = row$bf_prior_id,
        type = prior$type,
        alternative = prior$alternative,
        dpm = dp$dpm,
        dpsd = dp$dpsd,
        tail = row$tail,
        evidence_threshold = row$evidence_threshold,
        n = n1,
        df = pars$df,
        neff = pars$neff,
        fixture_prob = row$prob,
        fixture_mcse = row$mcse,
        ptbf01 = approx,
        exact_t = exact$probability,
        fixture_minus_exact_z =
            (row$prob - exact$probability) / row$mcse,
        ptbf01_minus_exact = approx - exact$probability,
        stringsAsFactors = FALSE)
}

case_specs <- list(
    paired_two_sided = list(
        design = "t-paired-dpoint-0p5-short",
        prior = paste0(
            "t-tbf01-cauchy-null0-loc0-scale0p35-df1-two-sided-on-",
            "t-paired-dpoint-0p5-short"),
        tail = "H1", threshold = 30,
        n = c(10, 20, 30, 50, 100)),
    one_sample_two_sided = list(
        design = "t-one-dpoint-0p5-short",
        prior = paste0(
            "t-tbf01-cauchy-null0-loc0-scale0p35-df1-two-sided-on-",
            "t-one-dpoint-0p5-short"),
        tail = "H1", threshold = 30,
        n = c(10, 20, 30, 50, 100)),
    two_sample_two_sided = list(
        design = "t-two-dpoint-0p5-short",
        prior = paste0(
            "t-tbf01-cauchy-null0-loc0-scale0p35-df1-two-sided-on-",
            "t-two-dpoint-0p5-short"),
        tail = "H1", threshold = 30,
        n = c(10, 20, 30, 50, 100)),
    paired_greater = list(
        design = "t-paired-dpoint-0p5-short",
        prior = paste0(
            "t-tbf01-cauchy-null0-loc0-scale0p70710678-df1-greater-on-",
            "t-paired-dpoint-0p5-short"),
        tail = "H1", threshold = 30,
        n = c(10, 20, 30, 50, 100)),
    two_sample_uncertain = list(
        design = "t-two-dnorm-0p5-s0p5-short",
        prior = paste0(
            "t-tbf01-cauchy-null0-loc0-scale0p35-df1-two-sided-on-",
            "t-two-dnorm-0p5-s0p5-short"),
        tail = "H1", threshold = 30,
        n = c(10, 20, 30, 50, 100))
)

rows <- list()
idx <- 0L
for (case_name in names(case_specs)) {
    spec <- case_specs[[case_name]]
    for (n in spec$n) {
        idx <- idx + 1L
        row <- find_row(spec$design, spec$prior, spec$tail,
                        spec$threshold, n)
        rows[[idx]] <- cbind(case = case_name,
                             exact_point_or_mixture_reference(row),
                             stringsAsFactors = FALSE)
    }
}

out <- do.call(rbind, rows)
out$abs_ptbf01_exact <- abs(out$ptbf01_minus_exact)
out$abs_fixture_exact_z <- abs(out$fixture_minus_exact_z)

print(out[, c("case", "type", "dpm", "dpsd", "n", "df", "neff",
              "fixture_prob", "fixture_mcse", "ptbf01", "exact_t",
              "ptbf01_minus_exact", "fixture_minus_exact_z")],
      row.names = FALSE, digits = 5)

cat("\nWorst ptbf01 vs exact rows:\n")
ord <- order(-out$abs_ptbf01_exact)
print(out[ord[seq_len(min(8L, nrow(out)))],
          c("case", "type", "dpm", "dpsd", "n", "ptbf01", "exact_t",
            "ptbf01_minus_exact", "fixture_prob", "fixture_minus_exact_z")],
      row.names = FALSE, digits = 5)

cat("\nWorst fixture vs exact rows:\n")
ord <- order(-out$abs_fixture_exact_z)
print(out[ord[seq_len(min(8L, nrow(out)))],
          c("case", "type", "dpm", "dpsd", "n", "fixture_prob",
            "fixture_mcse", "exact_t", "fixture_minus_exact_z")],
      row.names = FALSE, digits = 5)
