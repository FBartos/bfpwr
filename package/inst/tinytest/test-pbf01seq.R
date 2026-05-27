library(tinytest)
library(bfpwr)

## Tests pbf01seq identities and invariants for normal, normal-moment, and
## directional z designs. Fixed-sample sources are paper/bfssd.Rnw 354-369,
## 449-606, and 1672-1709; these sequential cases are package regressions rather
## than literal manuscript examples. Direct BFGSD Low-PV checks live in
## test-paper-sequential-regions.R.

## Helper for generic sequential invariants: each look has cumulative H1/H0/Inc,
## probabilities are monotone and exhaustive, and EN stays inside the schedule.
expect_pbf01seq_long_run_invariants <- function(res, n_looks, label) {
    expect_equal(length(res$cumpH1), n_looks,
                 info = paste(label, "has one cumulative H1 value per look"))
    expect_equal(length(res$cumpH0), n_looks,
                 info = paste(label, "has one cumulative H0 value per look"))
    expect_equal(length(res$cumpInc), n_looks,
                 info = paste(label, "has one cumulative inconclusive value per look"))
    expect_true(all(diff(res$cumpH1) >= -1e-12),
                info = paste(label, "cumulative H1 is nondecreasing"))
    expect_true(all(diff(res$cumpH0) >= -1e-12),
                info = paste(label, "cumulative H0 is nondecreasing"))
    expect_true(all(diff(res$cumpInc) <= 1e-12),
                info = paste(label, "cumulative inconclusive is nonincreasing"))
    expect_true(max(abs(res$cumpH1 + res$cumpH0 + res$cumpInc - 1)) <
                    1e-10,
                info = paste(label, "cumulative probabilities sum to one"))
    expect_true(is.finite(res$EN) && res$EN >= min(res$n) &&
                    res$EN <= max(res$n),
                info = paste(label, "expected sample size is within schedule"))
}

dirbf01_threshold <- function(k, se, null, pm, psd) {
    f <- function(x) {
        dirbf01(estimate = x, se = se, null = null, pm = pm, psd = psd,
                log = TRUE) - log(k)
    }
    center <- null
    width <- max(se, psd, 1)
    lower <- center - width
    upper <- center + width
    for (i in seq_len(60)) {
        f_lower <- f(lower)
        f_upper <- f(upper)
        if (is.finite(f_lower) && is.finite(f_upper) &&
            f_lower * f_upper <= 0) {
            return(stats::uniroot(f, lower = lower, upper = upper,
                                  tol = 1e-12)$root)
        }
        width <- width * 2
        lower <- center - width
        upper <- center + width
    }
    stop("could not bracket directional BF threshold")
}

pdirbf01_reference <- function(k, n, usd, null, pm, psd, dpm, dpsd,
                               lower.tail = TRUE) {
    se <- usd / sqrt(n)
    boundary <- dirbf01_threshold(k = k, se = se, null = null, pm = pm,
                                  psd = psd)
    design_sd <- sqrt(se^2 + dpsd^2)
    if (lower.tail) {
        stats::pnorm(boundary, mean = dpm, sd = design_sd,
                     lower.tail = FALSE)
    } else {
        stats::pnorm(boundary, mean = dpm, sd = design_sd,
                     lower.tail = TRUE)
    }
}

## One-look normal pbf01seq should reduce exactly to fixed-sample pbf01 from the
## BF01 power formulas (paper/bfssd.Rnw 449-606).
one_look_n <- 100
one_look_se <- 1 / sqrt(one_look_n)
one_look <- pbf01seq(k1 = 1 / 10,
                     k0 = 10,
                     se = one_look_se,
                     n = one_look_n,
                     pm = 0.2,
                     psd = 0,
                     dpm = 0.2,
                     dpsd = 0,
                     type = "normal",
                     strict = TRUE)
expect_equal(one_look$cumpH1,
             pbf01(k = 1 / 10, n = one_look_n, usd = 1, null = 0,
                   pm = 0.2, psd = 0, dpm = 0.2, dpsd = 0,
                   lower.tail = TRUE),
             tolerance = 1e-12,
             info = "one-look normal pbf01seq H1 agrees with pbf01")
expect_equal(one_look$cumpH0,
             pbf01(k = 10, n = one_look_n, usd = 1, null = 0,
                   pm = 0.2, psd = 0, dpm = 0.2, dpsd = 0,
                   lower.tail = FALSE),
             tolerance = 1e-12,
             info = "one-look normal pbf01seq H0 agrees with pbf01")

## One-look moment pbf01seq should reduce exactly to pnmbf01 from the
## normal-moment power formula (paper/bfssd.Rnw 1693-1709).
moment_one_look <- pbf01seq(k1 = 1 / 10,
                            k0 = 10,
                            se = one_look_se,
                            n = one_look_n,
                            psd = 1 / sqrt(2),
                            dpm = 0.2,
                            dpsd = 0.1,
                            type = "moment",
                            strict = TRUE)
expect_equal(moment_one_look$cumpH1,
             pnmbf01(k = 1 / 10, n = one_look_n, usd = 1, null = 0,
                     psd = 1 / sqrt(2), dpm = 0.2, dpsd = 0.1,
                     lower.tail = TRUE),
             tolerance = 1e-12,
             info = "one-look moment pbf01seq H1 agrees with pnmbf01")
expect_equal(moment_one_look$cumpH0,
             pnmbf01(k = 10, n = one_look_n, usd = 1, null = 0,
                     psd = 1 / sqrt(2), dpm = 0.2, dpsd = 0.1,
                     lower.tail = FALSE),
             tolerance = 1e-12,
             info = "one-look moment pbf01seq H0 agrees with pnmbf01")

## One-look directional pbf01seq has no manuscript source; it is checked against
## an independent directional threshold calculation in this file.
directional_one_look <- pbf01seq(k1 = 1 / 10,
                                 k0 = 10,
                                 se = one_look_se,
                                 n = one_look_n,
                                 pm = 0,
                                 psd = 1 / sqrt(2),
                                 dpm = 0.2,
                                 dpsd = 0.1,
                                 type = "directional",
                                 strict = TRUE)
expect_equal(directional_one_look$cumpH1,
             pdirbf01_reference(k = 1 / 10, n = one_look_n, usd = 1,
                                null = 0, pm = 0, psd = 1 / sqrt(2),
                                dpm = 0.2, dpsd = 0.1,
                                lower.tail = TRUE),
             tolerance = 1e-10,
             info = "one-look directional pbf01seq H1 agrees with direct directional reference")
expect_equal(directional_one_look$cumpH0,
             pdirbf01_reference(k = 10, n = one_look_n, usd = 1,
                                null = 0, pm = 0, psd = 1 / sqrt(2),
                                dpm = 0.2, dpsd = 0.1,
                                lower.tail = FALSE),
             tolerance = 1e-10,
             info = "one-look directional pbf01seq H0 agrees with direct directional reference")

## Nonzero-null one-look checks verify package recentering against fixed-sample
## pbf01/pnmbf01 and the independent directional reference; not manuscript rows.
shifted_null <- 0.2
shifted_n <- 80
shifted_se <- 1 / sqrt(shifted_n)
shifted_normal <- pbf01seq(k1 = 1 / 10,
                           k0 = 10,
                           se = shifted_se,
                           n = shifted_n,
                           pm = 0.5 - shifted_null,
                           psd = 0,
                           dpm = 0.35 - shifted_null,
                           dpsd = 0.05,
                           type = "normal",
                           strict = TRUE)
expect_equal(shifted_normal$cumpH1,
             pbf01(k = 1 / 10, n = shifted_n, usd = 1,
                   null = shifted_null, pm = 0.5, psd = 0,
                   dpm = 0.35, dpsd = 0.05, lower.tail = TRUE),
             tolerance = 1e-12,
             info = "nonzero-null normal pbf01seq recentering agrees with pbf01 H1")
expect_equal(shifted_normal$cumpH0,
             pbf01(k = 10, n = shifted_n, usd = 1,
                   null = shifted_null, pm = 0.5, psd = 0,
                   dpm = 0.35, dpsd = 0.05, lower.tail = FALSE),
             tolerance = 1e-12,
             info = "nonzero-null normal pbf01seq recentering agrees with pbf01 H0")

shifted_moment <- pbf01seq(k1 = 1 / 10,
                           k0 = 10,
                           se = shifted_se,
                           n = shifted_n,
                           psd = 0.5,
                           dpm = 0.35 - shifted_null,
                           dpsd = 0.05,
                           type = "moment",
                           strict = TRUE)
expect_equal(shifted_moment$cumpH1,
             pnmbf01(k = 1 / 10, n = shifted_n, usd = 1,
                     null = shifted_null, psd = 0.5,
                     dpm = 0.35, dpsd = 0.05, lower.tail = TRUE),
             tolerance = 1e-12,
             info = "nonzero-null moment pbf01seq recentering agrees with pnmbf01 H1")
expect_equal(shifted_moment$cumpH0,
             pnmbf01(k = 10, n = shifted_n, usd = 1,
                     null = shifted_null, psd = 0.5,
                     dpm = 0.35, dpsd = 0.05, lower.tail = FALSE),
             tolerance = 1e-12,
             info = "nonzero-null moment pbf01seq recentering agrees with pnmbf01 H0")

shifted_directional <- pbf01seq(k1 = 1 / 10,
                                k0 = 10,
                                se = shifted_se,
                                n = shifted_n,
                                pm = 0.5 - shifted_null,
                                psd = 0.5,
                                dpm = 0.35 - shifted_null,
                                dpsd = 0.05,
                                type = "directional",
                                strict = TRUE)
expect_equal(shifted_directional$cumpH1,
             pdirbf01_reference(k = 1 / 10, n = shifted_n, usd = 1,
                                null = shifted_null, pm = 0.5, psd = 0.5,
                                dpm = 0.35, dpsd = 0.05,
                                lower.tail = TRUE),
             tolerance = 1e-10,
             info = "nonzero-null directional pbf01seq recentering agrees with direct directional H1")
expect_equal(shifted_directional$cumpH0,
             pdirbf01_reference(k = 10, n = shifted_n, usd = 1,
                                null = shifted_null, pm = 0.5, psd = 0.5,
                                dpm = 0.35, dpsd = 0.05,
                                lower.tail = FALSE),
             tolerance = 1e-10,
             info = "nonzero-null directional pbf01seq recentering agrees with direct directional H0")

## Long-schedule checks are package-level invariants for cumulative probabilities
## and expected sample sizes over many looks.
long_n <- seq(10, 1000, by = 10)
normal_long <- pbf01seq(k1 = 1 / 10,
                        k0 = 10,
                        se = 1 / sqrt(long_n),
                        n = long_n,
                        pm = 0.2,
                        psd = 0,
                        dpm = 0.2,
                        dpsd = 0,
                        type = "normal",
                        strict = TRUE)
expect_pbf01seq_long_run_invariants(normal_long, length(long_n),
                                    "100-look strict normal pbf01seq")

directional_long <- pbf01seq(k1 = 1 / 10,
                             k0 = 10,
                             se = 1 / sqrt(long_n),
                             n = long_n,
                             pm = 0,
                             psd = 1 / sqrt(2),
                             dpm = 0.2,
                             dpsd = 0,
                             type = "directional",
                             strict = TRUE)
expect_pbf01seq_long_run_invariants(directional_long, length(long_n),
                                    "100-look strict directional pbf01seq")

moment_n <- seq(10, 50, by = 10)
moment_long <- pbf01seq(k1 = 1 / 10,
                        k0 = 10,
                        se = 1 / sqrt(moment_n),
                        n = moment_n,
                        psd = 1 / sqrt(2),
                        dpm = 0.2,
                        dpsd = 0,
                        type = "moment",
                        strict = TRUE)
expect_pbf01seq_long_run_invariants(moment_long, length(moment_n),
                                    "5-look strict moment pbf01seq")

## ## check with simulation that correct probabilities calculated
## simbenchmark <- function(nsim, k1, k0, usd, n, pm, psd, dpm, dpsd, type) {
##     ## simulate stage-wise BFs
##     smd <- rnorm(n = nsim, mean = dpm, sd = dpsd)
##     bfmat <- sapply(X = smd, FUN = function(smdi) {
##         y1 <- rnorm(n = max(n), mean = 0, sd = usd)
##         y2 <- rnorm(n = max(n), mean = smdi*usd, sd = usd)
##         est <- sapply(seq_along(n), FUN = function(i) {
##             (mean(y2[1:n[i]]) - mean(y1[1:n[i]]))/usd
##         })
##         se <- sapply(seq_along(n), FUN = function(i) {
##             sqrt(2/n[i])
##         })
##         bf <- sapply(seq_along(n), FUN = function(i) {
##             if (type == "normal") {
##                 bf01(estimate = est[i], se = se[i], null = 0, pm = pm,
##                      psd = psd)
##             } else if (type == "directional") {
##                 dirbf01(estimate = est[i], se = se[i], null = 0, pm = pm,
##                         psd = psd)
##             } else {
##                 nmbf01(estimate = est[i], se = se[i], null = 0, psd = psd)
##             }
##         })
##         return(bf)
##     })

##     ## estimate probabilities
##     stop <- apply(X = bfmat, MARGIN = 2, FUN = function(x) {
##         result <- "inconclusive"
##         for (xi in x) {
##             if (xi >= k0) {
##                 result <- "H0"
##                 break
##             }
##             if (xi <= k1) {
##                 result <- "H1"
##                 break
##             }
##         }
##         return(result)
##     })
##     pH1sim <- mean(stop == "H1")
##     pH0sim <- mean(stop == "H0")
##     pIncsim <- mean(stop == "inconclusive")

##     ## compute probabilities numerically
##     se <- sqrt(2/n)
##     res <- bfpwr::pbf01seq(k1 = k1, k0 = k0, se = se, n = n, pm = pm,
##                            psd = psd, dpm = dpm, dpsd = dpsd, type = type,
##                            strict = TRUE)

##     ## put everything together
##     out <- data.frame(method = c("simulation", "numerical"),
##                       pH1 = c(pH1sim, tail(res$cumpH1, n = 1)),
##                       pH0 = c(pH0sim, tail(res$cumpH0, n = 1)),
##                       pInc = c(pIncsim, tail(res$cumpInc, n = 1)))
##     return(out)
## }

## simgrid <- expand.grid(k1 = c(1/10, 0.3),
##                        k0 = c(3, 10),
##                        usd = c(sqrt(2)),
##                        pm = c(-0.25, 0.5),
##                        dpm = c(0, 0.5),
##                        psd = c(0.1, 1),
##                        dpsd = c(0, 0.1),
##                        type = c("normal", "moment", "directional"),
##                        stringsAsFactors = FALSE)
## n <- seq(25, 250, 25)
## set.seed(42)
## nsim <- 10000
## simres <- do.call("rbind", lapply(X = seq(1, nrow(simgrid)), FUN = function(i) {
##     res <- simbenchmark(nsim = nsim, k1 = simgrid$k1[i], k0 = simgrid$k0[i],
##                         usd = simgrid$usd[i], n = n, pm = simgrid$pm[i],
##                         psd = simgrid$psd[i], dpm = simgrid$dpm[i],
##                         dpsd = simgrid$dpsd[i], type = simgrid$type[i])
##     condition <- simgrid[i,]
##     data.frame(i = i, condition, res)
## }))

## for (i in seq(1, nrow(simgrid))) {
##     pH1sim <- simres[simres$i == i,]$pH1[1]
##     pH1num <- simres[simres$i == i,]$pH1[2]
##     print(expect_equal(pH1sim, pH1num, info = i, tolerance = 0.01))
## }

## library(ggplot2)
## simplt <- ggplot(data = simres, aes(x = method, y = pH0, color = method,
##                                     group = interaction(pm, dpm, psd, type, k1,
##                                                         k0, dpsd))) +
##     facet_grid(pm + dpm + psd + type ~ k1 + k0 + dpsd,
##                labeller = label_both) +
##     geom_line(color = 1, alpha = 0.8) +
##     geom_point(size = 1) +
##     theme_bw()
## ggsave(plot = simplt, filename = "simbenchmark.pdf", width = 20, height = 20)
