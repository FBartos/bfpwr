library(tinytest)
library(bfpwr)

## A two-look design one observation apart has Bernoulli stopping variance,
## regardless of the absolute sample sizes or the rarity of early stopping.
for (p in c(0, 1e-9, 0.25, 1)) {
    moments <- bfpwr:::.bfseq_sample_size_moments(
        pH1 = c(p, 0), pH0 = c(0, 0), n = c(1e5, 1e5 + 1)
    )
    expect_equal(moments$EN, 1e5 + 1 - p)
    expect_equal(moments$VarN, p*(1 - p), tolerance = 1e-14,
                 info = "sample-size variance avoids cancellation at large N")
}

logs <- log(c(0.2, 0.3))
expect_equal(
    bfpwr:::.bfpwr_logspace_sum(c(logs, -Inf)),
    log(0.5),
    tolerance = 1e-12,
    info = "log-space summation should ignore -Inf zero-mass terms"
)

expect_equal(
    bfpwr:::.bfpwr_logspace_sum(c(-Inf, -Inf)),
    -Inf,
    info = "log-space summation should return -Inf when all terms are zero mass"
)

expect_error(
    bfpwr:::.bfpwr_logspace_sum(c(logs, NA_real_)),
    "invalid non-finite",
    info = "log-space summation should reject NA values"
)

expect_error(
    bfpwr:::.bfpwr_logspace_sum(c(logs, NaN)),
    "invalid non-finite",
    info = "log-space summation should reject NaN values"
)

expect_error(
    bfpwr:::.bfpwr_logspace_sum(c(logs, Inf)),
    "invalid non-finite",
    info = "log-space summation should reject positive infinity"
)

expect_equal(
    bfpwr:::.bfpwr_lpnorm_interval(
        lower = c(-Inf, -1, 1),
        upper = c(0, 1, Inf),
        mean = c(0, 0.5, 0),
        sd = c(1, 2, 1)
    ),
    log(c(
        stats::pnorm(0),
        stats::pnorm(1, mean = 0.5, sd = 2) -
            stats::pnorm(-1, mean = 0.5, sd = 2),
        stats::pnorm(1, lower.tail = FALSE)
    )),
    tolerance = 1e-12,
    info = "normal interval helper should handle vector inputs"
)



## Interior beta interval masses are subtracted on the side of the interval.
## Here the interval lies in the far upper tail, where the lower-tail
## difference cancels 1 - 1 and used to return -36.04365 (reference by
## adaptive integration of the density: -36.240947703761).
expect_equal(
    bfpwr:::.bfpwr_lpbeta_interval(0.3630126, 0.67826436, 39.92865881,
                                   236.48482111),
    -36.240947703761, tolerance = 1e-10,
    info = "interval masses in the upper tail should not use the cancelling difference"
)
expect_equal(
    bfpwr:::.bfpwr_lpbeta_interval(1 - 0.67826436, 1 - 0.3630126,
                                   236.48482111, 39.92865881),
    -36.240947703761, tolerance = 1e-10,
    info = "mirrored interval masses in the lower tail should agree"
)
expect_equal(
    sum(exp(bfpwr:::.binbf01_pred_logpmf(
        x = 0:63, n = 63, dp = NA, da = 14.92865881, db = 198.48482111,
        dl = 0.3630126, du = 0.67826436
    ))),
    1, tolerance = 1e-10,
    info = "truncated beta predictive masses should sum to one"
)

for (shapes in list(c(0.5, 0.5), c(3, 4), c(40, 240), c(240, 40), c(2, 300))) {
    for (limits in list(c(0.2, 0.7), c(0.001, 0.01), c(0.4, 0.9), c(0.99, 0.999))) {
        direct <- log(diff(stats::pbeta(limits, shapes[1], shapes[2])))
        mirrored <- bfpwr:::.bfpwr_lpbeta_interval(
            lower = 1 - limits[2], upper = 1 - limits[1],
            shape1 = shapes[2], shape2 = shapes[1]
        )
        value <- bfpwr:::.bfpwr_lpbeta_interval(
            lower = limits[1], upper = limits[2],
            shape1 = shapes[1], shape2 = shapes[2]
        )
        expect_equal(value, mirrored, tolerance = 1e-9,
                     info = "interval masses should be symmetric under reflection")
        if (is.finite(direct) && direct > -4) {
            expect_equal(value, direct, tolerance = 1e-9,
                         info = "interval masses should match the direct difference where it is accurate")
        }
    }
}
