library(tinytest)
library(bfpwr)

## A change of measurement units preserves the predictive probability.
for (alternative in c("greater", "less")) {
    direction <- if (alternative == "greater") 1 else -1
    for (dpsd in c(0, 0.3)) {
        reference <- pnmbf01(k = .1, n = 25, usd = 1, null = 2, psd = .7,
            dpm = 2 + direction*.4, dpsd = dpsd, alternative = alternative,
            lower.tail = c(TRUE, FALSE))
        if (dpsd > 0) expect_equal(reference[1], .34439877128421587)
        for (scale in c(1e-200, 1, 1e200)) {
            actual <- pnmbf01(k = .1, n = 25, usd = scale, null = 2*scale,
                psd = .7*scale, dpm = (2 + direction*.4)*scale,
                dpsd = dpsd*scale, alternative = alternative,
                lower.tail = c(TRUE, FALSE))
            expect_equal(actual, reference, tolerance = 1e-12)
            expect_equal(sum(actual), 1, tolerance = 1e-12)
        }
    }
}
