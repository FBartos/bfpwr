library(tinytest)
library(bfpwr)

## Under a normal prior, the predictive t statistic is a scaled noncentral t.
n <- 30
t <- 2.5
m <- 0.2
s <- 0.5
scale <- sqrt(1 + n*s^2)
expected <- dt(t, df = n - 1)/(dt(t/scale, df = n - 1,
                                   ncp = sqrt(n)*m/scale)/scale)
expect_equal(tbf01(t = t, n = n, plocation = m, pscale = s, pdf = Inf,
                   type = "one.sample"), expected, tolerance = 1e-8)

for (alternative in c("less", "greater")) {
    bf <- tbf01(t = t, n = n, plocation = m, pscale = s, pdf = c(Inf, 1e8),
                type = "one.sample", alternative = alternative)
    expect_equal(bf[1], bf[2], tolerance = 1e-6, info = alternative)
}
for (pdf in list(-Inf, -1, 0, NA, NA_real_, NaN, "Inf")) {
    expect_error(tbf01(t = t, n = n, pdf = pdf))
}

## Fixed and sequential designs should agree with their large-df limits.
power <- ptbf01(k = 1/6, n = 50, pdf = c(Inf, 1e8), dpm = 0.5, dpsd = 0)
expect_equal(power[1], power[2], tolerance = 1e-6)
expect_equal(powertbf01(k = 1/6, n = 50, pdf = Inf, dpm = 0.5, dpsd = 0)$power,
             power[2], tolerance = 1e-6)

nfixed <- ntbf01(k = 1/6, power = 0.5, pdf = c(Inf, 1e8), dpm = 0.5,
                 dpsd = 0, nrange = c(30, 100))
expect_true(is.finite(nfixed[1]))
expect_equal(nfixed[1], nfixed[2])

seq_normal <- ptbf01seq(k1 = 1/6, k0 = 2, n = c(30, 50), pdf = Inf,
                        dpm = 0.5, dpsd = 0)
seq_t <- ptbf01seq(k1 = 1/6, k0 = 2, n = c(30, 50), pdf = 1e8,
                   dpm = 0.5, dpsd = 0)
expect_equal(seq_normal$cumpH1, seq_t$cumpH1, tolerance = 1e-6)
expect_equal(seq_normal$cumpH0, seq_t$cumpH0, tolerance = 1e-6)

nseq <- ntbf01seq(k1 = 1/6, k0 = 2, power = 0.5, pdf = c(Inf, 1e8),
                  dpm = 0.5, dpsd = 0, nrange = c(30, 100))
expect_true(is.finite(nseq[1]))
expect_equal(nseq[1], nseq[2])

## Integrate on the prior's own scale to resolve a narrow, shifted peak.
for (pdf in c(30, 1e8, Inf)) {
    ref <- integrate(function(z) dt(2, 98, ncp = 5*(0.5 + 0.005*z))*dt(z, pdf),
                     -20, 20, rel.tol = 1e-10)$value
    expect_equal(tbf01(2, 50, plocation = 0.5, pscale = 0.005, pdf = pdf,
                       log = TRUE), log(dt(2, 98)/ref), tolerance = 1e-8)
}

## The normal prior gives independent references in extreme likelihood tails.
scale <- sqrt(1 + 5000*0.1^2)
predictive <- integrate(function(p) {
    denominator <- sqrt(qchisq(p, df = 19998)/19998)
    denominator*dnorm(denominator, mean = sqrt(5000)*0.5, sd = scale)
}, 0, 1, rel.tol = 1e-10, abs.tol = 0)$value
expected <- dt(1, 19998, log = TRUE) - log(predictive)
expect_equal(tbf01(1, 10000, plocation = 0.5, pscale = 0.1, pdf = Inf,
                   log = TRUE), expected, tolerance = 1e-8)
expect_equal(tbf01(1, 10000, plocation = 0.5, pscale = 0.1, pdf = 1e8,
                   log = TRUE), expected, tolerance = 1e-6)
## Independent integration over the chi-square denominator of the t statistic.
expect_equal(tbf01(-20, 7880, plocation = 0.35, pscale = 0.02, pdf = Inf,
                   alternative = "greater", log = TRUE), 156.0745, tolerance = 1e-6)

## At n = 2 per group the centered normal BF01 cannot fall below 2/3.
## The first look has no H1 stopping region, but the next look can stop.
for (pdf in c(30, Inf)) {
    expect_equal(ptbf01(1/6, 2, pdf = pdf), 0)
    expect_equal(ptbf01(1/6, 2, pdf = pdf, lower.tail = FALSE), 1)
    design <- suppressWarnings(ptbf01seq(1/6, k0 = 2, n = c(2, 50), pdf = pdf,
                                         dpm = 0.5, dpsd = 0))
    expect_equal(design$cumpH1[1], 0)
    expect_true(design$cumpH1[2] > 0)
    searched <- suppressWarnings(powertbf01seq(k1 = 1/6, k0 = 2, power = 0.5,
        pdf = pdf, dpm = 0.5, dpsd = 0, looks = 2))
    expect_true(searched$solver$reached)
}
one_sided <- suppressWarnings(ptbf01seq(1/6, k0 = 3, n = c(2, 20), pdf = Inf,
    dpm = 1, dpsd = 3, alternative = "greater"))
expect_equal(one_sided$cumpH1[1], 0)
expect_true(one_sided$cumpH1[2] > 0)

## A normal tail far beyond its mode concentrates at the truncation point.
## The Cauchy tail instead has limiting density a/(a + d)^2 for d > 0.
expect_equal(tbf01(2, 50, plocation = -0.5, pscale = 1e-5, pdf = Inf,
                   alternative = "greater"), 1, tolerance = 1e-8)
ref <- integrate(function(d) dt(2, 98, ncp = 5*d)*0.5/(0.5 + d)^2,
                 0, 10, rel.tol = 1e-10)$value
expect_equal(tbf01(2, 50, plocation = -0.5, pscale = 1e-8, pdf = 1,
                   alternative = "greater"), dt(2, 98)/ref, tolerance = 1e-8)

## When the prior and sampling df coincide, BF01 can approach zero slowly.
## A remote finite crossing must still reproduce its threshold.
remote <- suppressWarnings(ptbf01seq(1/100, k0 = 2, n = 2, pdf = 2))
expect_true(all(is.finite(remote$zk1)))
expect_equal(tbf01(as.numeric(remote$zk1), n = 2, pdf = 2), rep(1/100, 2),
             tolerance = 1e-6)
expect_true(is.nan(suppressWarnings(ptbf01(1/6, 50, pdf = Inf, maxiter = 1,
                                           check.conv = TRUE))))

## A shifted two-sided prior can cross a threshold in just one tail.
for (pdf in c(1e8, Inf)) {
    design <- ptbf01seq(1/6, k0 = 3, n1 = 53, n2 = 57, plocation = 0.35,
                        pscale = 0.102, pdf = pdf, dpm = 0.35, dpsd = 0)
    expect_equal(c(design$zk1[1], design$zk0[1]), c(-Inf, -Inf))
    expect_equal(tbf01(c(design$zk1[2], design$zk0[2]), n1 = 53, n2 = 57,
                       plocation = 0.35, pscale = 0.102, pdf = pdf),
                 c(1/6, 3), tolerance = 1e-6)
}
