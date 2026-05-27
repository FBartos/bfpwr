#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(bfpwr))

fmt <- function(x, digits = 10) {
    paste(format(signif(x, digits = digits), scientific = FALSE, trim = TRUE),
          collapse = ", ")
}

status <- function(ok) {
    if (is.na(ok)) return("INFO")
    if (isTRUE(ok)) "PASS" else "FAIL"
}

row <- function(id, source, expected, value, tolerance, ok, notes = "") {
    data.frame(
        id = id,
        source = source,
        expected = expected,
        package_value = value,
        tolerance = tolerance,
        status = status(ok),
        notes = notes,
        stringsAsFactors = FALSE
    )
}

rows <- list()

## Fixed-sample formulas exercised by package/inst/tinytest/test-paper-fixed-formulas.R
bf_expected <- exp(stats::dnorm(0.23, mean = 0.1, sd = 0.14, log = TRUE) -
                   stats::dnorm(0.23, mean = 0.4,
                                sd = sqrt(0.14^2 + 0.25^2), log = TRUE))
bf_value <- bf01(estimate = 0.23, se = 0.14, null = 0.1, pm = 0.4,
                 psd = 0.25)
rows[[length(rows) + 1]] <- row(
    "bf01 density ratio",
    "paper/bfssd.Rnw eq. BF01; test-paper-fixed-formulas.R",
    fmt(bf_expected),
    fmt(bf_value),
    "1e-12",
    isTRUE(all.equal(bf_value, bf_expected, tolerance = 1e-12))
)

nmbf_expected <- {
    estimate <- 0.25
    se <- 0.12
    psd <- 0.5/sqrt(2)
    marginal_sd <- sqrt(se^2 + psd^2)
    post_var <- 1/(1/se^2 + 1/psd^2)
    post_mean_diff <- psd^2/(se^2 + psd^2) * estimate
    moment_factor <- (post_var + post_mean_diff^2)/psd^2
    exp(stats::dnorm(estimate, mean = 0, sd = se, log = TRUE) -
        (stats::dnorm(estimate, mean = 0, sd = marginal_sd, log = TRUE) +
         log(moment_factor)))
}
nmbf_value <- nmbf01(estimate = 0.25, se = 0.12, null = 0,
                     psd = 0.5/sqrt(2))
rows[[length(rows) + 1]] <- row(
    "nmbf01 normal-moment marginal likelihood",
    "paper/bfssd.Rnw eqs. nlBF/pnlBF; test-paper-fixed-formulas.R",
    fmt(nmbf_expected),
    fmt(nmbf_value),
    "1e-12",
    isTRUE(all.equal(nmbf_value, nmbf_expected, tolerance = 1e-12))
)

bin_expected <- exp(17*log(0.5) + (25 - 17)*log1p(-0.5) +
                    lbeta(2, 3) - lbeta(2 + 17, 3 + 25 - 17))
bin_value <- binbf01(x = 17, n = 25, p0 = 0.5, type = "point", a = 2, b = 3)
rows[[length(rows) + 1]] <- row(
    "binbf01 beta-binomial point null",
    "package binomial formula; test-paper-fixed-formulas.R",
    fmt(bin_expected),
    fmt(bin_value),
    "1e-12",
    isTRUE(all.equal(bin_value, bin_expected, tolerance = 1e-12)),
    "Binomial functions are package checks; they are not part of paper/bfssd.Rnw."
)

rows[[length(rows) + 1]] <- row(
    "fixed-sample power/size wrapper corpus",
    "test-paper-fixed-formulas.R",
    "independent closed-form, direct integration, and exact enumeration references",
    "27 tinytest expectations",
    "1e-12 to 1e-4",
    NA,
    paste(c("pbf01", "nbf01", "powerbf01", "pnmbf01", "nnmbf01",
            "powernmbf01", "tbf01", "ptbf01", "ntbf01", "powertbf01",
            "pbinbf01", "nbinbf01", "powerbinbf01"),
          collapse = ", ")
)

## Local manuscript examples from paper/bfssd.Rnw
pm <- -6
pow <- 0.8
sd <- 15
alpha <- 0.05
n_freq <- 2*sd^2*(qnorm(p = 1 - alpha/2) + qnorm(p = pow))^2/pm^2
est <- -1.74
ci <- c(-7.17, 3.69)
se <- (ci[2] - ci[1])/(2*qnorm(p = 0.975))
p <- 2*pnorm(q = abs(est/se), lower.tail = FALSE)
lr01 <- bf01(estimate = est, se = se, null = 0, pm = pm, psd = 0)
rows[[length(rows) + 1]] <- row(
    "mirtazapine frequentist n",
    "paper/bfssd.Rnw lines 1060-1099",
    "ceiling(2*sd^2*(z_(1-alpha/2)+z_power)^2/pm^2)",
    paste0("n=", fmt(n_freq), "; ceiling=", ceiling(n_freq)),
    "rounding only",
    ceiling(n_freq) == 99,
    paste0("p=", fmt(p), "; se=", fmt(se))
)
rows[[length(rows) + 1]] <- row(
    "mirtazapine BF01",
    "paper/bfssd.Rnw lines 1109-1116",
    "BF01 from normal point-prior likelihood ratio",
    fmt(lr01),
    "rounding only",
    round(lr01, 1) == 2.7,
    "Rnw Sexpr round(lr01, 1) prints 2.7."
)

null <- 0
usd <- sqrt(2)*sd
psd <- 0
k <- 1/10
power <- 0.8
dpm <- pm
dpsd1 <- 0
dpsd2 <- 2
nnum <- nbf01(k = k, power = power, usd = usd, null = null, pm = pm,
              psd = psd, dpm = dpm, dpsd = dpsd1)
zb <- qnorm(p = power)
nanalyt <- ceiling((zb + sqrt(zb^2 - log(k^2)*(pm + null - 2*dpm)/
                              (null - pm)))^2/
                   (pm + null - 2*dpm)^2*usd^2)
a <- ((null + pm)/2 - dpm)^2 - zb^2*dpsd2^2
b <- usd^2*((null + pm - 2*dpm)*log(k)/(null - pm) - zb^2)
c <- (usd^2*log(k)/(null - pm))^2
nanalyt2 <- ceiling((-b + sqrt(b^2 - 4*a*c))/(2*a))
rows[[length(rows) + 1]] <- row(
    "mirtazapine BF sample size",
    "paper/bfssd.Rnw lines 1135-1143 and 1218-1224",
    paste0("nnum=", fmt(nnum), "; nanalyt=", nanalyt,
           "; nanalyt2=", nanalyt2),
    paste0("nnum=", fmt(nnum), "; ceiling=", ceiling(nnum),
           "; nanalyt=", nanalyt, "; nanalyt2=", nanalyt2),
    "1e-8 / integer equality",
    isTRUE(all.equal(ceiling(nnum), nanalyt, tolerance = 1e-8))
)

null <- 0
sd <- sqrt(2)
pm <- 0
psd <- 1/sqrt(2)
k <- 1/6
power <- 0.95
dpm <- 0.5
dpsd <- 0
dpsd2 <- 0.1
n_normal <- nbf01(k = k, power = power, usd = sd, null = null, pm = pm,
                  psd = psd, dpm = dpm, dpsd = dpsd)
n_normal2 <- nbf01(k = k, power = power, usd = sd, null = null, pm = pm,
                   psd = psd, dpm = dpm, dpsd = dpsd2)
n_normal_h0 <- nbf01(k = 1/k, power = power, usd = sd, null = null, pm = pm,
                     psd = psd, dpm = null, dpsd = 0, lower.tail = FALSE)
rows[[length(rows) + 1]] <- row(
    "Schoenbrodt normal-prior sample sizes",
    "paper/bfssd.Rnw lines 1275-1292 and 1438-1448",
    "n=153; n2=211; nH0=6691",
    paste0("n=", fmt(n_normal), "; n2=", fmt(n_normal2),
           "; nH0=", fmt(n_normal_h0)),
    "integer equality",
    isTRUE(all.equal(c(n_normal, n_normal2, n_normal_h0),
                     c(153, 211, 6691), tolerance = 1e-12))
)

null <- 0
plocation <- 0
pscale <- 1/sqrt(2)
pdf <- 1
dpm <- 0.5
dpsd1 <- 0
dpsd2 <- 0.1
power <- 0.95
k <- 1/6
nex <- ntbf01(k = k, power = power, null = null, plocation = plocation,
              pscale = pscale, pdf = pdf, alternative = "greater",
              type = "two.sample", dpm = dpm, dpsd = c(dpsd1, dpsd2))
rows[[length(rows) + 1]] <- row(
    "one-sided JZS t sample size",
    "paper/bfssd.Rnw lines 1528-1532 and 1615-1616",
    "143, 195",
    fmt(nex),
    "integer equality",
    isTRUE(all.equal(as.numeric(nex), c(143, 195), tolerance = 1e-12))
)

nseq <- seq(1, 1100, 1)
dpm <- 0.5
null <- 0
k <- 1/6
sd <- 2
psd <- 0.5/sqrt(2)
dpriors <- cbind(dpm = c(0.5, 0.5, 0), dpsd = c(0, 0.1, 0))
power <- 0.95
nH1 <- ceiling(sapply(X = c(1, 2), FUN = function(i) {
    nnmbf01(k = k, power = power, usd = sd, null = null, psd = psd,
            dpm = dpriors[i, 1], dpsd = dpriors[i, 2])
}))
nH0 <- ceiling(sapply(X = 3, FUN = function(i) {
    nnmbf01(k = 1/k, power = power, usd = sd, null = null, psd = psd,
            dpm = dpriors[i, 1], dpsd = dpriors[i, 2],
            lower.tail = FALSE)
}))
nH0normal <- nbf01(k = 6, power = 0.95, usd = sqrt(2), null = 0, pm = 0,
                   psd = 1/sqrt(2), dpm = 0, dpsd = 0,
                   lower.tail = FALSE)
rows[[length(rows) + 1]] <- row(
    "normal-moment sample sizes",
    "paper/bfssd.Rnw lines 1731-1774 and 1796-1802",
    "nH1=302, 429; nH0=997; nH0normal=6691",
    paste0("nH1=", fmt(nH1), "; nH0=", fmt(nH0),
           "; nH0normal=", fmt(nH0normal)),
    "integer equality",
    isTRUE(all.equal(c(nH1, nH0, nH0normal),
                     c(302, 429, 997, 6691), tolerance = 1e-12))
)

## External BFGSD paper examples from SamCH93/bfgsd.
p0H1 <- 0.5
p1H1 <- 0.75
ORH1 <- (p1H1/(1 - p1H1))/(p0H1/(1 - p0H1))
pm <- log(ORH1)
psd <- 0
a1 <- 21
b1 <- 24 - 21
c1 <- 15
d1 <- 26 - 15
logOR1 <- log(a1*d1/(b1*c1))
selogOR1 <- sqrt(1/a1 + 1/b1 + 1/c1 + 1/d1)
bf1 <- bf01(estimate = logOR1, se = selogOR1, null = 0, pm = pm, psd = psd)
bf1_expected <- exp(stats::dnorm(logOR1, mean = 0, sd = selogOR1,
                                 log = TRUE) -
                    stats::dnorm(logOR1, mean = pm, sd = selogOR1,
                                 log = TRUE))
a2 <- 42
b2 <- 50 - 42
c2 <- 30
d2 <- 50 - 30
logOR2 <- log(a2*d2/(b2*c2))
selogOR2 <- sqrt(1/a2 + 1/b2 + 1/c2 + 1/d2)
bf2 <- bf01(estimate = logOR2, se = selogOR2, null = 0, pm = pm, psd = psd)
bf2_expected <- exp(stats::dnorm(logOR2, mean = 0, sd = selogOR2,
                                 log = TRUE) -
                    stats::dnorm(logOR2, mean = pm, sd = selogOR2,
                                 log = TRUE))
rows[[length(rows) + 1]] <- row(
    "Low-PV interim BF01 values",
    "SamCH93/bfgsd paper/BFGSD.R lines 340-379; BFGSD.Rnw lines 950-966",
    paste0("bf1=", fmt(bf1_expected), "; bf2=", fmt(bf2_expected),
           "; reciprocal rounds 9.2 and 27.9"),
    paste0("bf1=", fmt(bf1), " (1/", fmt(1/bf1, 5), "); bf2=", fmt(bf2),
           " (1/", fmt(1/bf2, 5), ")"),
    "1e-12 / rounded reciprocal",
    isTRUE(all.equal(c(bf1, bf2), c(bf1_expected, bf2_expected),
                     tolerance = 1e-12)) &&
        identical(round(c(1/bf1, 1/bf2), 1), c(9.2, 27.9))
)

k1 <- 1/10
k0 <- 10
n <- c(25, 50, 75)
se <- sqrt(1/(p0H1*(1 - p0H1)*n) + 1/(p1H1*(1 - p1H1)*n))
se0 <- sqrt(1/(p0H1*(1 - p0H1)*n) + 1/(p0H1*(1 - p0H1)*n))
bfdesignH1 <- pbf01seq(k1 = k1, k0 = k0, se = se, n = n, pm = pm,
                       psd = psd, dpm = pm, dpsd = 0, type = "normal")
bfdesignH0 <- pbf01seq(k1 = k1, k0 = k0, se = se0, n = n, pm = pm,
                       psd = psd, dpm = 0, dpsd = 0, type = "normal")
lowpv_expected <- c(0.3513772, 0.6701477, 0.8266490,
                    0.01464240, 0.02500717, 0.03001821,
                    48.47064,
                    0.4150487, 0.7304677, 0.8678707,
                    0.01551580, 0.02464182, 0.02853864,
                    45.35815)
lowpv_value <- c(bfdesignH1$cumpH1, bfdesignH1$cumpH0, bfdesignH1$EN,
                 bfdesignH0$cumpH0, bfdesignH0$cumpH1, bfdesignH0$EN)
rows[[length(rows) + 1]] <- row(
    "Low-PV three-look design",
    "SamCH93/bfgsd paper/BFGSD.R lines 392-412; BFGSD.Rnw lines 1005-1025",
    "pinned cumulative H1/H0 probabilities and expected n",
    paste0("H1 cumpH1=", fmt(bfdesignH1$cumpH1), "; H1 cumpH0=",
           fmt(bfdesignH1$cumpH0), "; EN1=", fmt(bfdesignH1$EN),
           "; H0 cumpH0=", fmt(bfdesignH0$cumpH0), "; H0 cumpH1=",
           fmt(bfdesignH0$cumpH1), "; EN0=", fmt(bfdesignH0$EN)),
    "5e-7",
    isTRUE(all.equal(lowpv_value, lowpv_expected, tolerance = 5e-7))
)

powerfun <- Vectorize(FUN = function(nmax, m, H) {
    n <- seq(nmax/m, nmax, length.out = m)
    if (H == "H1") {
        se <- sqrt(1/(p0H1*(1 - p0H1)*n) + 1/(p1H1*(1 - p1H1)*n))
        dpm <- pm
    } else {
        se <- sqrt(1/(p0H1*(1 - p0H1)*n) + 1/(p0H1*(1 - p0H1)*n))
        dpm <- 0
    }
    res <- pbf01seq(k1 = k1, k0 = k0, se = se, n = n, pm = pm, psd = psd,
                    dpm = dpm, dpsd = 0, type = "normal")
    if (H == "H1") res$cumpH1[m] else res$cumpH0[m]
})
m <- 3
pow <- 0.9
rootfun <- function(nmax, m, H) powerfun(nmax = nmax, m = m, H = H) - pow
nmax1 <- uniroot(f = rootfun, interval = c(10, 1000), m = m, H = "H1")$root
nmax0 <- uniroot(f = rootfun, interval = c(10, 1000), m = m, H = "H0")$root
rows[[length(rows) + 1]] <- row(
    "Low-PV maximum sample-size roots",
    "SamCH93/bfgsd paper/BFGSD.R lines 512-514; BFGSD.Rnw lines 1168-1179",
    "ceilings 102 and 87 for 90% correct-evidence probability under H1/H0",
    paste0("nmax1=", fmt(nmax1), " (ceiling ", ceiling(nmax1),
           "); nmax0=", fmt(nmax0), " (ceiling ", ceiling(nmax0), ")"),
    "integer equality",
    identical(c(ceiling(nmax1), ceiling(nmax0)), c(102, 87))
)

n <- seq(40, 100, 10)
plocation <- 0
pscale <- 1/sqrt(2)
pdf <- 1
type <- "two.sample"
alternative <- "greater"
dpm <- 0.5
dpsd <- 0.1
k0 <- 6
k1 <- 1/30
res1 <- ptbf01seq(k1 = k1, k0 = k0, n = n, plocation = plocation,
                  pscale = pscale, pdf = pdf, dpm = dpm, dpsd = dpsd,
                  type = type, alternative = alternative)
res0 <- ptbf01seq(k1 = k1, k0 = k0, n = n, plocation = plocation,
                  pscale = pscale, pdf = pdf, dpm = 0, dpsd = 0,
                  type = type, alternative = alternative)
appendix_expected <- c(
    0.2013309, 0.3033625, 0.3956815, 0.4778737, 0.5509074, 0.6147944,
    0.6690324,
    0.006282381, 0.008682759, 0.010268499, 0.011435333, 0.012141290,
    0.012837484, 0.013228923,
    73.94402,
    0.3092336, 0.4103508, 0.4857032, 0.5446327, 0.5932026, 0.6335704,
    0.6669929,
    0.0008085054, 0.0012818460, 0.0017501755, 0.0020809656,
    0.0023190945, 0.0025378291, 0.0027805852,
    70.12528)
appendix_value <- c(res1$cumpH1, res1$cumpH0, res1$EN1,
                    res0$cumpH0, res0$cumpH1, res0$EN1)
rows[[length(rows) + 1]] <- row(
    "BFGSD appendix one-sided JZS sequential design",
    "SamCH93/bfgsd paper/BFGSD.R lines 823-827; BFGSD.Rnw lines 1701-1704",
    "pinned cumulative H1/H0 probabilities and expected n",
    paste0("H1 cumpH1=", fmt(res1$cumpH1), "; H1 cumpH0=",
           fmt(res1$cumpH0), "; EN1=", fmt(res1$EN1),
           "; H0 cumpH0=", fmt(res0$cumpH0), "; H0 cumpH1=",
           fmt(res0$cumpH1), "; EN0=", fmt(res0$EN1)),
    "5e-7",
    isTRUE(all.equal(appendix_value, appendix_expected, tolerance = 5e-7))
)

out <- do.call(rbind, rows)

cat("| ID | Source | Expected | Package value | Tolerance | Status | Notes |\n")
cat("| --- | --- | --- | --- | --- | --- | --- |\n")
for (i in seq_len(nrow(out))) {
    cat("| ", out$id[i], " | ", out$source[i], " | ", out$expected[i],
        " | ", out$package_value[i], " | ", out$tolerance[i], " | ",
        out$status[i], " | ", out$notes[i], " |\n", sep = "")
}
