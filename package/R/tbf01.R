## Route clearly wrong-tail one-sided statistics through stable quadrature.
.tbf01_tail_quadrature_cutoff <- 4

.tbf01_valid_tail_nquad <- function(tail.nquad) {
    length(tail.nquad) == 1 && is.numeric(tail.nquad) &&
        is.finite(tail.nquad) && tail.nquad >= 2 &&
        tail.nquad == floor(tail.nquad)
}

.tbf01_valid_drange <- function(drange) {
    (is.numeric(drange) && length(drange) == 2 && all(is.finite(drange)) &&
     drange[2] > drange[1]) || (is.character(drange) &&
                                length(drange) == 1 && !is.na(drange) &&
                                drange == "adaptive")
}

.tbf01_pars <- function(n1, n2, type) {
    ## Effective sample size for the noncentrality parameter sqrt(neff)*d.
    if (type == "two.sample") {
        list(df = n1 + n2 - 2, neff = 1/(1/n1 + 1/n2))
    } else {
        list(df = n1 - 1, neff = n1)
    }
}

.tbf01_prior_region <- function(plocation, pscale, pdf, alternative) {
    ## One-sided alternatives truncate and renormalize the analysis prior.
    q0 <- (0 - plocation)/pscale
    if (alternative == "two.sided") {
        list(lower = -Inf, upper = Inf, log_norm_const = 0)
    } else if (alternative == "greater") {
        list(lower = 0, upper = Inf,
             log_norm_const = stats::pt(q = q0, df = pdf, lower.tail = FALSE,
                                        log.p = TRUE))
    } else {
        list(lower = -Inf, upper = 0,
             log_norm_const = stats::pt(q = q0, df = pdf, lower.tail = TRUE,
                                        log.p = TRUE))
    }
}

## exp(x) - 1 - x loses its quadratic term near zero (large df).
.tbf01_expm1_minus_x <- function(x, h = expm1(x)) {
    value <- h - x
    small <- abs(x) < 1e-4
    z <- x[small]
    value[small] <- z^2/2*(1 + z/3*(1 + z/4))
    value
}

## For lambda = sqrt(neff)*SMD and u = t/sqrt(df + t^2), the likelihood
## ratio is E[exp(-lambda^2/2 + u*lambda*Y)], Y ~ chi(df + 1). It remains
## well defined at u = +/-1, the two limits needed by boundary searches.
.tbf01_log_ratio <- function(u, v, df, lambda) {
    if (all(u == 0)) return(-lambda^2/2)
    k <- df + 1
    shift <- u*lambda
    logMode <- asinh(shift/(2*sqrt(k)))
    mode <- sqrt(k)*exp(logMode)
    curvature <- sqrt(k) + mode
    curvature <- curvature*sqrt((sqrt(k)/curvature)^2 + (mode/curvature)^2)

    ## In log(Y), the tilted chi density peaks at Y = mode. Integrate in
    ## units of its local width. The long left panel resolves small df.
    rule <- .bfpwr_gauss_legendre(48)
    nodes <- c(36*rule$x - 44, 8*rule$x - 8, 9*rule$x)
    weights <- c(36*rule$w, 8*rule$w, 9*rule$w)
    x <- outer(1/curvature, nodes)
    h <- expm1(x)
    logDensity <- -k*.tbf01_expm1_minus_x(x, h) - (mode*h)^2/2
    correction <- log(drop(exp(logDensity) %*% weights)/curvature)
    logNorm <- -log(k) - stats::dgamma(k/2, shape = k/2, log = TRUE)

    ## For a positive shift, use v^2 = 1 - u^2 without subtracting two
    ## nearly equal squares. This also covers the polynomial tail at v = 0.
    exponent <- ifelse(shift >= 0, -(lambda*v)^2/2 + (k/2)*(shift/mode),
                       -lambda^2/2 + shift*mode/2)
    k*logMode + exponent + correction - logNorm
}

## A normal prior on the noncentrality gives a scaled noncentral t.
## For a centered half-normal prior, the additional factor is a t CDF.
.tbf01_log_normal_bf10 <- function(u, v, df, m, logG, side) {
    ## g is the normal prior variance. Work on log scale because even a
    ## moderate t prior mixes over a wide range of g.
    logDen <- if (v == 0) 0 else .bfpwr_log1pexp(logG + 2*log(v))
    logVar <- logG - logDen                  # log(g / (1 + g*v^2))
    logBF10 <- df/2*.bfpwr_log1pexp(2*log(abs(u)) + logVar) - logDen/2
    if (m != 0) {
        logScale <- .bfpwr_log1pexp(logG)/2  # log(sqrt(1 + g))
        delta <- m*exp(-logScale)
        logBF10 <- logBF10 + .tbf01_log_ratio(
            u*exp(-logDen/2), if (v == 0) 0 else exp(log(v) + logScale - logDen/2),
            df, delta)
    } else if (side != 0) {
        arg <- side*sign(u)*exp(log(abs(u)) + (log(df + 1) + logVar)/2)
        logBF10 <- logBF10 + log(2) + stats::pt(arg, df + 1, log.p = TRUE)
    }
    logBF10
}

## A t prior is a normal scale mixture with precision W ~ chi-squared(pdf)/pdf.
## Integrate on the log-precision scale to keep the mixing peak resolved as
## pdf increases. Centered priors need only central t probabilities.
.tbf01_log_scale_mixture <- function(u, v, df, m, s, pdf, side, ...) {
    shape <- pdf/2
    rootShape <- sqrt(shape)
    logConstant <- stats::dgamma(shape, shape, log = TRUE) + log(rootShape)
    logDensity <- function(x) {
        y <- x/rootShape
        logG <- 2*log(s) - y
        logConstant - shape*.tbf01_expm1_minus_x(y) +
            .tbf01_log_normal_bf10(u, v, df, m, logG, side)
    }
    centers <- 0
    widths <- 1
    if (v > 0 && pdf < df) {
        ## A large statistic can put a second peak at a much smaller W.
        peak <- rootShape*(2*log(s*v) + log((pdf + 1)/(df - pdf)))
        if (peak < -4) {
            centers <- c(centers, peak)
            widths <- c(widths, sqrt(pdf/(pdf + 1)))
        }
    }
    if (v == 0 && pdf > df) {
        ## The df-th moment weights the lower mixing tail by W^(-df/2).
        centers <- c(centers, rootShape*log1p(-df/pdf))
        widths <- c(widths, sqrt(pdf/(pdf - df)))
    }
    -.tbf01_log_integral(logDensity, centers, widths, ...)
}

.tbf01_log_bf <- function(t, df, neff, plocation, pscale, pdf, region,
                          rel.tol = getOption("bfpwr.rel.tol", 1e-8),
                          abs.tol = rel.tol,
                          subdivisions = getOption("bfpwr.subdivisions", 1000),
                          tail.nquad = getOption("bfpwr.tail.nquad", 512), ...) {
    u <- sign(t)/sqrt(1 + df/t^2)
    v <- if (abs(t) > sqrt(df)) abs(u)*(sqrt(df)/abs(t)) else 1/sqrt(1 + t^2/df)
    m <- sqrt(neff)*plocation
    s <- sqrt(neff)*pscale
    lower <- region$lower
    upper <- region$upper
    side <- if (lower == 0) 1 else if (upper == 0) -1 else 0
    wrongTail <- side != 0 && side*t <= -.tbf01_tail_quadrature_cutoff

    if (is.infinite(pdf) && (side == 0 || m == 0)) {
        return(-.tbf01_log_normal_bf10(u, v, df, m, 2*log(s), side))
    }
    ## On the matching tail the likelihood ratio grows as |lambda|^df.
    ## Its limit averages to infinity when the prior lacks that moment.
    if (is.infinite(t) && pdf <= df && (t == lower || t == upper)) return(-Inf)
    if (m == 0 || (side == 0 && is.infinite(t))) {
        return(.tbf01_log_scale_mixture(u, v, df, m, s, pdf, side,
                   fixed = wrongTail, nquad = tail.nquad, rel.tol = rel.tol,
                   abs.tol = abs.tol, subdivisions = subdivisions, ...))
    }

    ## Partition the effect axis at the prior and likelihood peaks. In each
    ## part an arctangent change of variable puts the peak at one end of a
    ## unit interval and resolves it on its own scale, even for narrow priors.
    width <- sqrt((1 + v^2)/2)/v
    fraction <- stats::plogis(2*(log(s) - log(width)))
    centers <- c(m, (1 - fraction)*m + fraction*t, t)
    widths <- c(s, min(s, width), width)
    if (side != 0) {
        centers <- c(centers, 0)
        widths <- c(widths, min(s, 1)/(1 + abs(u)*sqrt(df + 1)))
    }
    logDensity <- function(d) {
        .tbf01_log_ratio(u, v, df, d) +
            stats::dt((d - m)/s, df = pdf, log = TRUE) - log(s) - region$log_norm_const
    }
    -.tbf01_log_integral(logDensity, centers, widths, lower, upper,
                         fixed = wrongTail, nquad = tail.nquad,
                         rel.tol = rel.tol, abs.tol = abs.tol,
                         subdivisions = subdivisions, ...)
}

.tbf01_log_integral <- function(logDensity, centers, widths, lower = -Inf, upper = Inf,
                                fixed, nquad, ...) {
    centers <- pmax(lower, pmin(upper, centers))
    sorted <- order(centers, widths)
    keep <- sorted[!duplicated(centers[sorted])]
    centers <- centers[keep]
    widths <- widths[keep]
    left <- seq_len(length(centers) - 1L)
    mid <- centers[left] + diff(centers)*widths[left]/(widths[left] + widths[left + 1L])
    angles <- c(atan((c(lower, mid) - centers)/widths),
                atan((c(mid, upper) - centers)/widths))
    keep <- which(angles != 0)
    angles <- angles[keep]
    centers <- rep(centers, 2)[keep]
    widths <- rep(widths, 2)[keep]
    logParts <- function(x) {
        angle <- outer(angles, x)
        d <- centers + widths*tan(angle)
        matrix(logDensity(as.vector(d)), nrow = length(angles)) +
            log(abs(angles)*widths) - 2*log(cos(angle))
    }
    if (!fixed) {
        offset <- max(logParts(c(0, 0.25, 0.5, 0.75, 0.95)))
        integral <- try(stats::integrate(function(x) colSums(exp(logParts(x) - offset)),
                                         0, 1, ...)$value, silent = TRUE)
        if (is.numeric(integral) && is.finite(integral) && integral > 0) {
            return(log(integral) + offset)
        }
    }
    rule <- .bfpwr_gauss_legendre(nquad)
    .bfpwr_logspace_sum(as.vector(sweep(logParts(rule$x), 2, log(rule$w), "+")))
}

tbf01. <- function(t, n, n1 = n, n2 = n, plocation = 0, pscale = 1/sqrt(2),
                   pdf = 1, type = c("two.sample", "one.sample",  "paired"),
                   alternative = c("two.sided", "less", "greater"), log = FALSE,
                   tail.nquad = getOption("bfpwr.tail.nquad", 512),
                   ...) {
    ## input checks
    stopifnot(
        length(t) == 1,
        is.numeric(t),
        is.finite(t),

        length(n1) == 1,
        is.numeric(n1),
        is.finite(n1),
        0 < n1,

        length(n2) == 1,
        is.numeric(n2),
        is.finite(n2),
        0 < n2,

        length(plocation) == 1,
        is.numeric(plocation),
        is.finite(plocation),

        length(pscale) == 1,
        is.numeric(pscale),
        is.finite(pscale),
        0 < pscale,

        length(pdf) == 1,
        is.numeric(pdf),
        0 < pdf,

        length(log) == 1,
        is.logical(log),
        !is.na(log),

        .tbf01_valid_tail_nquad(tail.nquad)
    )
    type <- match.arg(type)
    alternative <- match.arg(alternative)
    .bfpwr_validate_controls(list(...), root = FALSE)
    if (type != "two.sample") {
        if (n1 != n2) {
            warning(paste0('different n1 and n2 supplied but type set to "', type,
                           '", using n = n1'))
        }
    }

    pars <- .tbf01_pars(n1 = n1, n2 = n2, type = type)
    region <- .tbf01_prior_region(plocation = plocation, pscale = pscale,
                                  pdf = pdf, alternative = alternative)

    log_bf <- .tbf01_log_bf(t, pars$df, pars$neff, plocation, pscale, pdf,
                            region, tail.nquad = tail.nquad, ...)

    if (log) return(log_bf)
    else return(exp(log_bf))
}


#' @title t-test Bayes factor
#'
#' @description This function computes the Bayes factor that forms the basis of
#'     the informed Bayesian \eqn{t}-test from Gronau et al. (2020). The Bayes
#'     factor quantifies the evidence that the data provide for the null
#'     hypothesis that the standardized mean difference (SMD) is zero against
#'     the alternative that the SMD is non-zero. A location-scale
#'     \eqn{t}-distribution is assumed for the SMD under the alternative
#'     hypothesis. The Jeffreys-Zellner-Siow (JZS) Bayes factor (Rouder et al.,
#'     2009) is obtained as a special case by setting the location of the prior
#'     to zero and the prior degrees of freedom to one, which is the default.
#'     Infinite prior degrees of freedom give a normal analysis prior.
#'
#'  The data are summarized by \eqn{t}-statistics and sample sizes. The
#'     following types of \eqn{t}-statistics are accepted:
#'
#' - Two-sample \eqn{t}-test where the SMD represents the standardized
#'   mean difference between two group means (assuming equal variances in
#'   both groups)
#' - One-sample \eqn{t}-test where the SMD represents the standardized
#'    mean difference to the null value
#' - Paired \eqn{t}-test where the SMD represents the standardized mean
#'   change score
#'
#' @md
#'
#' @details The Bayes factor averages the likelihood of the observed
#'     \eqn{t}-statistic over the analysis prior, truncated for one-sided
#'     alternatives (Gronau et al., 2020, equation 5). A two-sided normal
#'     prior gives a scaled noncentral \eqn{t} density. Centered \eqn{t}
#'     priors can be integrated as scale mixtures of normals. For shifted
#'     \eqn{t} priors, numerical integration resolves the prior and the
#'     likelihood on their respective scales.
#'
#'     At a fixed sample size, a normal prior or a \eqn{t} prior with more
#'     degrees of freedom than the test statistic can give a positive limiting
#'     Bayes factor as the statistic grows. Some evidence thresholds are
#'     therefore unattainable at small sample sizes.
#'
#' @param t \eqn{t}-statistic
#' @param n Sample size (per group)
#' @param n1 Sample size in group 1 (only required for two-sample \eqn{t}-test
#'     with unequal group sizes)
#' @param n2 Sample size in group 2 (only required for two-sample \eqn{t}-test
#'     with unequal group sizes)
#' @param plocation \eqn{t} prior location. Defaults to \code{0}
#' @param pscale \eqn{t} prior scale. Defaults to \code{1/sqrt(2)}
#' @param pdf \eqn{t} prior degrees of freedom. Defaults to \code{1} (a Cauchy
#'     prior). Set to \code{Inf} for a normal prior with mean \code{plocation}
#'     and standard deviation \code{pscale} before any one-sided truncation.
#' @param type Type of \eqn{t}-test. Can be \code{"two.sample"} (default),
#'     \code{"one.sample"}, or \code{"paired"}
#' @param alternative Direction of the test. Can be either \code{"two.sided"}
#'     (default), \code{"less"}, or \code{"greater"}. The latter two truncate
#'     the analysis prior to negative and positive effects, respectively.
#' @param log Logical indicating whether the natural logarithm of the Bayes
#'     factor should be returned. Defaults to \code{FALSE}
#' @param tail.nquad Number of Gauss-Legendre quadrature nodes used for stable
#'     wrong-tail one-sided calculations. Larger values are more accurate but
#'     slower. Defaults to \code{512} nodes per integration interval.
#' @param ... Additional arguments passed to \code{stats::integrate} for the
#'     direct one-dimensional Bayes factor integral. Defaults are
#'     \code{rel.tol = 1e-8}, \code{abs.tol = rel.tol}, and
#'     \code{subdivisions = 1000}.
#'     Session defaults can be changed with \link{bfpwrOptions}.
#'
#' @inherit bf01 return
#'
#' @author Samuel Pawel, František Bartoš
#'
#' @references Rouder, J. N., Speckman, P. L., Sun, D., Morey, R. D., Iverson,
#'     G. (2009). Bayesian \eqn{t} tests for accepting and rejecting the null
#'     hypothesis. Psychonomic Bulletin & Review, 16(2):225-237.
#'     \doi{10.3758/PBR.16.2.225}
#'
#' Gronau, Q. F., Ly., A., Wagenmakers, E.J. (2020). Informed Bayesian
#'     \eqn{t}-Tests. The American Statistician, 74(2):137-143.
#'     \doi{10.1080/00031305.2018.1562983}
#'
#' @seealso \link{powertbf01}, \link{ptbf01}, \link{ntbf01}
#'
#' @examples
#' ## analyses from Rouder et al. (2009):
#' ## values from Table 1
#' tbf01(t = c(0.69, 3.20), n = 100, pscale = 1, type = "one.sample")
#' ## examples from p. 232
#' tbf01(t = c(2.24, 2.03), n = 80, pscale = 1, type = "one.sample")
#'
#' ## analyses from Gronau et al. (2020) section 3.2:
#' ## informed prior
#' tbf01(t = -0.90, n1 = 53, n2 = 57, plocation = 0.350, pscale = 0.102, pdf = 3,
#'       alternative = "greater", type = "two.sample")
#' ## default (one-sided) prior
#' tbf01(t = -0.90, n1 = 53, n2 = 57, plocation = 0, pscale = 1/sqrt(2), pdf = 1,
#'       alternative = "greater", type = "two.sample")
#'
#' @export
tbf01 <- Vectorize(FUN = tbf01.,
                   vectorize.args = c("t", "n", "n1", "n2", "plocation",
                                      "pscale", "pdf", "type", "alternative",
                                      "log"))

## ## verify results with BayesFactor package
## library(BayesFactor)
## n <- 100
## set.seed(44)
## y <- rnorm(n = n)
## t <- unname(t.test(y)$statistic)
## tbf01(t = t, n = n, r = 1/sqrt(2), type = "one.sample")
## 1/exp(ttest.tstat(t = t, n1 = n, r = 1/sqrt(2))$bf)

## set.seed(4456)
## n1 <- 100
## n2 <- 50
## x <- rnorm(n = n1)
## y <- rnorm(n = n2, mean = 0.5)
## t <- unname(t.test(x, y)$statistic)
## tbf01(t = t, n1 = n1, n2 = n2, r = 1/sqrt(2), type = "two.sample")
## 1/exp(ttest.tstat(t = t, n1 = n1, n2 = n2, r = 1/sqrt(2))$bf)

## set.seed(100)
## n <- 100
## y1 <- rnorm(n = n, mean = 0.5)
## y2 <- rnorm(n = n, mean = 0)
## t <- unname(t.test(y1, y2, paired = TRUE)$statistic)
## tbf01(t = t, n = n, r = 1/sqrt(2), type = "paired")
## 1/exp(ttest.tstat(t = t, n1 = n, r = 1/sqrt(2))$bf)
