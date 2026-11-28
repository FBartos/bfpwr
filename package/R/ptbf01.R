ptbf01. <- function(k, n, n1 = n, n2 = n, null = 0, plocation = 0,
                    pscale = 1/sqrt(2), pdf = 1, dpm = plocation, dpsd = pscale,
                    type = c("two.sample", "one.sample", "paired"),
                    alternative = c("two.sided", "less", "greater"),
                    lower.tail = TRUE, drange = "adaptive",
                    tail.eps = getOption("bfpwr.tail.eps", 1e-6),
                    tail.nquad = getOption("bfpwr.tail.nquad", 512), ...) {
    ## input checks
    stopifnot(
        length(k) == 1,
        is.numeric(k),
        is.finite(k),
        0 < k,

        length(n1) == 1,
        is.numeric(n1),
        is.finite(n1),
        1 < n1,

        length(n2) == 1,
        is.numeric(n2),
        is.finite(n2),
        1 < n2,

        length(null) == 1,
        is.numeric(null),
        is.finite(null),

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

        length(dpm) == 1,
        is.numeric(dpm),
        is.finite(dpm),

        length(dpsd) == 1,
        is.numeric(dpsd),
        is.finite(dpsd),
        0 <= dpsd,

        length(lower.tail) == 1,
        is.logical(lower.tail),
        !is.na(lower.tail),

        length(tail.eps) == 1,
        is.numeric(tail.eps),
        is.finite(tail.eps),
        tail.eps > 0,
        tail.eps < 0.5,

        .tbf01_valid_tail_nquad(tail.nquad),

        .tbf01_valid_drange(drange)
    )
    type <- match.arg(type)
    alternative <- match.arg(alternative)
    if (type != "two.sample") {
        if (n1 != n2) {
            warning(paste0('different n1 and n2 supplied but type set to "', type,
                           '", using n = n1'))
            n2 <- n1
        }
    }

    ## Compute the same t boundaries as a single sequential look, then
    ## integrate the known-variance predictive distribution of the estimate.
    dots <- list(...)
    .bfpwr_validate_controls(dots)
    se <- 1/sqrt(.tbf01_pars(n1, n2, type)$neff)
    estsd <- sqrt(se^2 + dpsd^2)
    trange <- if (is.numeric(drange)) (drange - null)/se else drange
    limits <- if (alternative != "two.sided" && !is.numeric(trange)) {
        .bfpwr_one_sided_tail_limits(origin = 0, step_scale = 1,
            mean = (dpm - null)/se, sd = estsd/se, tail.eps = tail.eps)
    }
    result <- do.call(.bfpwr_tcrit_result, c(list(
        k = k, n1 = n1, n2 = n2, plocation = plocation - null,
        pscale = pscale, pdf = pdf, type = type, alternative = alternative,
        trange = trange, search_limit = limits, tail.nquad = tail.nquad
    ), dots))
    if (result$status == "impossible") {
        if (alternative == "two.sided") return(as.numeric(lower.tail))
        warning("BF01 = k appears unattainable for this one-sided t test")
        bf0 <- do.call(tbf01, c(list(t = 0, n1 = n1, n2 = n2, type = type,
            alternative = alternative, plocation = plocation - null, pscale = pscale,
            pdf = pdf, log = TRUE, tail.nquad = tail.nquad), .bfpwr_integrate_dots(dots)))
        return(as.numeric((bf0 < log(k)) == lower.tail))
    }
    if (result$status == "tail_cutoff") {
        warning(paste0("Adaptive t power-boundary search reached the predictive ",
                       "tail cutoff (absolute error <= ", format(tail.eps),
                       "). Pass a numeric 'drange' interval to search exact bounds."))
    } else if (result$status != "ok") {
        warning("Numerical problems finding critical value")
        return(NaN)
    }
    bounds <- null + se*result$value
    if (alternative != "two.sided") {
        return(stats::pnorm(bounds, dpm, estsd,
                            lower.tail = (alternative == "less") == lower.tail))
    }
    if (!lower.tail) return(exp(.bfpwr_lpnorm_interval(bounds[1], bounds[2], dpm, estsd)))
    exp(.bfpwr_logspace_sum(c(
        stats::pnorm(bounds[1], dpm, estsd, log.p = TRUE),
        stats::pnorm(bounds[2], dpm, estsd, lower.tail = FALSE, log.p = TRUE)
    )))
}


#' @title Cumulative distribution function of the t-test Bayes factor
#'
#' @description This function computes the probability of obtaining a
#'     \eqn{t}-test Bayes factor (\link{tbf01}) more extreme than a threshold
#'     \code{k} with a specified sample size. The future standardized estimate
#'     is modeled with its known-variance normal distribution. Thus, although
#'     the analysis Bayes factor uses the finite-sample \eqn{t} distribution,
#'     the returned power is an asymptotic approximation rather than exact
#'     noncentral-\eqn{t} power; differences can be material at small sample
#'     sizes.
#'
#' @inheritParams tbf01
#' @inheritParams pbf01
#' @param null Standardized mean difference under the point null hypothesis.
#'     Defaults to \code{0}
#' @param alternative Direction of the test. Can be either \code{"two.sided"}
#'     (default), \code{"less"}, or \code{"greater"}. The latter two truncate
#'     the analysis prior to effects below and above \code{null}, respectively.
#' @param dpm Mean of the normal design prior assigned to the standardized mean
#'     difference. Defaults to the analysis prior location
#' @param dpsd Standard deviation of the normal design prior assigned to the
#'     standardized mean difference. Set to \code{0} to obtain a point prior at
#'     the design prior mean. Defaults to the analysis prior scale
#' @param drange Standardized mean difference search range over which the
#'     critical values are searched for. Can be either set to a numerical range
#'     or to \code{"adaptive"} (default) which determines the range in an
#'     adaptive way from the other input parameters. For a two-sided test, a
#'     numeric interval supplies initial brackets that can be extended to a
#'     crossing. A one-sided numeric search stays inside the supplied interval.
#' @param tail.eps One-sided adaptive searches stop once the remaining
#'     predictive probability in the searched tail is at most this value. If no
#'     boundary is found before that point, the boundary-free probability
#'     implied by the searched side is returned and the omitted tail probability
#'     is bounded by \code{tail.eps}. Smaller values search farther and can
#'     recover extremely remote boundaries at additional computational cost.
#'     Defaults to \code{1e-6}
#' @param tail.nquad Number of Gauss-Legendre quadrature nodes used by
#'     \code{\link{tbf01}} for stable wrong-tail one-sided calculations. Larger
#'     values are more accurate but slower. Defaults to \code{512}.
#' @param ... Optional numerical controls for all boundary searches.
#'     \code{rel.tol} (default \code{1e-8}), \code{abs.tol} (default
#'     \code{rel.tol}), and \code{subdivisions} (default \code{1000}),
#'     together with \code{stop.on.error} and \code{keep.xy}, are passed to
#'     \code{stats::integrate} for BF evaluation. Root searches use
#'     \code{tol} (default \code{1e-8}), \code{maxiter}, \code{trace}, and
#'     \code{check.conv} from \code{stats::uniroot}. One-sided scouting uses
#'     a looser integral unless \code{rel.tol} is supplied; returned roots are
#'     checked against the final integral with the requested accuracy.
#'     Session defaults can be changed with \link{bfpwrOptions}.
#'     These controls reduce numerical error within the normal predictive
#'     approximation; they do not remove that approximation.
#'
#' @inherit pbf01 return
#'
#' @author Samuel Pawel, František Bartoš
#'
#' @seealso \link{tbf01}, \link{ntbf01}, \link{powertbf01}
#'
#' @examples
#' ## example from Schönbrodt and Wagenmakers (2018, p. 135)
#' ptbf01(k = 1/6, n = 146, dpm = 0.5, dpsd = 0, alternative = "greater")
#' ptbf01(k = 6, n = 146, dpm = 0, dpsd = 0, alternative = "greater",
#'        lower.tail = FALSE)
#'
#' ## two-sided
#' ptbf01(k = 1/6, n = 146, dpm = 0.5, dpsd = 0)
#' ptbf01(k = 6, n = 146, dpm = 0, dpsd = 0, lower.tail = FALSE)
#'
#' ## one-sample test
#' ptbf01(k = 1/6, n = 146, dpm = 0.5, dpsd = 0, alternative = "greater",
#'        type = "one.sample")
#'
#' @export
ptbf01 <- Vectorize(FUN = ptbf01.,
                    vectorize.args = c("k", "n", "n1", "n2", "null",
                                       "plocation", "pscale", "pdf", "type",
                                       "alternative", "dpm", "dpsd",
                                       "lower.tail"))

## ## verify with simulation
## set.seed(10000)
## nsim <- 10000
## dpm <- 0.5
## dpsd <- 0.1
## n1 <- 30
## n2 <- 120
## r <- 1/sqrt(2)
## plocation <- 0.3
## bfs <- replicate(n = nsim, expr = {
##     if (dpsd == 0) m <- dpm
##     else m <- rnorm(n = 1, mean = dpm, sd = dpsd)
##     x <- rnorm(n = n1, mean = 0, sd = 1)
##     y <- rnorm(n = n2, mean = m, sd = 1)
##     t <- t.test(y, x)$statistic
##     ## se <- sqrt(1/n1 + 1/n2)
##     ## est <- rnorm(n = 1, mean = dpm, sd = sqrt(se^2 + dpsd^2))
##     ## t <- est/se
##     tbf01(t = t, n1 = n1, n2 = n2, pscale = r, plocation = plocation,
##           type ="two.sample", alternative = "two.sided")
## })
## mean(bfs < 1/10)
## ptbf01(k = 1/10, n1 = n1, n2 = n2, pscale = r, plocation = plocation, dpm = dpm,
##        dpsd = dpsd, type = "two.sample", alternative = "two.sided", drange = "adaptive")
