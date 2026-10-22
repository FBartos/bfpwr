nmbf01. <- function(estimate, se, null = 0, psd, log = FALSE,
                    alternative = c("two.sided", "less", "greater")) {
    ## input checks
    stopifnot(
        length(estimate) == 1,
        is.numeric(estimate),
        is.finite(estimate),

        length(se) == 1,
        is.numeric(se),
        is.finite(se),
        0 < se,

        length(null) == 1,
        is.numeric(null),
        is.finite(null),

        length(psd) == 1,
        is.numeric(psd),
        is.finite(psd),
        0 < psd,

        length(log) == 1,
        is.logical(log),
        !is.na(log)
    )

    alternative <- match.arg(alternative)
    if (alternative != "two.sided") {
        direction <- if (alternative == "greater") 1 else -1
        logbf <- .nmbf01_log_one_sided(z = direction*(estimate - null)/se,
                                      r = psd/se)
        if (log) return(logbf)
        else return(exp(logbf))
    }

    logbf <- 1.5*log(1 + psd^2/se^2) - 0.5*(estimate - null)^2/se^2/(1 + se^2/psd^2) -
        log(1 + (estimate - null)^2/se^2/(1 + se^2/psd^2))

    if (log) return(logbf)
    else return(exp(logbf))
}


## Point-null vs. positive normal moment prior, in standard-error units.
## With q = z*r/sqrt(1 + r^2), the half-line marginal involves
## I2(q) = integral_0^Inf t^2*exp(q*t - t^2/2) dt.
.nmbf01_log_one_sided <- function(z, r) {
    q <- z/sqrt(1 + 1/r^2)
    if (q < -5) {
        ## The closed form subtracts nearly equal normal-tail terms here.
        ## For a = -q, backward recurrence gives R_j = I_j/I_(j-1) =
        ## j/(a + R_(j+1)), and I0 = 1/(a + R1). Fifty terms suffice for a > 5.
        ratio <- 0
        for (j in 50:1) {
            ratio <- j/(-q + ratio)
            if (j == 2) ratio2 <- ratio
        }
        logIntegral <- log(ratio) + log(ratio2) - log(-q + ratio)
        return(1.5*log1p(r^2) + log(2*pi)/2 - log(2) - logIntegral)
    }
    if (abs(q) < 1e-5) {
        ## Preserve BF01 - 1 when a narrow prior makes q almost zero.
        phi0 <- stats::dnorm(0)
        logHalfMoment <- log1p(q*(4*phi0 + q*(1 + q*2*phi0/3)))
    } else {
        halfMoment <- (1 + q^2)*stats::pnorm(q) + q*stats::dnorm(q)
        logHalfMoment <- log(2*halfMoment)
    }
    1.5*log1p(r^2) - q^2/2 - logHalfMoment
}


#' @title Normal moment prior Bayes factor
#'
#' @description This function computes the Bayes factor that quantifies the
#'     evidence that the data (in the form of an asymptotically normally
#'     distributed parameter estimate with standard error) provide for a point
#'     null hypothesis with a normal moment prior assigned to the parameter
#'     under the alternative.
#'
#' @details A normal moment prior has density \eqn{f(x \mid \code{null},
#'     \code{psd}) = N(x \mid \code{null}, \code{psd}^2) \times (x -
#'     \code{null})^2/ \code{psd}^2}{f(x|\code{null},\code{psd}) =
#'     N(x|\code{null},\code{psd})*(x - \code{null})^2/\code{psd}^2} with
#'     \eqn{N(x \mid m, v)}{N(x|m,v)} the normal density with mean \eqn{m} and
#'     variance \eqn{v} evaluated at \eqn{x}.
#'
#' @param alternative Direction of the alternative hypothesis, one of
#'     \code{"two.sided"} (default), \code{"less"}, or \code{"greater"}.
#'     For a one-sided alternative, the normal moment analysis prior is
#'     restricted to values below \code{null} (\code{"less"}) or above
#'     \code{null} (\code{"greater"}) and renormalized. The point null
#'     hypothesis is unchanged; \code{psd} is the scale before truncation.
#'
#' @inheritParams bf01
#' @param psd Spread of the normal moment prior assigned to the parameter under
#'     the alternative. The modes of the prior are located at
#'     \eqn{\pm\sqrt{2}\,\code{psd}}{+-sqrt(2)*\code{psd}}
#'
#' @inherit bf01 return
#'
#' @author Samuel Pawel
#'
#' @references Johnson, V. E. and Rossell, D. (2010). On the use of non-local
#'     prior densities in Bayesian hypothesis tests. Journal of the Royal
#'     Statistical Society: Series B (Statistical Methodology), 72(2):143–170.
#'     \doi{10.1111/j.1467-9868.2009.00730.x}
#'
#' Pramanik, S. and Johnson, V. E. (2024). Efficient alternatives for Bayesian
#'     hypothesis tests in psychology. Psychological Methods, 29(2):243–261.
#'     \doi{10.1037/met0000482}
#'
#' @examples
#' nmbf01(estimate = 0.25, se = 0.05, null = 0, psd = 0.5/sqrt(2)) # mode at 0.5
#'
#' ## positive normal moment alternative
#' nmbf01(estimate = 0.25, se = 0.05, psd = 0.5/sqrt(2), alternative = "greater")
#'
#' @seealso \link{nmbf01}, \link{pnmbf01}, \link{nnmbf01}, \link{powernmbf01}
#'
#' @export
nmbf01 <- Vectorize(FUN = nmbf01.)
