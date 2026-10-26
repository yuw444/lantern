# ============================================================================
# Cauchy combination test
# ============================================================================

#' Cauchy combination test for multiple p-values
#'
#' Combines \eqn{K} p-values from correlated or independent tests into a single
#' meta-analysis p-value using the Cauchy combination test (Liu & Xie 2020).
#' Unlike Fisher's method, the test remains valid under arbitrary dependency
#' structures among p-values.
#'
#' @param p_values Numeric vector of p-values (\eqn{K \ge 2}).
#'   Values must be in \eqn{(0, 1]}.  Exact zeros are clipped to
#'   \code{.Machine$double.eps}.
#' @param weights Optional non-negative numeric weight vector of the same
#'   length as \code{p_values}.  Defaults to equal weights.
#'
#' @return A single numeric p-value in \eqn{(0, 1]}.
#'
#' @references
#' Liu, Y. and Xie, J. (2020). Cauchy Combination Test: A Powerful Test With
#' Analytic p-Value Calculation Under Arbitrary Dependency Structures.
#' \emph{J. Am. Stat. Assoc.} \strong{115}(529), 393--402.
#' \doi{10.1080/01621459.2018.1554485}
#'
#' @examples
#' # Combine African-ancestry and European-ancestry p-values for one gene
#' cauchy_combine(c(0.01, 0.04))
#'
#' # Three-population meta-analysis
#' cauchy_combine(c(0.01, 0.04, 0.20))
#'
#' # Weighted combination (double weight on AFR)
#' cauchy_combine(c(0.01, 0.04), weights = c(2, 1))
#'
#' # Apply across a results table (one gene per row)
#' # mapply(cauchy_combine, split(cbind(p_aa, p_ee), seq_len(nrow(res))))
#'
#' @seealso \code{\link{write_dosage_gds}}, \code{\link{ancestry_split}},
#'   \code{\link{ancestry_smmat}}
#'
#' @export
cauchy_combine <- function(p_values, weights = NULL) {
  if (!is.numeric(p_values)) stop("p_values must be numeric")
  if (!is.null(weights) && !is.numeric(weights)) stop("weights must be numeric")
  p <- as.numeric(p_values)
  if (length(p) < 2) stop("p_values must have length >= 2")
  if (any(is.na(p)))  stop("p_values contains NA")
  if (any(p < 0 | p > 1)) stop("p_values must be in (0, 1]")
  p <- pmax(p, .Machine$double.eps)

  if (is.null(weights)) {
    T_stat <- mean(tan((0.5 - p) * pi))
  } else {
    w <- as.numeric(weights)
    if (length(w) != length(p)) stop("weights must be the same length as p_values")
    if (any(w < 0))             stop("weights must be non-negative")
    w      <- w / sum(w)
    T_stat <- sum(w * tan((0.5 - p) * pi))
  }

  pcauchy(T_stat, lower.tail = FALSE)
}
