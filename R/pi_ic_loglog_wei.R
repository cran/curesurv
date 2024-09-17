#' @title pi_ic_loglog_wei function
#'
#' @description calculates the confidence intervals of the cure proportion pi
#' using variances of log(-log(pi)) by delta method
#'
#'
#' @param object ouput from a model implemented in curesurv
#'
#' @param z_ucured covariates matrix acting on survival function of uncured
#'
#' @param z_pcured covariates matrix acting on cure proportion
#'
#' @param x time at which the estimates are predicted
#'
#' @param level (1-alpha/2)-order quantile of a normal distribution
#'
#' @param cumLexctopred a pre prediction obtained with cumLexc_alphaweibull_topred
#'
#' @param Dpi Partial derivative of Pi calculated by dpidtheta_wei function
#'
#' @keywords internal





pi_ic_loglog_wei <- function(object,
                          z_pcured = z_pcured,
                          z_ucured = z_ucured,
                          x,
                          level,cumLexctopred,Dpi)
{
  n=1+ncol(z_pcured)

  cured <- cumLexctopred$cured
  Dpi_loglog<-sweep(do.call("cbind",Dpi), 1/(log(cured)*cured), MARGIN = 1, '*' )

    varloglogpi <- diag(Dpi_loglog %*%
                       object$varcov_star[1:n,1:n] %*%
                       t(Dpi_loglog))


  upper_bound <-  exp(-exp(log(-log(cured)) - stats::qnorm(level) * sqrt(varloglogpi)))
  lower_bound <-  exp(-exp(log(-log(cured)) + stats::qnorm(level) * sqrt(varloglogpi)))

  IC <- list(t = x, lower_bound = lower_bound, upper_bound = upper_bound)

  return(IC)

}
