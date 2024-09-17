#' @title logCumHaz_ic_wei function
#'
#' @description calculates the confidence intervals of the log of cumulative excess hazard log(-log(Sn(t)))
#' using variances of log(-log(Sn(t))). In this formula, the variance of log of
#' the cumulative excess hazard is obtained using delta method at the scale of log-log of net survival, and with
#' the expression Var(log(log(Sn(t)))) = (dlog(-log(Sn(t)))/dtheta)Var(theta)(dlog(-log(Sn(t)))/dtheta)^T;
#' the Var(theta) is the variance-covariance matrix of theta (estimated parameters).
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
#' @param cumLexctopred a pre prediction obtained with cumLexc_alphaweibull
#'
#' @keywords internal


logCumHaz_ic_wei <- function(object,
                             z_pcured = z_pcured,
                             z_ucured = z_ucured,
                             x,
                             level,cumLexctopred)
{

  cumHazE <- cumLexctopred$cumhaz

  varlogcum <- diag(var_logCumHaz_wei(object,
                                      z_pcured = z_pcured,
                                      z_ucured = z_ucured,
                                      x = x,cumLexctopred
                                      ))
  lower_bound <-  log(cumHazE) - stats::qnorm(level) * sqrt(varlogcum)
  upper_bound <-  log(cumHazE) + stats::qnorm(level) * sqrt(varlogcum)

  IC <- list(t = x, lower_bound = lower_bound, upper_bound = upper_bound)

  return(IC)

}
