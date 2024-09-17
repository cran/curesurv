#' @title sn_ic_log_wei function
#'
#' @description calculates the confidence intervals of the net survival (Sn(t))
#' using variances of log(Sn(t)). In this formula, the variance of net survival
#' is obtained using delta method at the scale of log of net survival, and with
#' the expression Var(log(Sn(t))) = (dlog(Sn)/dtheta)Var(theta)(dlog(Sn)/dtheta)^T;
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
#' @param cumLexctopred a pre prediction obtained with cumLexc_alphaweibull_topred
#'
#' @param Dsn Partial derivative of Sn calculated by dsndtheta_wei function
#'
#' @keywords internal





sn_ic_log_wei <- function(object,
                          z_pcured = z_pcured,
                          z_ucured = z_ucured,
                          x,
                          level,cumLexctopred,Dsn)
{

  SurvE <- cumLexctopred$SurvE
  Dsn_log<-sweep(Dsn, 1/SurvE, MARGIN = 1, '*' )

  if(object$pophaz.alpha) {
    varlogsn <- diag(Dsn_log %*%
                    object$varcov_star[1:(ncol(object$varcov_star)-1),1:(ncol(object$varcov_star)-1)] %*%
                    t(Dsn_log))
  }else{
    varlogsn <- diag(Dsn_log %*%
                    object$varcov_star %*%
                    t(Dsn_log))
  }

  lower_bound <-  exp(log(SurvE) - stats::qnorm(level) * sqrt(varlogsn))
  upper_bound <-  exp(log(SurvE) + stats::qnorm(level) * sqrt(varlogsn))

  IC <- list(t = x, lower_bound = lower_bound, upper_bound = upper_bound)

  return(IC)

}
