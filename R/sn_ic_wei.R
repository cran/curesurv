#' @title sn_ic_wei function
#'
#' @description calculates the confidence intervals of the net survival (Sn(t))
#' using variances of Sn(t). In this formula, the variance of net survival
#' is obtained using the expression Var(Sn(t)) = (dSn/dtheta)Var(theta)(dSn/dtheta)^T
#' where Var(theta) is the variance-covariance matrix of theta.
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
#' @param cumLexctopred pre prediction obtained with cumLexc_alphaweibull_topred
#'
#' @param Dsn Partial derivative of Sn calculated by dsndtheta_wei function
#'
#' @keywords internal


sn_ic_wei <- function(object,
                     z_pcured = z_pcured,
                     z_ucured = z_ucured,
                     x,
                     level,
                     cumLexctopred,Dsn)
{
  SurvE <- cumLexctopred$SurvE

  if(object$pophaz.alpha) {
    varsn <- diag(Dsn %*%
      object$varcov_star[1:(ncol(object$varcov_star)-1),1:(ncol(object$varcov_star)-1)] %*%
        t(Dsn))
  }else{
    varsn <- diag(Dsn %*%
      object$varcov_star %*%
      t(Dsn))
  }



  lower_bound <-  SurvE - stats::qnorm(level) * sqrt(varsn)
  upper_bound <-  SurvE + stats::qnorm(level) * sqrt(varsn)


  IC <- list(t = x,
             netsurv = SurvE,
             lower_bound = lower_bound,
             upper_bound = upper_bound)

  return(IC)

}
