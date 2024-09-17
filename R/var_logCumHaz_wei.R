#' @title var_logCumHaz_wei function
#'
#' @description Calculates the variance of log of cumulative excess hazard with
#'  delta method on the log-log of net survival scale.
#'
#' @param object ouput from a model implemented in curesurv
#'
#' @param z_ucured covariates matrix acting on survival function of uncured
#'
#' @param z_pcured covariates matrix acting on cure proportion
#'
#' @param x time at which the estimates are predicted
#'
#' @param cumLexctopred pre prediction obtained with cumLexc_alphaweibull_topred
#'
#' @keywords internal


var_logCumHaz_wei <- function(object, z_ucured =  z_ucured,
                              z_pcured = z_pcured,
                              x = x,cumLexctopred)
{
  if (!inherits(object, "curesurv"))
    stop("Primary argument much be a curesurv object")
  logCum <- dlogCumHazdtheta_wei(z_ucured =  z_ucured, z_pcured = z_pcured,
                                 x = x, theta = object$coefficients,cumLexctopred)

  if(object$pophaz.alpha) {
    var_logCumHaz <- (logCum) %*%
      object$varcov_star[1:(ncol(object$varcov_star)-1),1:(ncol(object$varcov_star)-1)] %*%
      t(logCum)
  }else{
    var_logCumHaz <- (logCum) %*%
      object$varcov_star %*%
      t(logCum)
  }


  return(var_logCumHaz)
}
