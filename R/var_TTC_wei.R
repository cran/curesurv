#' @title var_TTC_wei function
#'
#' @description Calculates the variance of TTC with delta method. The expression
#' of this variance is expressed as Var(TTC) = (dTTC/dtheta)Var(theta)(dTTC/dtheta)^T
#' where Var(theta) is the variance-covariance matrix of theta.
#'
#'
#' @param object ouput from a model implemented in curesurv
#'
#' @param z_ucured covariates matrix acting on survival function of uncured
#'
#' @param z_pcured covariates matrix acting on cure proportion
#'
#' @param epsilon  value fixed by user to estimate the TTC \eqn{\text{Pi}(t)\geq (1-\epsilon)}.
#'   By default \eqn{\epsilon = 0.05}.
#'
#' @param TTC time to cure previously calculated using TTC_wei
#'
#' @keywords internal

var_TTC_wei <- function(object,
                        z_ucured =  z_ucured,
                        z_pcured = z_pcured,
                        epsilon = epsilon, TTC)
{
  if (!inherits(object, "curesurv"))
  stop("Primary argument much be a curesurv object")
  dTTC_dtheta <-  dTTCdtheta_wei(z_ucured =  z_ucured,
                                 z_pcured = z_pcured,
                                 theta = object$coefficients,
                                 epsilon = epsilon, TTC)

  if(object$pophaz.alpha) {
    varTTC <- dTTC_dtheta %*%
      object$varcov_star[1:(ncol(object$varcov_star)-1),1:(ncol(object$varcov_star)-1)] %*%
      t(dTTC_dtheta)
  }else{
    varTTC <- dTTC_dtheta %*%
      object$varcov_star %*%
      t(dTTC_dtheta)
  }


  return(varTTC)
}
