#' @title varlogTTC_wei function
#'
#' @description Calculates the variance of \code{log(TTC)}  with delta method.
#' The expression of this variance is expressed as:
#' \code{Var(log(TTC)) = (dlog(TTC)/dtheta)Var(theta)(dlog(TTC)/dtheta)^T}
#'
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
#' @param TTC time to cure previsouly estimated by TTC_wei
#'
#' @keywords internal


varlogTTC_wei <- function(object = object,
                          z_ucured =  z_ucured,
                          z_pcured = z_pcured,
                          epsilon = epsilon,TTC) {


  dlogTTCdtheta <- dlogTTCdtheta_wei(z_ucured =  z_ucured,
                                     z_pcured = z_pcured,
                                     object = object,
                                     epsilon = epsilon,TTC)
  if(object$pophaz.alpha) {
    var_logTTC <- (dlogTTCdtheta) %*% object$varcov_star[1:(ncol(object$varcov_star)-1),1:(ncol(object$varcov_star)-1)] %*% t(dlogTTCdtheta)

  }else{
    var_logTTC <- (dlogTTCdtheta) %*% object$varcov_star %*% t(dlogTTCdtheta)

  }


  return(var_logTTC)
}
