#' @title TTC_ic_wei function
#'
#' @description calculates the confidence interval of the time TTC.
#' Note that this function is for mixture cure model with Weibull distribution
#' considered for uncured patients.
#'
#'
#' @param z_ucured covariates matrix acting on survival function of uncured
#'
#'
#' @param z_pcured covariates matrix acting on cure proportion
#'
#'
#' @param object An object of class \code{curesurv}.
#'
#'
#' @param epsilon  value fixed by user to estimate the TTC \eqn{\text{Pi}(t)\geq (1-\epsilon)}.
#'   By default \eqn{\epsilon = 0.05}.
#'
#' @param level \code{1-alpha/2}-order quantile of a normal distribution
#'
#' @param TTC time to cure calculated by TTC_wei
#'
#' @param Dttc partial derivates of TTC by dTTCdtheta_wei
#'
#' @keywords internal


TTC_ic_wei <- function(object, z_pcured = z_pcured,
                       z_ucured = z_ucured, epsilon = 0.05, level,TTC=NULL,Dttc=NULL) {
  if (!inherits(object, "curesurv"))
    stop("Primary argument much be a curesurv object")
  if(is.null(TTC)){
    TTC <- TTC_wei(z_pcured = z_pcured,
                                z_ucured = z_ucured,
                                theta = object$coefficients,
                                epsilon = epsilon)

  }
  if(is.null(Dttc)){
    Dttc<-dTTCdtheta_wei(z_ucured =  z_ucured,
                         z_pcured = z_pcured,
                         theta = object$coefficients,
                         epsilon = epsilon, TTC)
  }

  if(object$pophaz.alpha) {
    varTTC <- diag(Dttc %*%
      object$varcov_star[1:(ncol(object$varcov_star)-1),1:(ncol(object$varcov_star)-1)] %*%
      t(Dttc))
  }else{
    varTTC <- diag(Dttc %*%
      object$varcov_star %*%
      t(Dttc))
  }

  lower_bound <-  TTC - stats::qnorm(level) * sqrt(varTTC)
  upper_bound <-  TTC + stats::qnorm(level) * sqrt(varTTC)

  IC <- list(t = TTC,
             TTC = TTC,
             lower_bound = lower_bound,
             upper_bound = upper_bound,
             varTTC = varTTC)

  return(IC)
}
