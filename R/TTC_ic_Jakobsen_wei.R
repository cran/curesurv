#'@title TTC_ic_Jakobsen_wei function
#'
#' @description calculates the confidence interval of the time TTC using the
#' Jakobsen's approach.
#' Note that this function is for mixture cure model with Weibull distribution
#' considered for uncured patients.
#'
#' @param object ouput from a model implemented in curesurv
#'
#' @param z_ucured covariates matrix acting on survival function of uncured
#'
#' @param z_pcured covariates matrix acting on cure proportion
#'
#'
#' @param epsilon  value fixed by user to estimate the TTC \eqn{\text{Pi}(t)\geq (1-\epsilon)}.
#'   By default \eqn{\epsilon = 0.05}.
#'
#' @param level \code{1-alpha/2}-order quantile of a normal distribution
#'
#' @keywords internal



TTC_ic_Jakobsen_wei <- function(object, z_pcured = z_pcured,
                                z_ucured = z_ucured,
                                epsilon = 0.05,
                                level = 0.975) {
  if (!inherits(object, "curesurv"))
    stop("Primary argument much be a curesurv object")

  time_to_cure_ttc <- TTC_wei(z_pcured = z_pcured,
                              z_ucured = z_ucured,
                              theta = object$coefficients,
                              epsilon = epsilon)



  varTTC <- diag(var_TTC_Jakobsen_wei(object,
                                      z_ucured =  z_ucured,
                                      z_pcured = z_pcured,
                                      epsilon = epsilon))

  lower_bound <-  time_to_cure_ttc - stats::qnorm(level) * sqrt(varTTC)
  upper_bound <-  time_to_cure_ttc + stats::qnorm(level) * sqrt(varTTC)

  IC <- list(time_to_cure_ttc = time_to_cure_ttc,
             lower_bound = lower_bound,
             upper_bound = upper_bound)

  return(IC)


}
