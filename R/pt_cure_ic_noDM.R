#' @title pt_cure_ic_noDM function
#' @description Variance of p(t) without delta method
#' @param x fit ouput from curesurv
#' @param z_pcured covariates
#' @param t time
#' @param level level of confidence
#' @param theta estimated parameters
#'
#' @keywords internal


pt_cure_ic_noDM = function(x, z_pcured, t, level=0.975, theta)
{

  z <- stats::qnorm(level)
  lower_bound <-  x$pt_cure - z * sqrt(var(x$pt_cure))
  upper_bound <-  x$pt_cure + z * sqrt(var(x$pt_cure))


  IC <- list(t = t,
             lower_bound = lower_bound,
             upper_bound = upper_bound)

  return(IC)

}
