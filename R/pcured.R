#' @title pcured function
#'
#' @description calculates cure fraction of mixture and non mixture cure models.
#'
#' @param object ouput from a model implemented in curesurv
#'
#' @param z_ucured covariates matrix acting on survival function of uncured
#'
#' @param z_pcured covariates matrix acting on cure proportion
#'
#' @param z_alpha Covariates matrix acting on parameter alpha of the density of
#'  time-to-null excess hazard model
#'
#' @param z_tau Covariates matrix acting on time-to-null parameter.
#'
#' @param x time at which the predictions are provided
#'
#' @keywords pcured
#'
#' @return An object of class \code{c("cure_fraction", "data.frame")}.
#' This object is a list containing the following components:
#'
#'
#' \item{time}{time in the input new data}
#'
#' \item{netsurv}{predicted net survival at the time provided in the new data}
#'
#' \item{pi}{pi or net survival at time tau}
#'
#' @keywords internal


pcured <- function(z_pcured = z_pcured,
                   z_ucured = z_ucured,
                   z_tau = z_tau,
                   z_alpha = z_alpha,
                   x = x,
                   object) {

  if (object$model == "mixture") {
    if (object$dist == "weib") {
    netsurv <- exp(-cumLexc_alphaweibull_topred(z_pcured = z_pcured,
                                                z_ucured = z_ucured,
                                                x = x,
                                                object$coefficients)$cumhaz)

    cured <- (cumLexc_alphaweibull_topred(z_pcured = z_pcured,
                                          z_ucured = z_ucured,
                                          x = time,
                                          object$coefficients)$cured)

    pi_and_netsurv <- data.frame(time = x,
                                 netsurv = netsurv,
                                 pi = cured)

    }

    } else if (object$model == "nmixture") {
      if ((object$dist == "tneh" &
           object$link_tau == "loglinear")) {
        netsurv <- cumLexc_mul_topred(z_tau,
                                      z_alpha,
                                      x = time,
                                      object$coefficients)$netsurv
        netsurv_tau <-  cumLexc_mul_topred(z_tau,
                                           z_alpha,
                                           x = time,
                                           object$coefficients)$pi

      }else if ((object$dist == "tneh" &
                 object$link_tau == "linear")) {
        netsurv <- cumLexc_ad_topred(z_tau,
                                     z_alpha,
                                     x = time,
                                     object$coefficients)$netsurv
        netsurv_tau <-  cumLexc_ad_topred(z_tau,
                                          z_alpha,
                                          x = time,
                                          object$coefficients)$pi
      }


      pi_and_netsurv <- data.frame(time = time,
                                   netsurv = netsurv,
                                   pi = netsurv_tau)
    }


    class(pi_and_netsurv) <- c("cure_fraction",
                               "data.frame")
    return(pi_and_netsurv)

}
