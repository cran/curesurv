#' @title dptdtheta_multneh function
#'
#' @description Partial derivatives of probability to be cure by theta from
#' non-mixture model with distribution "tneh".
#'
#'
#' @param object ouput from model implemented in curesurv
#'
#' @param z_alpha Covariates matrix acting on parameter alpha of the density of
#'  time-to-null excess hazard model
#'
#' @param z_tau Covariates matrix acting on time-to-null parameter.
#'
#' @param x time at which the estimates are predicted
#'
#' @param cumLexctopred pre prediction obtained from cumLexc_mul_topred, if NULL will be calculated
#'
#' @param Dpi partial derivative of pi according to theta, if NULL will be calculated
#'
#' @param Dsn partial derivative of net survival according to theta , if NULL will be calculated
#'
#' @keywords internal



dptdtheta_multneh <- function(z_alpha = z_alpha,
                              z_tau = z_tau,
                              x = x,
                              object,
                              cumLexctopred=NULL,
                              Dpi=NULL,
                              Dsn=NULL) {

  if (!inherits(object, "curesurv"))
    stop("Primary argument much be a curesurv object")

  theta <- object$coefficients
  if(is.null(cumLexctopred)){
    cumLexctopred<-cumLexc_mul_topred(z_tau,z_alpha,x,theta)
  }
  if(is.null(Dpi)){
    Dpi <- dpidtheta_multneh(z_tau = z_tau,
                             z_alpha = z_alpha,
                             x = x,
                             object,cumLexctopred)
  }
  if(is.null(Dsn)){
    Dsn <- dsndtheta_multneh(z_tau = z_tau,
                             z_alpha = z_alpha,
                             x = x,
                             object,cumLexctopred,Dpi)
  }

  n_z_tau <- ncol(z_tau)
  n_z_alpha <- ncol(z_alpha)
  n_z_tau_ad <- n_z_tau - 1
  n_z_alpha_ad <- n_z_alpha - 1
  alpha0 <- object$coefficients[1]
  cumLexc <- cumLexctopred$cumhaz
  pi <- cumLexctopred$pi
  sn <- cumLexctopred$netsurv


  D <- matrix(0, nrow = length(x), ncol = length(object$coefficients))


  if (n_z_tau == 0 & n_z_alpha == 0) {

    D[, 1] <- (1/sn^2) * (Dpi[,1] * sn - pi * Dsn[,1])
    D[, 2] <- (1/sn^2) * (Dpi[,2] * sn - pi * Dsn[,2])
    D[, 3] <- (1/sn^2) * (Dpi[,3] * sn - pi * Dsn[,3])

  } else if (n_z_tau > 0 & n_z_alpha > 0) {

    D[, 1] <- (1/sn^2) * (Dpi[,1] * sn - pi * Dsn[,1])
    D[, 2:(n_z_alpha + 1)] <- (1/matrix(rep(sn,n_z_alpha),ncol=n_z_alpha)^2) * (Dpi[,2:(n_z_alpha + 1)] * matrix(rep(sn,n_z_alpha),ncol=n_z_alpha) - matrix(rep(pi,n_z_alpha),ncol=n_z_alpha) * Dsn[,2:(n_z_alpha + 1)])

    D[, (n_z_alpha + 2)] <- (1/sn^2) * (Dpi[,(n_z_alpha + 2)] * sn - pi * Dsn[,(n_z_alpha + 2)])
    D[, (n_z_alpha + 3)] <- (1/sn^2) * (Dpi[,(n_z_alpha + 3)] * sn - pi * Dsn[,(n_z_alpha + 3)])
    D[, (n_z_alpha + 4):(n_z_alpha + 3 + n_z_tau)] <- (1/matrix(rep(sn,n_z_tau),ncol=n_z_tau)^2) * (Dpi[,(n_z_alpha + 4):(n_z_alpha + 3 + n_z_tau)] * matrix(rep(sn,n_z_tau),ncol=n_z_tau) - matrix(rep(pi,n_z_tau),ncol=n_z_tau) * Dsn[,(n_z_alpha + 4):(n_z_alpha + 3 + n_z_tau)]) # Dérivée partielle selon tau*

  }  else if (n_z_tau > 0 & n_z_alpha == 0) {

    D[, 1] <- (1/sn^2) * (Dpi[,1] * sn - pi * Dsn[,1])
    D[, (n_z_alpha + 2)] <- (1/sn^2) * (Dpi[,(n_z_alpha + 2)] * sn - pi * Dsn[,(n_z_alpha + 2)])
    D[, (n_z_alpha + 3)] <- (1/sn^2) * (Dpi[,(n_z_alpha + 3)] * sn - pi * Dsn[,(n_z_alpha + 3)])
    D[, (n_z_alpha + 4):(n_z_alpha + 3 + n_z_tau)] <- (1/sn^2) * (Dpi[,(n_z_alpha + 4):(n_z_alpha + 3 + n_z_tau)] * sn - pi * Dsn[,(n_z_alpha + 4):(n_z_alpha + 3 + n_z_tau)]) # Dérivée partielle selon tau*

  } else if (n_z_tau == 0 & n_z_alpha > 0) {

    D[, 1] <- (1/sn^2) * (Dpi[,1] * sn - pi * Dsn[,1])
    D[, 2:(n_z_alpha + 1)] <- (1/sn^2) * (Dpi[,2:(n_z_alpha + 1)] * sn - pi * Dsn[,2:(n_z_alpha + 1)])
    D[, (n_z_alpha + 2)] <- (1/sn^2) * (Dpi[,(n_z_alpha + 2)] * sn - pi * Dsn[,(n_z_alpha + 2)])
    D[, (n_z_alpha + 3)] <- (1/sn^2) * (Dpi[,(n_z_alpha + 3)] * sn - pi * Dsn[,(n_z_alpha + 3)])

  }

  return(D)
}


