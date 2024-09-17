#' dpidtheta_adtneh2 function
#'
#' @description Partial derivatives of pi (net survival at tau) by theta
#'
#' #' @description Partial derivatives of cure fraction (or net survival at tau)
#' by theta from non-mixture model with distribution "tneh".
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
#' @param cumLexctopred pre prediction (obtained from cumLexc_ad2_topred), if NULL then it is calculated
#'
#' @keywords internal

dpidtheta_adtneh2 <- function(z_tau = z_tau,
                              z_alpha = z_alpha,
                              x = x,
                              object,
                              cumLexctopred=NULL) {
  if(is.null(cumLexctopred)){
    cumLexctopred<-cumLexc_ad2_topred(z_tau,z_alpha,x,object$coefficient)
  }

  cumLexc <- cumLexctopred$cumhaz
  pi <- cumLexctopred$pi

  theta <- object$coefficients

  n_z_tau <- ncol(z_tau)
  n_z_alpha <- ncol(z_alpha)
  n_z_tau_ad <- n_z_tau - 1
  n_z_alpha_ad <- n_z_alpha - 1
  alpha0 <- theta[1]
  if (n_z_tau == 0 & n_z_alpha == 0) {

    alpha <- theta[1]
    beta <- theta[2]
    tau <- theta[3]

    aux <- -beta(alpha, beta) * pi
    D <- matrix(0, length(x), length(theta))

    D[, 1] <- aux * tau * (digamma(alpha) - digamma(alpha + beta))
    D[, 2] <- aux * tau * (digamma(beta) - digamma(alpha + beta))
    D[, 3] <- aux
  } else if (n_z_tau > 0 & n_z_alpha > 0) {

    alpha_k <- theta[2:(n_z_alpha + 1)]
    alpha <- alpha0 + z_alpha %*% alpha_k
    beta <- theta[n_z_alpha + 2]
    tau0 <- theta[n_z_alpha + 2 + 1]
    tau_z <- theta[(n_z_alpha + 2 + 1 + 1):(n_z_alpha + 2 + n_z_tau + 1)]
    tau <- tau0 + z_tau %*% tau_z

    aux <- -beta(alpha, beta) * pi
    D <- matrix(0, length(x),length(theta))

    D[, 1] <- aux * tau * (digamma(alpha) - digamma(alpha + beta))
    D[, 2:(n_z_alpha + 1)] <- D[, 1] * z_alpha
    D[, (n_z_alpha + 2)] <- aux * tau * (digamma(beta) - digamma(alpha + beta))
    D[, (n_z_alpha + 3)] <- aux
    D[, (n_z_alpha + 4):(n_z_alpha + 3 + n_z_tau)] <- D[, (n_z_alpha + 3)] * z_tau

  }


  else if (n_z_tau > 0 & n_z_alpha == 0) {


    beta <- theta[n_z_alpha + 2]
    tau0 <- theta[n_z_alpha + 2 + 1]
    tau_z <- theta[(n_z_alpha + 2 + 1 + 1):(n_z_alpha + 2 + n_z_tau + 1)]
    alpha <- alpha0
    tau <- tau0 + z_tau %*% tau_z
    beta2 <- beta

    D <- matrix(0, length(x), length(theta))

    aux <- -beta(alpha, beta2) * pi

    D[, 1] <- aux * tau * (digamma(alpha) - digamma(alpha + beta2))
    D[, (n_z_alpha + 2)] <- aux * tau * (digamma(beta2) - digamma(alpha + beta2))
    D[, (n_z_alpha + 3)] <- aux
    D[, (n_z_alpha + 4):(n_z_alpha + 3 + n_z_tau)] <- D[, (n_z_alpha + 3)] * z_tau

  }


  else if (n_z_tau == 0 & n_z_alpha > 0) {

    alpha_k <- theta[2:(n_z_alpha + 1)]
    alpha <- alpha0 + z_alpha %*% alpha_k
    beta <- theta[n_z_alpha + 2]
    tau <- theta[n_z_alpha + 2 + 1]

    D <- matrix(0, length(x), (n_z_alpha + 3 + n_z_tau))

    aux <- -beta(alpha, beta) * pi

    D[, 1] <- aux * tau * (digamma(alpha) - digamma(alpha + beta))
    D[, 2:(n_z_alpha + 1)] <- D[, 1] * z_alpha
    D[, (n_z_alpha + 2)] <- aux * tau * (digamma(beta) - digamma(alpha + beta))
    D[, (n_z_alpha + 3)] <- aux

  }



  return(D)
}
