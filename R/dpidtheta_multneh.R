#' dpidtheta_multneh function
#'
#'  @description Partial derivatives of cure fraction (or net survival at tau)
#' by theta from non-mixture model with distribution "tneh" when link_tau="loglinear".
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
#' @param cumLexctopred pre prediction obtained from cumLexc_mul_topred, calculated if NULL
#'
#' @keywords internal


dpidtheta_multneh <- function(z_tau = z_tau,
                              z_alpha = z_alpha,
                              x = x,
                              object,
                              cumLexctopred=NULL) {
  theta <- object$coefficients
  if(is.null(cumLexctopred)){
    cumLexctopred<-cumLexc_mul_topred(z_tau,z_alpha,x,theta)
  }

  cumLexc <- cumLexctopred$cumhaz
  pi <- cumLexctopred$pi



  n_z_tau <- ncol(z_tau)
  n_z_alpha <- ncol(z_alpha)
  n_z_tau_ad <- n_z_tau - 1
  n_z_alpha_ad <- n_z_alpha - 1
  alpha0 <- theta[1]
  if (n_z_tau == 0 & n_z_alpha == 0) {

    alpha <- theta[1]
    beta <- theta[2]
    tau <- exp(theta[3])

    aux <- -tau*beta(alpha, beta)*pi
    D <- matrix(0, length(x), length(theta))

    D[, 1] <- aux* (digamma(alpha) - digamma(alpha + beta))
    D[, 2] <- aux*(digamma(beta) - digamma(alpha + beta))
    D[, 3] <- aux
  } else if (n_z_tau > 0 & n_z_alpha > 0) {

    alpha_k <- theta[2:(n_z_alpha + 1)]
    alpha <- exp(alpha0 + z_alpha %*% alpha_k)
    beta <- exp(theta[n_z_alpha + 2])+1
    tau0 <- theta[n_z_alpha + 2 + 1]
    tau_z <- theta[(n_z_alpha + 2 + 1 + 1):(n_z_alpha + 2 + n_z_tau + 1)]
    tau <- exp(tau0 + z_tau %*% tau_z)

    aux <- -tau*beta(alpha, beta)*pi
    D <- matrix(0, length(x),length(theta))

    D[, 1] <- aux * (digamma(alpha) - digamma(alpha + beta))*alpha
    D[, 2:(n_z_alpha + 1)] <- D[, 1] * z_alpha
    D[, (n_z_alpha + 2)] <- aux * (digamma(beta) - digamma(alpha + beta))*(beta-1)
    D[, (n_z_alpha + 3)] <- aux
    D[, (n_z_alpha + 4):(n_z_alpha + 3 + n_z_tau)] <- D[, (n_z_alpha + 3)] * z_tau

  }


  else if (n_z_tau > 0 & n_z_alpha == 0) {


    beta <- exp(theta[n_z_alpha + 2])+1
    tau0 <- theta[n_z_alpha + 2 + 1]
    tau_z <- theta[(n_z_alpha + 2 + 1 + 1):(n_z_alpha + 2 + n_z_tau + 1)]
    alpha <- exp(alpha0)
    tau <- exp(tau0 + z_tau %*% tau_z)
    beta2 <- beta

    D <- matrix(0, length(x), length(theta))

    aux <- -tau*beta(alpha, beta)*pi

    D[, 1] <- aux * (digamma(alpha) - digamma(alpha + beta2))*alpha
    D[, (n_z_alpha + 2)] <- aux * (digamma(beta2) - digamma(alpha + beta2))*(beta-1)
    D[, (n_z_alpha + 3)] <- aux
    D[, (n_z_alpha + 4):(n_z_alpha + 3 + n_z_tau)] <- D[, (n_z_alpha + 3)] * z_tau

  }


  else if (n_z_tau == 0 & n_z_alpha > 0) {

    alpha_k <- theta[2:(n_z_alpha + 1)]
    alpha <- exp(alpha0 + z_alpha %*% alpha_k)
    beta <- exp(theta[n_z_alpha + 2])+1
    tau <- exp(theta[n_z_alpha + 2 + 1])

    D <- matrix(0, length(x), (n_z_alpha + 3 + n_z_tau))

    aux <- -tau*beta(alpha, beta)*pi

    D[, 1] <- aux * (digamma(alpha) - digamma(alpha + beta))*alpha
    D[, 2:(n_z_alpha + 1)] <- D[, 1] * z_alpha
    D[, (n_z_alpha + 2)] <- aux * (digamma(beta) - digamma(alpha + beta))*(beta-1)
    D[, (n_z_alpha + 3)] <- aux

  }



  return(D)
}
