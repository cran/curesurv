#' dexhazdtheta_multneh function
#'
#'  @description Partial derivatives of excess hazard
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
#' @param cumLexctopred a pre-prediction parameter, calculated if NULL
#'
#' @keywords internal

dexhazdtheta_multneh <- function(z_tau = z_tau,
                                 z_alpha = z_alpha,
                                 x = x,
                                 object,
                                 cumLexctopred=NULL) {
  theta <- object$coefficients
  if(is.null(cumLexctopred)){
    cumLexctopred<-cumLexc_mul_topred(z_tau,z_alpha,x,theta)
  }

  cumLexc <- cumLexctopred$cumhaz
  ex_haz <- lexc_mul(z_tau, z_alpha, x, theta)

  n_z_tau <- ncol(z_tau)
  n_z_alpha <- ncol(z_alpha)
  n_z_tau_ad <- n_z_tau - 1
  n_z_alpha_ad <- n_z_alpha - 1
  alpha0 <- theta[1]
  if (n_z_tau == 0 & n_z_alpha == 0) {

    alpha <- exp(theta[1])
    beta <- exp(theta[2])+1
    tau <- exp(theta[3])
    u<-x/tau
    D <- matrix(0, length(x), length(theta))


    D[, 1] <- ifelse(x<tau,alpha*log(u)*(u)^(alpha-1)*(1-u)^(beta-1),0)
    u2<-ifelse(x<tau,u,1)
    D[, 2] <- ifelse(x<tau,(beta-1)*log(1-u2)*(u2)^(alpha-1)*(1-u2)^(beta-1),0)
    D[, 3] <- ifelse(x<tau,
                     -((alpha-1)*(u)*(u)^(alpha-2)*(1-u)^(beta-1))+
                       (u)^(alpha-1)*(beta-1)*(u)*(1-u)^(beta-2),
                     0)
  } else if (n_z_tau > 0 & n_z_alpha > 0) {

    alpha_k <- theta[2:(n_z_alpha + 1)]
    alpha <- exp(alpha0 + z_alpha %*% alpha_k)
    beta <- exp(theta[n_z_alpha + 2])+1
    tau0 <- theta[n_z_alpha + 2 + 1]
    tau_z <- theta[(n_z_alpha + 2 + 1 + 1):(n_z_alpha + 2 + n_z_tau + 1)]
    tau <-exp(tau0 + z_tau %*% tau_z)
    u<-x/tau
    D <- matrix(0, length(x),length(theta))

    D[, 1] <- ifelse(x<tau,alpha*log(u)*(u)^(alpha-1)*(1-u)^(beta-1),0)
    D[, 2:(n_z_alpha + 1)] <- D[, 1] * z_alpha
    u2<-ifelse(x<tau,u,1)
    D[, (n_z_alpha + 2)] <- ifelse(x<tau,(beta-1)*log(1-u2)*(u2)^(alpha-1)*(1-u2)^(beta-1),0)

    D[, (n_z_alpha + 3)] <- ifelse(x<tau,
                                   -((alpha-1)*(u)*(u)^(alpha-2)*(1-u)^(beta-1))+
                                     (u)^(alpha-1)*(beta-1)*(u)*(1-u)^(beta-2),
                                   0)
    D[, (n_z_alpha + 4):(n_z_alpha + 3 + n_z_tau)] <- D[, (n_z_alpha + 3)] * z_tau

  }


  else if (n_z_tau > 0 & n_z_alpha == 0) {


    beta <- exp(theta[n_z_alpha + 2])+1
    tau0 <- theta[n_z_alpha + 2 + 1]
    tau_z <- theta[(n_z_alpha + 2 + 1 + 1):(n_z_alpha + 2 + n_z_tau + 1)]
    alpha <- exp(alpha0)
    tau <- exp(tau0 + z_tau %*% tau_z)
    beta2 <- beta
    u<-x/tau
    D <- matrix(0, length(x), length(theta))
    D[, 1] <- ifelse(x<tau,alpha*log(u)*(u)^(alpha-1)*(1-u)^(beta-1),0)
    u2<-ifelse(x<tau,u,1)
    D[, (n_z_alpha + 2)] <- ifelse(x<tau,(beta-1)*log(1-u2)*(u2)^(alpha-1)*(1-u2)^(beta-1),0)
    D[, (n_z_alpha + 3)] <- ifelse(x<tau,
                                   -((alpha-1)*(u)*(u)^(alpha-2)*(1-u)^(beta-1))+
                                     (u)^(alpha-1)*(beta-1)*(u)*(1-u)^(beta-2),
                                   0)
    D[, (n_z_alpha + 4):(n_z_alpha + 3 + n_z_tau)] <- D[, (n_z_alpha + 3)] * z_tau

  }


  else if (n_z_tau == 0 & n_z_alpha > 0) {

    alpha_k <- theta[2:(n_z_alpha + 1)]
    alpha <- exp(alpha0 + z_alpha %*% alpha_k)
    beta <- exp(theta[n_z_alpha + 2])+1
    tau <- exp(theta[n_z_alpha + 2 + 1])
    u<-x/tau
    D <- matrix(0, length(x), (n_z_alpha + 3 + n_z_tau))

    D[, 1] <- ifelse(x<tau,alpha*log(u)*(u)^(alpha-1)*(1-u)^(beta-1),0)
    D[, 2:(n_z_alpha + 1)] <- D[, 1] * z_alpha
    u2<-ifelse(x<tau,u,1)
    D[, (n_z_alpha + 2)] <- ifelse(x<tau,(beta-1)*log(1-u2)*(u2)^(alpha-1)*(1-u2)^(beta-1),0)
    D[, (n_z_alpha + 3)] <- ifelse(x<tau,
                                   -((alpha-1)*(u)*(u)^(alpha-2)*(1-u)^(beta-1))+
                                     (u)^(alpha-1)*(beta-1)*(u)*(1-u)^(beta-2),
                                   0)

  }


  return(D)
}

