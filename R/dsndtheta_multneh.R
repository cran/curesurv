#' @title dsndtheta_multneh function
#'
#' @description Partial derivatives of sn (net survival) by theta for the link_tau="loglinear" model
#'
#' @param object ouput from a model implemented in curesurv
#'
#' @param z_alpha Covariates matrix acting on parameter alpha of the density of
#'  time-to-null excess hazard model
#'
#' @param z_tau Covariates matrix acting on time-to-null parameter.
#'
#' @param x time at which the predictions are provided
#'
#' @param cumLexctopred pre prediction obtained from cumLexc_mul_topred, estimated if NULL
#'
#' @param Dpi partial derivative of pi according to theta, if NULL will be calculated
#'
#' @author Juste Goungounga, Judith Breaud, Eugenie Blandin, Olayidé Boussari, Valérie Jooste
#'
#' @references Boussari O, Bordes L, Romain G, Colonna M, Bossard N, Remontet L,
#'  Jooste V. Modeling excess hazard with time-to-cure as a parameter.
#'   Biometrics. 2021 Dec;77(4):1289-1302. doi: 10.1111/biom.13361.
#'    Epub 2020 Sep 12. PMID: 32869288.
#' (\href{https://pubmed.ncbi.nlm.nih.gov/32869288/}{pubmed})
#'
#'
#'  Boussari O, Romain G, Remontet L, Bossard N, Mounier M, Bouvier AM,
#'  Binquet C, Colonna M, Jooste V. A new approach to estimate time-to-cure from
#'  cancer registries data. Cancer Epidemiol. 2018 Apr;53:72-80.
#'  doi: 10.1016/j.canep.2018.01.013. Epub 2018 Feb 4. PMID: 29414635.
#' (\href{https://pubmed.ncbi.nlm.nih.gov/29414635/}{pubmed})
#'
#' @keywords internal

dsndtheta_multneh <- function(z_tau, z_alpha, x, object,cumLexctopred=NULL,Dpi=NULL) {
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
  cumLexc <- cumLexctopred$cumhaz
  sn <- cumLexctopred$netsurv

  n_z_tau <- ncol(z_tau)
  n_z_alpha <- ncol(z_alpha)
  n_z_tau_ad <- n_z_tau - 1
  n_z_alpha_ad <- n_z_alpha - 1
  alpha0 <- theta[1]
  D <- matrix(0, length(x), length(theta))
  Dpi <- dpidtheta_multneh(z_tau = z_tau,
                           z_alpha = z_alpha,
                           x = x,
                           object,cumLexctopred)

  if (n_z_tau == 0 & n_z_alpha == 0) {

    alpha <- exp(theta[1])
    beta <- exp(theta[2])+1
    tau <- exp(theta[3])

    u <- x / (tau)

    aux <- -sn*tau

    funca <- function(y,i) {
      y ^ (alpha - 1) * (1 - y) ^ (beta - 1) * log(y)
    }
    funcb <- function(y,i) {
      y ^ (alpha - 1) * (1 - y) ^ (beta - 1) * log(1 - y)
    }
    dLa <- sapply(1:length(u), function(i) {
      if (0<u[i]&u[i] < 1) {
        integrate(funca, 0, u[i], i)$value
      } else {
        NA
      }
    })

    dLb <- sapply(1:length(u), function(i) {
      if (0<u[i]&u[i] < 1) {
        integrate(funcb, 0, u[i], i)$value
      } else {
        NA
      }
    })

    dLT<-pbeta(u,alpha,beta)*beta(alpha,beta) - dbeta(u, alpha + 1, beta) * beta(alpha + 1, beta)

    D[, 1] <- ifelse(u<1,(aux * alpha * dLa),Dpi[,1])
    D[, 2] <- ifelse(u<1,(aux * (beta-1) * dLb),Dpi[,2])
    D[, 3] <- ifelse(u<1,(aux * dLT),Dpi[,3])


  }


  else if (n_z_tau > 0 & n_z_alpha > 0) {

    alpha_k <- theta[2:(n_z_alpha + 1)]
    alpha <- exp(alpha0 + z_alpha %*% alpha_k)
    beta <- exp(theta[n_z_alpha + 2])+1
    tau0 <- theta[n_z_alpha + 2 + 1]
    tau_z <- theta[(n_z_alpha + 2 + 1 + 1):(n_z_alpha + 2 + n_z_tau + 1)]
    tau <- exp(tau0 + z_tau %*% tau_z)

    u <- x / (tau)

    aux <- -sn*tau

    funca <- function(y, i) {
      y ^ (alpha[i] - 1) * (1 - y) ^ (beta - 1) * log(y)
    }

    funcb <- function(y, i) {
      y ^ (alpha[i] - 1) * (1 - y) ^ (beta - 1) * log(1 - y)
    }

    dLa <- sapply(1:length(u), function(i) {
      if (0<u[i]&u[i] < 1) {
        integrate(funca, 0, u[i], i)$value
      } else {
        NA
      }
    })

    dLb <- sapply(1:length(u), function(i) {
      if (0<u[i]&u[i] < 1) {
        integrate(funcb, 0, u[i], i)$value
      } else {
        NA
      }
    })
    dLT <- pbeta(u,alpha,beta)*beta(alpha,beta) - dbeta(u, alpha + 1, beta) * beta(alpha + 1, beta)
    #D
    D[, 1] <- ifelse(u<1,(aux * alpha * dLa),Dpi[, 1])
    D[, 2:(n_z_alpha + 1)] <- (D[, 1] * z_alpha)
    D[, (n_z_alpha + 2)] <- ifelse(u<1,(aux * (beta-1) * dLb),Dpi[, (n_z_alpha + 2)])
    D[, (n_z_alpha + 3)] <- ifelse(u<1,(aux * dLT),Dpi[, (n_z_alpha + 3)])
    D[, (n_z_alpha + 4):(n_z_alpha + 3 + n_z_tau)] <- (D[, (n_z_alpha + 3)] * z_tau)


  }


  else if (n_z_tau > 0 & n_z_alpha == 0) {


    beta <- exp(theta[n_z_alpha + 2])+1
    tau0 <- theta[n_z_alpha + 2 + 1]
    tau_z <- theta[(n_z_alpha + 2 + 1 + 1):(n_z_alpha + 2 + n_z_tau + 1)]
    alpha <- exp(alpha0)
    tau <- exp(tau0 + z_tau %*% tau_z)

    u <- x / (tau)


    aux <- -sn*tau


    funca <- function(y,i) {
      y ^ (alpha - 1) * (1 - y) ^ (beta - 1) * log(y)
    }
    funcb <- function(y,i) {
      y ^ (alpha - 1) * (1 - y) ^ (beta - 1) * log(1 - y)
    }


    dLa <- sapply(1:length(u), function(i) {
      if (0<u[i]&u[i] < 1) {
        integrate(funca, 0, u[i], i)$value
      } else {
        NA
      }
    })

    dLb <- sapply(1:length(u), function(i) {
      if (0<u[i]&u[i] < 1) {
        integrate(funcb, 0, u[i], i)$value
      } else {
        NA
      }
    })
    dLT <- pbeta(u,alpha,beta)*beta(alpha,beta) - dbeta(u, alpha + 1, beta) * beta(alpha + 1, beta)
    D[, 1] <- ifelse(u<1,(aux * alpha * dLa),Dpi[,1])
    D[, (n_z_alpha + 2)] <- ifelse(u<1,(aux * (beta-1) * dLb),Dpi[,(n_z_alpha + 2)])
    D[, (n_z_alpha + 3)] <- ifelse(u<1,(aux * dLT),Dpi[, (n_z_alpha + 3)])
    D[, (n_z_alpha + 4):(n_z_alpha + 3 + n_z_tau)] <- (D[, (n_z_alpha + 3)] * z_tau)

  }


  else if (n_z_tau == 0 & n_z_alpha > 0) {

    alpha_k <- theta[2:(n_z_alpha + 1)]
    alpha <- exp(alpha0 + z_alpha %*% alpha_k)
    beta <- exp(theta[n_z_alpha + 2])+1
    tau <- exp(theta[n_z_alpha + 2 + 1])

    u <- x / (tau)

    aux <- -sn*tau
    funca <- function(y, i) {
      y ^ (alpha[i] - 1) * (1 - y) ^ (beta - 1) * log(y)
    }

    funcb <- function(y, i) {
      y ^ (alpha[i] - 1) * (1 - y) ^ (beta - 1) * log(1 - y)
    }


    dLa <- sapply(1:length(u), function(i) {
      if (0<u[i]&u[i] < 1) {
        integrate(funca, 0, u[i], i)$value
      } else {
        NA
      }
    })

    dLb <- sapply(1:length(u), function(i) {
      if (0<u[i]&u[i] < 1) {
        integrate(funcb, 0, u[i], i)$value
      } else {
        NA
      }
    })
    dLT <- pbeta(u,alpha,beta)*beta(alpha,beta) - dbeta(u, alpha + 1, beta) * beta(alpha + 1, beta)
    D[, 1] <- ifelse(u<1,(aux * alpha * dLa),Dpi[, 1])
    D[, 2:(n_z_alpha + 1)] <- (D[, 1] * z_alpha)
    D[, (n_z_alpha + 2)] <- ifelse(u<1,(aux * (beta-1) * dLb),Dpi[, (n_z_alpha + 2)])
    D[, (n_z_alpha + 3)] <- ifelse(u<1,(aux * dLT),Dpi[, (n_z_alpha + 3)])

  }

  return(D)
}
