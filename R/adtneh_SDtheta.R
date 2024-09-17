
adtneh_SDtheta <- function(theta,
                           x = x ,
                           d = d,
                           riskpop = riskpop
) {


  n_z_tau <- ncol(z_tau)
  n_z_alpha <- ncol(z_alpha)
  n_z_tau_ad <- n_z_tau - 1
  n_z_alpha_ad <- n_z_alpha - 1
  alpha0 <- (theta[1])
  if (n_z_tau > 0 & n_z_alpha > 0) {


    alpha_k <- theta[2:(n_z_alpha + 1)]
    beta <- (theta[n_z_alpha + 2])
    tau0 <- theta[n_z_alpha + 2 + 1]
    tau_z <- theta[(n_z_alpha + 2 + 1 + 1):(n_z_alpha + 2 + n_z_tau + 1)]
    if (any((alpha0 + z_alpha %*% alpha_k) < 0)) {
      alpha_k <- abs(alpha_k)
    } else if (any((alpha0 + z_alpha %*% alpha_k) == 0)) {
      alpha_k <- abs(alpha_k) + 1e-09
    } else {
      alpha_k  <- alpha_k
    }

    alpha <- (alpha0 + z_alpha %*% alpha_k)
    if (any((tau0 + z_tau %*% tau_z) < 0)) {
      tau_z <- abs(tau_z)
    } else if (any((tau0 + z_tau %*% tau_z) == 0)) {
      tau_z <- abs(tau_z) + 1e-09
    }
    tau <- (tau0 + z_tau %*% tau_z)
    u <- x / (tau)
  } else if (n_z_tau > 0 & n_z_alpha == 0)
  {
    beta <- (theta[2])
    tau0 <- theta[2 + 1]
    tau_z <- theta[(2 + 1 + 1):(2 + n_z_tau + 1)]

    alpha <- (alpha0)
    if (any((tau0 + z_tau %*% tau_z) < 0)) {
      tau_z <- abs(tau_z)
    } else if (any(tau0 + z_tau %*% tau_z == 0)) {
      tau_z <- abs(tau_z) + 1e-09
    }
    tau <- (tau0 + z_tau %*% tau_z)
    u <- x / (tau)


  } else if (n_z_tau == 0 & n_z_alpha > 0) {

    alpha_k <- theta[2:(n_z_alpha + 1)]
    beta <- (theta[n_z_alpha + 2])
    tau0 <- theta[n_z_alpha + 2 + 1]
    if (any((alpha0 + z_alpha %*% alpha_k) < 0)) {
      alpha_k <- abs(alpha_k)
    } else if (any((alpha0 + z_alpha %*% alpha_k) == 0)) {
      alpha_k <- abs(alpha_k) + 1e-09
    }

    alpha <- (alpha0 + z_alpha %*% alpha_k)

    tau <- (tau0)
    u <- x / (tau)


  } else
    if (n_z_tau == 0 & n_z_alpha == 0) {

      beta <- (theta[2])
      tau0 <- theta[3]

      if (alpha0 < 0) {
        alpha <- abs(alpha0)
      }else if (alpha0 == 0) {
        alpha0 <- alpha0 + 1e-09
      }

      if (tau0 < 0) {
        tau0 <- abs(tau0)
      } else if (tau0 == 0) {
        tau0 <- tau0 + 1e-09
      }

      tau <- (tau0)
      u <- x / (tau)
      beta2 <- beta
    }




  n <- length(x)
  npar <- length(theta)
  A <- matrix(0, n, npar)




  u <- ifelse(u > 1, 1, u)
  dLalpha <- sapply(1:length(u),
                    function(i) {
                      u <- u[i]
                      alpha <- alpha[i]
                      func_Lalpha <- function(y) {
                        c(y^(alpha - 1)*(1 - y)^(beta - 1)*log(y))
                      }
                      integrandval <- stats::integrate(func_Lalpha,
                                                0, u)$value
                      return(integrandval)
                    })



  dLbeta <- sapply(1:length(u),
                   function(i) {
                     u <- u[i]
                     alpha <- alpha[i]
                     func_Lbeta <- function(y) {
                       c(y^(alpha - 1)*(1 - y)^(beta - 1)*log(1 - y))
                     }
                     integrandvalb <- stats::integrate(func_Lbeta,
                                                0, u)$value
                     return(integrandvalb)
                   })


  dLtau <- cumLexc_ad2(z_tau,
                      z_alpha,
                      x = x,
                      theta)/tau - stats::dbeta(u,
                                         alpha + 1,
                                         beta) * beta(alpha + 1,
                                                      beta)

  A[,1] <-  ((d * c(log(u))*lexc_ad(z_tau,
                                    z_alpha,
                                    x = x,
                                    theta)) / (riskpop +
                                                 lexc_ad(z_tau,
                                                         z_alpha,
                                                         x = x,
                                                         theta))) -
    tau * dLalpha

  if (n_z_tau > 0 & n_z_alpha > 0) {
    A[,2:(n_z_alpha + 1)] <-  ((z_alpha * c(d *
                                              c(log(u)) *
                                              lexc_ad(z_tau,
                                                      z_alpha,
                                                      x = x,
                                                      theta))) /
                                 c(riskpop + lexc_ad(z_tau,
                                                     z_alpha,
                                                     x = x,
                                                     theta))) - z_alpha * c(tau * dLalpha) #dalphaK

    log_1_u <- ifelse(is.na(log(1 - u)),
                      log(1e-19),
                      log(1 - u))

    A[,(1 + n_z_alpha + 1)] <- ((d * c(log_1_u) *
                                   lexc_ad(z_tau,
                                           z_alpha,
                                           x = x,
                                           theta)) /
                                  (riskpop + lexc_ad(z_tau,
                                                     z_alpha,
                                                     x = x,
                                                     theta))) - tau * dLbeta #dbeta

    A[,(n_z_alpha + 3)] <- ((d * (1/tau) * (2 - alpha - beta +
                                              (beta - 1/1 - u)) *
                               lexc_ad(z_tau,
                                       z_alpha,
                                       x = x,
                                       theta)) /
                              (riskpop + lexc_ad(z_tau,
                                                 z_alpha,
                                                 x = x,
                                                 theta))) - dLtau

    seq_tauz <- (n_z_alpha + 4):(n_z_alpha + n_z_tau + 3)

    A[,seq_tauz] <- ((d * (z_tau/c(tau)) * c(2 - alpha - beta +
                                               (beta - 1/1 - u)) *
                        c(lexc_ad(z_tau,
                                  z_alpha,
                                  x = x,
                                  theta))) /
                       c(riskpop +
                           lexc_ad(z_tau,
                                   z_alpha,
                                   x = x,
                                   theta))) - z_tau * c(dLtau)

  } else if (n_z_tau > 0 & n_z_alpha == 0) {

    log_1_u <- ifelse(is.na(log(1 - u)),
                      log(1e-19),
                      log(1 - u))
    A[,2] <- ((d * c(log_1_u) *
                 lexc_ad(z_tau,
                         z_alpha,
                         x = x,
                         theta)) /
                (riskpop + lexc_ad(z_tau,
                                   z_alpha,
                                   x = x,
                                   theta))) -  tau * dLbeta#dbeta

    A[,3] <- ((d * (1/tau) * c(2 - alpha - beta +
                                 (beta - 1/1 - u)) *
                 lexc_ad(z_tau,
                         z_alpha,
                         x = x,
                         theta)) /
                (riskpop + lexc_ad(z_tau,
                                   z_alpha,
                                   x = x,
                                   theta))) - dLtau

    seq_tauz <- (2 + 1 + 1):(2 + n_z_tau + 1)

    A[,seq_tauz] <- (d * (z_tau/c(tau)) * c(2 - alpha - beta +
                                              (beta - 1/1 - u)) *
                       c(lexc_ad(z_tau,
                                 z_alpha,
                                 x = x,
                                 theta))) /
      c(riskpop +
          lexc_ad(z_tau,
                  z_alpha,
                  x = x,
                  theta)) - z_tau * c(dLtau)


  } else if (n_z_tau == 0 & n_z_alpha > 0) {


    A[,2:(n_z_alpha + 1)] <-  ((d * z_alpha *
                                  c(log(u)) *
                                  c(lexc_ad(z_tau,
                                            z_alpha,
                                            x = x,
                                            theta))) /
                                 c(riskpop + lexc_ad(z_tau,
                                                     z_alpha,
                                                     x = x,
                                                     theta))) -  z_alpha * c(tau * dLalpha) #dalphaK
    log_1_u <- ifelse(is.na(log(1 - u)),
                      log(1e-19),
                      log(1 - u))
    A[,(1 + n_z_alpha + 1)] <- ((d * c(log_1_u) *
                                   lexc_ad(z_tau,
                                           z_alpha,
                                           x = x,
                                           theta)) /
                                  (riskpop + lexc_ad(z_tau,
                                                     z_alpha,
                                                     x = x,
                                                     theta))) -  tau * dLbeta#dbeta

    A[,(n_z_alpha + 3)] <- ((d * (1/tau) * (2 - alpha - beta +
                                              (beta - 1/1 - u)) *
                               lexc_ad(z_tau,
                                       z_alpha,
                                       x = x,
                                       theta)) /
                              (riskpop + lexc_ad(z_tau,
                                                 z_alpha,
                                                 x = x,
                                                 theta))) - dLtau

  } else if (n_z_tau == 0 & n_z_alpha == 0) {


    A[,2] <- ((d * c(log(1 - u)) *
                 lexc_ad(z_tau,
                         z_alpha,
                         x = x,
                         theta)) /
                (riskpop + lexc_ad(z_tau,
                                   z_alpha,
                                   x = x,
                                   theta))) - tau * dLbeta #dbeta

    A[, 3] <- ((d * (1/tau) * (2 - alpha - beta +
                                 (beta - 1/1 - u)) *
                  lexc_ad(z_tau,
                          z_alpha,
                          x = x,
                          theta)) /
                 (riskpop + lexc_ad(z_tau,
                                    z_alpha,
                                    x = x,
                                    theta))) - dLtau
  }

  model_grad <- colSums(A)
  if (anyNA(model_grad)) {
    idna <- which(is.na(model_grad))
    model_grad[idna] <- -1 - 1.0e-19

  }
  return(model_grad)
}
