
adtneh_vartheta <- function(z_tau, z_alpha,
                            time, theta,
                            d, riskpop) {
  x <- time
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
  if (any(alpha0 + z_alpha %*% alpha_k < 0)) {
    alpha_k <- abs(alpha_k)
  } else if (any(alpha0 + z_alpha %*% alpha_k == 0)) {
    alpha_k <- abs(alpha_k) + 1e-09
  }

  alpha <- (alpha0 + z_alpha %*% alpha_k)
  if (any(tau0 + z_tau %*% tau_z < 0)) {
    tau_z <- abs(tau_z)
  } else if (any(tau0 + z_tau %*% tau_z == 0)) {
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
      if (any(tau0 + z_tau %*% tau_z < 0)) {
        tau_z <- abs(tau_z)
      } else if (any(tau0 + z_tau %*% tau_z == 0)) {
        tau_z <- abs(tau_z) + 1e-09
      }
      tau <- (tau0 + z_tau %*% tau_z)
      u <- x / (tau)


  } else if (n_z_tau == 0 & n_z_alpha > 0){

    alpha_k <- theta[2:(n_z_alpha + 1)]
    beta <- (theta[n_z_alpha + 2])
    tau0 <- theta[n_z_alpha + 2 + 1]
    if (any(alpha0 + z_alpha %*% alpha_k < 0)) {
      alpha_k <- abs(alpha_k)
    } else if (any(alpha0 + z_alpha %*% alpha_k == 0)) {
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
    u <- time / (tau)
    beta2 <- beta
  }




    n <- length(time)
    npar <- length(theta)
    A <- matrix(0, n, npar)

    A[,1] <-  (d * c(log(u))*lexc_ad(z_tau,
                                   z_alpha,
                                   x = time,
                                   theta)) / (riskpop +
                                                lexc_ad(z_tau,
                                                        z_alpha,
                                                        x = time,
                                                        theta))


    if (n_z_tau > 0 & n_z_alpha > 0) {
      A[,2:(n_z_alpha + 1)] <-  (d * z_alpha *
                                   c(log(u)) *
                                   lexc_ad(z_tau,
                                           z_alpha,
                                           x = time,
                                           theta)) /
                                (riskpop + lexc_ad(z_tau,
                                                   z_alpha,
                                                   x = time,
                                                   theta))
      A[,(1 + n_z_alpha + 1)] <- (d * c(log(1 - u)) *
                                    lexc_ad(z_tau,
                                            z_alpha,
                                            x = time,
                                            theta)) /
                                 (riskpop + lexc_ad(z_tau,
                                                    z_alpha,
                                                    x = time,
                                                    theta))

      A[,(1 + n_z_alpha + 1 + 1)] <- (d * (1/tau) * (2 - alpha - beta +
                                                      (beta - 1/1 - u)) *
                                  lexc_ad(z_tau,
                                          z_alpha,
                                          x = time,
                                          theta)) /
                                  (riskpop + lexc_ad(z_tau,
                                                     z_alpha,
                                                     x = time,
                                                     theta))

      seq_tauz <- (n_z_alpha + 2 + 1 + 1):(n_z_alpha + 2 + n_z_tau + 1)

      A[,seq_tauz] <- (d * (z_tau/tau) * (2 - alpha - beta +
                                            (beta - 1/1 - u)) *
                                   lexc_ad(z_tau,
                                           z_alpha,
                                           x = time,
                                           theta)) /
                                    (riskpop +
                                       lexc_ad(z_tau,
                                               z_alpha,
                                               x = time,
                                               theta))
    } else if (n_z_tau > 0 & n_z_alpha == 0) {


      A[,2] <- (d * c(log(1 - u)) *
                                    lexc_ad(z_tau,
                                            z_alpha,
                                            x = time,
                                            theta)) /
        (riskpop + lexc_ad(z_tau,
                           z_alpha,
                           x = time,
                           theta))

      A[,3] <- (d * (1/tau) * (2 - alpha - beta +
                                 (beta - 1/1 - u)) *
                  lexc_ad(z_tau,
                          z_alpha,
                          x = time,
                          theta)) /
        (riskpop + lexc_ad(z_tau,
                           z_alpha,
                           x = time,
                           theta))


      seq_tauz <- (2 + 1 + 1):(2 + n_z_tau + 1)

      A[,seq_tauz] <- (d * (z_tau/tau) * (2 - alpha - beta +
                                            (beta - 1/1 - u)) *
                         lexc_ad(z_tau,
                                 z_alpha,
                                 x = time,
                                 theta)) /
        (riskpop +
           lexc_ad(z_tau,
                   z_alpha,
                   x = time,
                   theta))


    } else if (n_z_tau == 0 & n_z_alpha > 0) {

      A[,2:(n_z_alpha + 1)] <-  (d * z_alpha *
                                   c(log(u)) *
                                   lexc_ad(z_tau,
                                           z_alpha,
                                           x = time,
                                           theta)) /
        (riskpop + lexc_ad(z_tau,
                           z_alpha,
                           x = time,
                           theta))
      A[,(1 + n_z_alpha + 1)] <- (d * c(log(1 - u)) *
                                    lexc_ad(z_tau,
                                            z_alpha,
                                            x = time,
                                            theta)) /
        (riskpop + lexc_ad(z_tau,
                           z_alpha,
                           x = time,
                           theta))

      A[,(1 + n_z_alpha + 1 + 1)] <- (d * (1/tau) * (2 - alpha - beta +
                                                       (beta - 1/1 - u)) *
                                        lexc_ad(z_tau,
                                                z_alpha,
                                                x = time,
                                                theta)) /
        (riskpop + lexc_ad(z_tau,
                           z_alpha,
                           x = time,
                           theta))

    } else if (n_z_tau == 0 & n_z_alpha == 0) {


      A[,2] <- (d * c(log(1 - u)) *
                  lexc_ad(z_tau,
                          z_alpha,
                          x = time,
                          theta)) /
        (riskpop + lexc_ad(z_tau,
                           z_alpha,
                           x = time,
                           theta))

      A[, 3] <- (d * (1/tau) * (2 - alpha - beta +
                                                       (beta - 1/1 - u)) *
                                        lexc_ad(z_tau,
                                                z_alpha,
                                                x = time,
                                                theta)) /
        (riskpop + lexc_ad(z_tau,
                           z_alpha,
                           x = time,
                           theta))
      }
    solved <- solve(t(A) %*% (A), tol = 1e-30)
    return(solved)
    }
