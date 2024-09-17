ttc_tneh2 <- function(object = object,
                      fx = fx,
                      xmax = xmax,
                      z_tau,
                      z_alpha, epsilon = epsilon) {
    if(object$link_tau=="linear"){
      fx <- function(z_alpha, z_tau, x, object, epsilon = epsilon) {
        theta <- object$coefficients
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


          alpha <- (alpha0 + z_alpha %*% alpha_k)
          tau <- (tau0 + z_tau %*% tau_z)
          u <- x / (tau)
          beta2 <- beta
          cumhaz <- ifelse((x <= (tau)),
                           (tau) * beta(alpha, beta2) * stats::pbeta(u, alpha, beta2),
                           (tau) * beta(alpha, beta2) * stats::pbeta(1, alpha, beta2)
          )
          cumhaz2 <- (tau) * beta(alpha, beta2) * stats::pbeta(1, alpha, beta2)

        } else if (n_z_tau == 0 & n_z_alpha > 0) {

          alpha_k <- theta[2:(n_z_alpha + 1)]
          beta <- (theta[n_z_alpha + 2])
          tau0 <- theta[n_z_alpha + 2 + 1]
          alpha <- (alpha0 + z_alpha %*% alpha_k)
          tau <- (tau0 )

          u <- x / (tau)

          beta2 <- beta
          cumhaz <- ifelse((x <= (tau)),
                           (tau) * beta(alpha, beta2) * stats::pbeta(u, alpha, beta2),
                           (tau) * beta(alpha, beta2) * stats::pbeta(1, alpha, beta2)
          )
          cumhaz2 <- (tau) * beta(alpha, beta2) * stats::pbeta(1, alpha, beta2)

        }else if (n_z_tau > 0 & n_z_alpha == 0) {

          beta <- (theta[2])
          tau0 <- theta[n_z_alpha + 2 + 1]
          tau_z <- theta[(n_z_alpha + 2 + 1 + 1):(n_z_alpha + 2 + n_z_tau + 1)]

          alpha <- (alpha0)

          tau <- (tau0 + z_tau %*% tau_z)

          u <- x / (tau)
          u <- u*(u < 1)

          beta2 <- beta
          cumhaz <- ifelse((x <= (tau)),
                           (tau) * beta(alpha, beta2) * stats::pbeta(u, alpha, beta2),
                           (tau) * beta(alpha, beta2) * stats::pbeta(1, alpha, beta2)
          )
          cumhaz2 <- (tau) * beta(alpha, beta2) * stats::pbeta(1, alpha, beta2)

        }
        else if (n_z_tau == 0 & n_z_alpha == 0) {
          beta <- (theta[2])
          tau0 <- theta[3]
          alpha <- alpha0
          tau <- (tau0)

          u <- x / (tau)
          u <- u*(u < 1)

          beta2 <- beta
          cumhaz <- ifelse((x <= (tau)),
                           (tau) * beta(alpha, beta2) * stats::pbeta(u, alpha, beta2),
                           (tau) * beta(alpha, beta2) * stats::pbeta(1, alpha, beta2)
          )
          cumhaz2 <- (tau) * beta(alpha, beta2) * stats::pbeta(1, alpha, beta2)
        }

        pt <- (exp(-cumhaz2) / exp(-cumhaz))

        fx <-  ifelse(exp(-cumhaz2) >= c(1 - epsilon), 0, c(pt - (c(1 - epsilon))))
        return(fx)
      }



      n_z_tau <- ncol(z_tau)
      n_z_alpha <- ncol(z_alpha)
      n_z_tau_ad <- n_z_tau - 1
      n_z_alpha_ad <- n_z_alpha - 1


      if (n_z_tau > 0 & n_z_alpha > 0) {
        res_root <- lapply(1:nrow(z_alpha), function(i) {
          z_alpha_i <- (matrix(z_alpha[i,], ncol = n_z_alpha))
          colnames(z_alpha_i) <- colnames(z_alpha)
          z_tau_i <- (matrix(z_tau[i,], ncol = n_z_tau))
          colnames(z_tau_i) <- colnames(z_tau)
          uniroot(fx,
                  c(0, xmax),
                  z_alpha = z_alpha_i,
                  z_tau = z_tau_i,
                  object = object, epsilon = epsilon)
        })
      } else if (n_z_tau == 0 & n_z_alpha > 0) {
        res_root <- lapply(1:nrow(z_alpha), function(i) {
          z_alpha_i <- (matrix(z_alpha[i,], ncol = n_z_alpha))
          colnames(z_alpha_i) <- colnames(z_alpha)
          uniroot(fx,
                  c(0, xmax),
                  z_alpha = z_alpha_i,
                  z_tau = z_tau,
                  object = object, epsilon = epsilon)
        })


      }else if (n_z_tau > 0 & n_z_alpha == 0) {
        res_root <- lapply(1:nrow(z_tau), function(i) {
          z_tau_i <- (matrix(z_tau[i,], ncol = n_z_tau))
          colnames(z_tau_i) <- colnames(z_tau)
          uniroot(fx,
                  c(0, xmax),
                  z_alpha = z_alpha,
                  z_tau = z_tau_i,
                  object = object, epsilon = epsilon)
        })

      }
      else if (n_z_tau == 0 & n_z_alpha == 0) {
        res_root <- lapply(1:nrow(z_tau), function(i) {
          z_tau_i <- (matrix(z_tau[i,], ncol = n_z_tau))
          colnames(z_tau_i) <- colnames(z_tau)
          stats::uniroot(fx,
                         c(0, xmax),
                         z_alpha = z_alpha,
                         z_tau = z_tau_i,
                         object = object, epsilon = epsilon)
        })

      }

      return(res_root)
    }else if(object$link_tau=="loglinear"){
      fx <- function(z_alpha, z_tau, x, object, epsilon = epsilon) {
      theta <- object$coefficients
      n_z_tau <- ncol(z_tau)
      n_z_alpha <- ncol(z_alpha)
      n_z_tau_ad <- n_z_tau - 1
      n_z_alpha_ad <- n_z_alpha - 1
      alpha0 <- (theta[1])
      if (n_z_tau > 0 & n_z_alpha > 0) {
        alpha_k <- theta[2:(n_z_alpha + 1)]
        beta <- exp(theta[n_z_alpha + 2])+1
        tau0 <- theta[n_z_alpha + 2 + 1]
        tau_z <- theta[(n_z_alpha + 2 + 1 + 1):(n_z_alpha + 2 + n_z_tau + 1)]


        alpha <- exp(alpha0 + z_alpha %*% alpha_k)
        tau <- exp(tau0 + z_tau %*% tau_z)
        u <- x / (tau)
        beta2 <- beta
        cumhaz <- ifelse((x <= (tau)),
                         (tau) * beta(alpha, beta2) * stats::pbeta(u, alpha, beta2),
                         (tau) * beta(alpha, beta2) * stats::pbeta(1, alpha, beta2)
        )
        cumhaz2 <- (tau) * beta(alpha, beta2) * stats::pbeta(1, alpha, beta2)

      } else if (n_z_tau == 0 & n_z_alpha > 0) {

        alpha_k <- theta[2:(n_z_alpha + 1)]
        beta <- exp(theta[n_z_alpha + 2])+1
        tau0 <- theta[n_z_alpha + 2 + 1]
        alpha <- exp(alpha0 + z_alpha %*% alpha_k)
        tau <- exp(tau0 )

        u <- x / (tau)

        beta2 <- beta
        cumhaz <- ifelse((x <= (tau)),
                         (tau) * beta(alpha, beta2) * stats::pbeta(u, alpha, beta2),
                         (tau) * beta(alpha, beta2) * stats::pbeta(1, alpha, beta2)
        )
        cumhaz2 <- (tau) * beta(alpha, beta2) * stats::pbeta(1, alpha, beta2)

      }else if (n_z_tau > 0 & n_z_alpha == 0) {

        beta <- exp(theta[2])+1
        tau0 <- theta[n_z_alpha + 2 + 1]
        tau_z <- theta[(n_z_alpha + 2 + 1 + 1):(n_z_alpha + 2 + n_z_tau + 1)]

        alpha <- exp(alpha0)

        tau <- exp(tau0 + z_tau %*% tau_z)

        u <- x / (tau)
        u <- u*(u < 1)

        beta2 <- beta
        cumhaz <- ifelse((x <= (tau)),
                         (tau) * beta(alpha, beta2) * stats::pbeta(u, alpha, beta2),
                         (tau) * beta(alpha, beta2) * stats::pbeta(1, alpha, beta2)
        )
        cumhaz2 <- (tau) * beta(alpha, beta2) * stats::pbeta(1, alpha, beta2)

      }
      else if (n_z_tau == 0 & n_z_alpha == 0) {
        beta <- exp(theta[2])+1
        tau0 <- theta[3]
        alpha <- exp(alpha0)
        tau <- exp(tau0)

        u <- x / (tau)
        u <- u*(u < 1)

        beta2 <- beta
        cumhaz <- ifelse((x <= (tau)),
                         (tau) * beta(alpha, beta2) * stats::pbeta(u, alpha, beta2),
                         (tau) * beta(alpha, beta2) * stats::pbeta(1, alpha, beta2)
        )
        cumhaz2 <- (tau) * beta(alpha, beta2) * stats::pbeta(1, alpha, beta2)
      }

      pt <- (exp(-cumhaz2) / exp(-cumhaz))

      fx <-  ifelse(exp(-cumhaz2) >= c(1 - epsilon), 0, c(pt - (c(1 - epsilon))))
      return(fx)
    }



    n_z_tau <- ncol(z_tau)
    n_z_alpha <- ncol(z_alpha)
    n_z_tau_ad <- n_z_tau - 1
    n_z_alpha_ad <- n_z_alpha - 1


    if (n_z_tau > 0 & n_z_alpha > 0) {
      res_root <- lapply(1:nrow(z_alpha), function(i) {
        z_alpha_i <- (matrix(z_alpha[i,], ncol = n_z_alpha))
        colnames(z_alpha_i) <- colnames(z_alpha)
        z_tau_i <- (matrix(z_tau[i,], ncol = n_z_tau))
        colnames(z_tau_i) <- colnames(z_tau)
        uniroot(fx,
                c(0, xmax),
                z_alpha = z_alpha_i,
                z_tau = z_tau_i,
                object = object, epsilon = epsilon)
      })
    } else if (n_z_tau == 0 & n_z_alpha > 0) {
      res_root <- lapply(1:nrow(z_alpha), function(i) {
        z_alpha_i <- (matrix(z_alpha[i,], ncol = n_z_alpha))
        colnames(z_alpha_i) <- colnames(z_alpha)
        uniroot(fx,
                c(0, xmax),
                z_alpha = z_alpha_i,
                z_tau = z_tau,
                object = object, epsilon = epsilon)
      })


    }else if (n_z_tau > 0 & n_z_alpha == 0) {
      res_root <- lapply(1:nrow(z_tau), function(i) {
        z_tau_i <- (matrix(z_tau[i,], ncol = n_z_tau))
        colnames(z_tau_i) <- colnames(z_tau)
        uniroot(fx,
                c(0, xmax),
                z_alpha = z_alpha,
                z_tau = z_tau_i,
                object = object, epsilon = epsilon)
      })

    }
    else if (n_z_tau == 0 & n_z_alpha == 0) {
      res_root <- lapply(1:nrow(z_tau), function(i) {
        z_tau_i <- (matrix(z_tau[i,], ncol = n_z_tau))
        colnames(z_tau_i) <- colnames(z_tau)
        stats::uniroot(fx,
                       c(0, xmax),
                       z_alpha = z_alpha,
                       z_tau = z_tau_i,
                       object = object, epsilon = epsilon)
      })

    }

    return(res_root)

    }

}
