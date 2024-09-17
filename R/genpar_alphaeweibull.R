
genpar_alphaeweibull <- function(theta, z_pcured, pophaz.alpha,
                                z_pophaz.alpha, z_pophaz.alpha_id) {
  n_z_pcured <- ncol(z_pcured)

  if (pophaz.alpha) {

    if (n_z_pcured > 0) {
      if (length(z_pophaz.alpha_id) > 0) {
        beta0 <- theta[1]
        beta_k <- theta[(2):(1 + n_z_pcured)]
        delta_k <- theta[(1 + n_z_pcured + 2 + 1):(1 + n_z_pcured + 2 + n_z_pcured)]
        lambda <- theta[(1 + n_z_pcured + 1)]
        gamma <- theta[(1 + n_z_pcured + 2)]
        thetapw <- theta[(1 + n_z_pcured + 3)]
        pophaz.alpha0 <- theta[1 + n_z_pcured + 2 + n_z_pcured + 2]
        pophaz.alpha_k <- theta[(1 + n_z_pcured + 2 + n_z_pcured + 3):(1 + n_z_pcured + 2 + n_z_pcured + 2 + length(z_pophaz.alpha_id))]

        return(c(beta0 = as.numeric(beta0),
                 beta_k = as.numeric(beta_k),
                 lambda = as.numeric(lambda),
                 gamma = as.numeric(gamma),
                 thetapw = as.numeric(thetapw),
                 delta_k = as.numeric(delta_k),
                 pophaz.alpha0 = as.numeric(pophaz.alpha0),
                 pophaz.alpha_k = as.numeric(pophaz.alpha_k)))

      }else{
        beta0 <- theta[1]
        beta_k <- theta[(2):(1 + n_z_pcured)]
        lambda <- theta[(1 + n_z_pcured + 1)]
        gamma <- theta[(1 + n_z_pcured + 2)]
        thetapw <- theta[(1 + n_z_pcured + 3)]

        delta_k <- theta[(1 + n_z_pcured + 4):(1 + n_z_pcured + 3 + n_z_pcured)]
        pophaz.alpha0 <- theta[1 + n_z_pcured + 3 + n_z_pcured + 1]

        return(c(beta0 = as.numeric(beta0),
                 beta_k = as.numeric(beta_k),
                 lambda = as.numeric(lambda),
                 gamma = as.numeric(gamma),
                 thetapw = as.numeric(thetapw),
                 delta_k = as.numeric(delta_k),
                 pophaz.alpha0 = as.numeric(pophaz.alpha0)))
        }


    } else{

      if (length(z_pophaz.alpha_id) > 0) {
        beta0 <- theta[1]
        lambda <- theta[(1 + n_z_pcured + 1)]
        gamma <- theta[(1 + n_z_pcured + 2)]
        thetapw <- theta[(1 + n_z_pcured + 3)]

        pophaz.alpha0 <- theta[1 + n_z_pcured + 3 + n_z_pcured + 1]

        return(c(beta0 = as.numeric(beta0),
                 lambda = as.numeric(lambda),
                 gamma = as.numeric(gamma),
                 thetapw = as.numeric(thetapw),
                 pophaz.alpha0 = as.numeric(pophaz.alpha0)))

      }else{
        beta0 <- theta[1]
        lambda <- theta[(1 + n_z_pcured + 1)]
        gamma <- theta[(1 + n_z_pcured + 2)]
        thetapw <- theta[(1 + n_z_pcured + 3)]

        pophaz.alpha0 <- theta[1 + n_z_pcured + 3 + n_z_pcured + 1]

        return(c(beta0 = as.numeric(beta0),
                 lambda = as.numeric(lambda),
                 gamma = as.numeric(gamma),
                 thetapw = as.numeric(thetapw),
                 pophaz.alpha0 = as.numeric(pophaz.alpha0)))
      }

      }

    }else {

      if (n_z_pcured > 0) {
          beta0 <- theta[1]
          beta_k <- theta[(2):(1 + n_z_pcured)]
          delta_k <- theta[(1 + n_z_pcured + 2 + 1):(1 + n_z_pcured + 2 + n_z_pcured)]
          lambda <- theta[(1 + n_z_pcured + 1)]
          gamma <- theta[(1 + n_z_pcured + 2)]
          thetapw <- theta[(1 + n_z_pcured + 3)]

          return(c(beta0 = as.numeric(beta0),
                   beta_k = as.numeric(beta_k),
                   lambda = as.numeric(lambda),
                   gamma = as.numeric(gamma),
                   thetapw = as.numeric(thetapw),
                   delta_k = as.numeric(delta_k)))
      } else{
        beta0 <- theta[1]
        lambda <- theta[(1 + n_z_pcured + 1)]
        gamma <- theta[(1 + n_z_pcured + 2)]
        thetapw <- theta[(1 + n_z_pcured + 3)]

        return(c(beta0 = as.numeric(beta0),
                 lambda = as.numeric(lambda),
                 gamma = as.numeric(gamma),
                 thetapw = as.numeric(thetapw)))
      }

    }
  }




