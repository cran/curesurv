
genpar_alphaweibull <- function(theta, z_pcured, hazpop.alpha,
                                z_hazpop.alpha, z_hazpop.alpha_id) {
  n_z_pcured <- ncol(z_pcured)

  if (hazpop.alpha) {

    if (n_z_pcured > 0) {
      if (length(z_hazpop.alpha_id) > 0) {
        beta0 <- theta[1]
        beta_k <- theta[(2):(1 + n_z_pcured)]
        delta_k <- theta[(1 + n_z_pcured + 2 + 1):(1 + n_z_pcured + 2 + n_z_pcured)]
        lambda <- theta[(1 + n_z_pcured + 1)]
        gamma <- theta[(1 + n_z_pcured + 2)]
        hazpop.alpha0 <- theta[1 + n_z_pcured + 2 + n_z_pcured + 1]
        hazpop.alpha_k <- theta[(1 + n_z_pcured + 2 + n_z_pcured + 2):(1 + n_z_pcured + 2 + n_z_pcured + 1 + length(z_hazpop.alpha_id))]

        return(c(beta0 = as.numeric(beta0),
                 beta_k = as.numeric(beta_k),
                 lambda = as.numeric(lambda),
                 gamma = as.numeric(gamma),
                 delta_k = as.numeric(delta_k),
                 hazpop.alpha0 = as.numeric(hazpop.alpha0),
                 hazpop.alpha_k = as.numeric(hazpop.alpha_k)))

      }else{
        beta0 <- theta[1]
        beta_k <- theta[(2):(1 + n_z_pcured)]
        delta_k <- theta[(1 + n_z_pcured + 2 + 1):(1 + n_z_pcured + 2 + n_z_pcured)]
        lambda <- theta[(1 + n_z_pcured + 1)]
        gamma <- theta[(1 + n_z_pcured + 2)]
        hazpop.alpha0 <- theta[1 + n_z_pcured + 2 + n_z_pcured + 1]

        return(c(beta0 = as.numeric(beta0),
                 beta_k = as.numeric(beta_k),
                 lambda = as.numeric(lambda),
                 gamma = as.numeric(gamma),
                 delta_k = as.numeric(delta_k),
                 hazpop.alpha0 = as.numeric(hazpop.alpha0)))
        }


    } else{

      if (length(z_hazpop.alpha_id) > 0) {
        beta0 <- theta[1]
        lambda <- theta[(1 + n_z_pcured + 1)]
        gamma <- theta[(1 + n_z_pcured + 2)]
        hazpop.alpha0 <- theta[1 + n_z_pcured + 2 + n_z_pcured + 1]

        return(c(beta0 = as.numeric(beta0),
                 lambda = as.numeric(lambda),
                 gamma = as.numeric(gamma),
                 hazpop.alpha0 = as.numeric(hazpop.alpha0)))

      }else{
        beta0 <- theta[1]
        lambda <- theta[(1 + n_z_pcured + 1)]
        gamma <- theta[(1 + n_z_pcured + 2)]
        hazpop.alpha0 <- theta[1 + n_z_pcured + 2 + n_z_pcured + 1]

        return(c(beta0 = as.numeric(beta0),
                 lambda = as.numeric(lambda),
                 gamma = as.numeric(gamma),
                 hazpop.alpha0 = as.numeric(hazpop.alpha0)))
      }

      }

    }else {

      if (n_z_pcured > 0) {
          beta0 <- theta[1]
          beta_k <- theta[(2):(1 + n_z_pcured)]
          delta_k <- theta[(1 + n_z_pcured + 2 + 1):(1 + n_z_pcured + 2 + n_z_pcured)]
          lambda <- theta[(1 + n_z_pcured + 1)]
          gamma <- theta[(1 + n_z_pcured + 2)]

          return(c(beta0 = as.numeric(beta0),
                   beta_k = as.numeric(beta_k),
                   lambda = as.numeric(lambda),
                   gamma = as.numeric(gamma),
                   delta_k = as.numeric(delta_k)))
      } else{
        beta0 <- theta[1]
        lambda <- theta[(1 + n_z_pcured + 1)]
        gamma <- theta[(1 + n_z_pcured + 2)]

        return(c(beta0 = as.numeric(beta0),
                 lambda = as.numeric(lambda),
                 gamma = as.numeric(gamma)))
      }

    }
  }




