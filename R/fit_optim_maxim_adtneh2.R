fit.opt.maxim.adtneh2<- function(x = x,
                                      d = d,
                                      z_tau = z_tau,
                                      z_alpha = z_alpha,
                                      theta_init = theta_init,
                                      theta_lower = theta_lower,
                                      theta_upper = theta_upper,
                                      pophaz = pophaz,
                                      cumpophaz = cumpophaz,
                                      method_opt = method_opt,
                                      pophaz.alpha = pophaz.alpha,
                                      maxit_opt = maxit_opt,
                                      gradient = gradient,
                                      hessian_varcov = hessian_varcov,
                                      ncoor_des = ncoor_des,
                                      optim_func = optim_func,
                                      iter_eps = iter_eps,
                                      trace = trace,
                                      optim_fixed = optim_fixed,
                                      optimizer = optimizer) {

  fixed <- optim_fixed
  theta <- theta2 <- theta_init
  n_z_tau <- ncol(z_tau)
  n_z_alpha <- ncol(z_alpha)
  n_z_tau_ad <- n_z_tau - 1
  n_z_alpha_ad <- n_z_alpha - 1

  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))

  if (pophaz.alpha) {
    if (n_z_tau >= 1 & n_z_alpha >= 1) {
      names(theta) <- c('alpha0',
                        paste('alpha', 1:n_z_alpha, colnames(z_alpha) ,sep = "_"),
                        'beta',
                        'tau0',
                        paste('tau', 1:n_z_tau, colnames(z_tau), sep = "_"),
                        'pophaz.alpha')
      names(theta2) <- c('alpha0',
                         paste('alpha', 1:n_z_alpha,sep = "_"),
                         'beta',
                         'tau0',
                         paste('tau', 1:n_z_tau, sep = "_"),
                         'alpha')


    }
    else if (n_z_tau >= 1 & n_z_alpha == 0) {
      names(theta) <- c('alpha0',
                        'beta',
                        'tau0',
                        paste('tau', 1:n_z_tau, colnames(z_tau), sep = "_"),
                        'pophaz.alpha')
      names(theta2) <- c('alpha0',
                         'beta',
                         'tau0',
                         paste('tau', 1:n_z_tau, sep = "_"),
                         'alpha')


    } else if (n_z_tau == 0 & n_z_alpha >= 1) {
      names(theta) <- c('alpha0',
                        paste('alpha', 1:n_z_alpha, colnames(z_alpha) ,sep = "_"),
                        'beta',
                        'tau0',
                        'pophaz.alpha')
      names(theta2) <- c('alpha0',
                         paste('alpha', 1:n_z_alpha ,sep = "_"),
                         'beta',
                         'tau0',
                         'alpha')



    }
    else if (n_z_tau == 0 & n_z_alpha == 0) {
      names(theta) <- c('alpha0','beta','tau0', 'pophaz.alpha')
      names(theta2) <- c('alpha0','beta','tau0', 'alpha')

    }
  }else{
    if (n_z_tau >= 1 & n_z_alpha >= 1) {
      names(theta) <- c('alpha0',
                        paste('alpha', 1:n_z_alpha, colnames(z_alpha) ,sep = "_"),
                        'beta',
                        'tau0',
                        paste('tau', 1:n_z_tau, colnames(z_tau), sep = "_"))
      names(theta2) <- c('alpha0',
                         paste('alpha', 1:n_z_alpha ,sep = "_"),
                         'beta',
                         'tau0',
                         paste('tau', 1:n_z_tau, sep = "_"))


    }else  if (n_z_tau == 0 & n_z_alpha >= 1) {
      names(theta) <- c('alpha0',
                        paste('alpha', 1:n_z_alpha, colnames(z_alpha) ,sep = "_"),
                        'beta',
                        'tau0')
      names(theta2) <- c('alpha0',
                         paste('alpha', 1:n_z_alpha,sep = "_"),
                         'beta',
                         'tau0')


    } else if (n_z_tau >= 1 & n_z_alpha == 0) {
      names(theta) <- c('alpha0',
                        'beta',
                        'tau0',
                        paste('tau', 1:n_z_tau, colnames(z_tau), sep = "_"))
      names(theta2) <- c('alpha0',
                         'beta',
                         'tau0',
                         paste('tau', 1:n_z_tau, sep = "_"))


    } else if (n_z_tau == 0 & n_z_alpha == 0) {
      names(theta) <- c('alpha0','beta','tau0')
      names(theta2) <- c('alpha0','beta','tau0')

    }
  }






  if (gradient == TRUE) {
    if (n_z_alpha > 0 & n_z_tau > 0) {
      adtneh_ll_fn_deriv_i <- Deriv::Deriv( ~ -sum(
        d * log(c(exp(alpha_pophaz)) * pophaz +
                  ifelse((x <= (tau)),
                         ((x / (tau))^ ((alpha) - 1) * (1 - (x / (tau))) ^ ((beta - 1))),
                         0)) -
          ifelse((x <= (tau)),
                 (tau) * beta(alpha, beta) * stats::pbeta(u, alpha, beta),
                 (tau) * beta(alpha, beta) * stats::pbeta(1, alpha, beta) -
                   cumpophaz * c(exp(alpha_pophaz)), na.rm = TRUE)),
        c("alpha0","alpha_k","beta","tau0","tau_z"),
        combine = "cbind",
        function(alpha0, alpha_k, beta, tau0, tau_z, x, d, pophaz, cumpophaz,
                 alpha_pophaz, z_alpha, z_tau, n_z_alpha, n_z_tau)NULL)


      adtneh_ll_fn_deriv_all <- function(theta
      ) {
        alpha_k <- theta[2:(n_z_alpha + 1)]
        beta0 <- theta[1]

        beta <- (theta[n_z_alpha + 2])
        tau0 <- theta[n_z_alpha + 2 + 1]
        tau_z <- theta[(n_z_alpha + 2 + 1 + 1):(n_z_alpha + 2 + n_z_tau + 1)]
        alpha <- (alpha0 + z_alpha %*% alpha_k)
        tau <- (tau0 + z_tau %*% tau_z)
        u <- x / (tau)



        alpha_pophaz <- c(matrix(1, ncol = 1, nrow = length(pophaz)) %*% 0)
        pophaz <- c(pophaz)
        cumpophaz <- c(cumpophaz)
        adtneh_ll_fn_deriv_i(alpha0, alpha_k, beta, tau0, x, d,
                             pophaz, cumpophaz,
                             alpha_pophaz, z_alpha,
                             n_z_alpha, n_z_tau, out)
        gradval <- colSums(adtneh_ll_fn_deriv_i(alpha0,
                                                alpha_k,
                                                beta,
                                                tau0,
                                                x,
                                                d,
                                                pophaz,
                                                cumpophaz,
                                                alpha_pophaz,
                                                z_alpha,
                                                n_z_alpha,
                                                n_z_tau))


        names(gradval) <- NULL
        return(gradval)
      }




    }else if (n_z_alpha == 0 & n_z_tau > 0) {
      #is it normal that there is no code here ?

    } else if (n_z_alpha > 0 & n_z_tau == 0) {



      out <- statmod::gauss.quad(12, "legendre")

      adtneh_ll_fn_deriv_i1 <-
        function(alpha0,
                 alpha_k=NULL,
                 beta,
                 tau0,
                 x,
                 d,
                 pophaz,
                 cumpophaz,
                 alpha_pophaz,
                 z_alpha=NULL,
                 n_z_alpha,
                 n_z_tau,
                 out) {
          .e1 <- x/2
          .e2 <- x/tau0
          .e3 <- c(alpha0 + z_alpha %*% t(alpha_k))
          .e4 <- .e3 - 1
          .e5 <- beta - 1
          .e7 <- outer(.e1, out$nodes, "^") + .e1
          .e8 <- .e7/tau0
          .e9 <- 1 - .e2
          .e10 <- .e9^.e5
          .e11 <- .e2^.e4
          .e12 <- 1 - .e8
          .e13 <- .e8^.e4
          .e14 <- .e12^.e5
          .e16 <- .e10 * .e11 + pophaz * exp(alpha_pophaz)
          .e17 <- log(tau0)
          .e18 <- .e13 * .e14
          .e19 <- .e3 - 2
          .e20 <- beta - 2
          .e21 <- d * .e10
          .e23 <- log(.e7) - .e17
          .e25 <- log(x) - .e17
          .e26 <- tau0^2
          cbind(alpha0 = -(.e21 * .e25 * .e11/.e16 - x *
                             rowSums(x = .e18 * .e23 * out$weights)/2),
                alpha_k = -(d * z_alpha * .e10 * .e25 * .e11/.e16 - x *
                              rowSums(x = z_alpha * .e13 * .e14 *
                                        .e23 * out$weights)/2),
                beta = -(.e21 * log(.e9) * .e11/.e16 -
                           x * rowSums(x = .e18 * log(.e12) *
                                         out$weights)/2),
                tau0 = -(x * (d * (.e9^.e20 * .e5 * .e11 -
                                     .e10 * .e4 * .e2^.e19)/(.e26 * .e16) -
                                rowSums(x = (.e13 * .e12^.e20 * .e5 -
                                               .e8^.e19 * .e14 * .e4) *
                                          out$weights * .e7/.e26)/2)))
        }

      adtneh_ll_fn_deriv_i2 <-
        function(alpha0,
                 alpha_k=NULL,
                 beta,
                 tau0,
                 x,
                 d,
                 pophaz,
                 cumpophaz,
                 alpha_pophaz,
                 z_alpha=NULL,
                 n_z_alpha,
                 n_z_tau) {
          .e2 <- c(z_alpha %*% t(alpha_k))
          .e3 <- alpha0 + .e2
          if (any(.e3 < 0)) {
            stop("change initial values for alpha_k or the scale of z_alpha")
          }
          .e4 <- beta(.e3, beta)
          .e7 <- digamma(alpha0 + beta + .e2)
          .e8 <- tau0 * .e4
          .e11 <- .e8 * (digamma(.e3) - .e7)
          cbind(alpha0 = .e11,
                alpha_k = .e11 * z_alpha,
                beta = .e8 * (digamma(beta) - .e7),
                tau0 = .e4)

        }

      adtneh_ll_fn_deriv_i <-
        function(alpha0,
                 alpha_k=NULL,
                 beta,
                 tau0,
                 x,
                 d,
                 pophaz,
                 cumpophaz,
                 alpha_pophaz,
                 z_alpha=NULL,
                 n_z_alpha,
                 n_z_tau,
                 out){
          gradi <- rbind(
            adtneh_ll_fn_deriv_i1(
              alpha0,
              alpha_k,
              beta,
              tau0,
              x[which((x <= tau0))],
              d[which((x <= tau0))],
              pophaz[which((x <= tau0))],
              cumpophaz[which((x <= tau0))],
              alpha_pophaz[which((x <= tau0))],
              z_alpha[which((x <= tau0)),] ,
              n_z_alpha,
              n_z_tau,
              out
            ),
            adtneh_ll_fn_deriv_i2(
              alpha0,
              alpha_k,
              beta,
              tau0,
              x[which((x > tau0))],
              d[which((x > tau0))],
              pophaz[which((x > tau0))],
              cumpophaz[which((x > tau0))],
              alpha_pophaz[which((x > tau0))],
              z_alpha[which((x > tau0)),],
              n_z_alpha,
              n_z_tau
            ))


          return(gradi)

        }




      adtneh_ll_fn_deriv_all <- function(theta
      ) {
        alpha0 <- theta[1]
        alpha_k <- theta[2:(n_z_alpha + 1)]

        beta <- (theta[n_z_alpha + 2])
        tau0 <- theta[n_z_alpha + 2 + 1]
        alpha <- (alpha0 + z_alpha %*% alpha_k)
        tau <- (tau0)
        u <- x / (tau)
        alpha_pophaz <- c(matrix(1, ncol = 1, nrow = length(pophaz)) %*% 0)
        pophaz <- c(pophaz)
        cumpophaz <- c(cumpophaz)

        gradval1 <- adtneh_ll_fn_deriv_i(alpha0, alpha_k, beta, tau0, x, d,
                                         pophaz, cumpophaz,
                                         alpha_pophaz, z_alpha,
                                         n_z_alpha, n_z_tau, out)


        gradval <- colSums(gradval1, na.rm = TRUE)
        return(gradval)
      }








    } else if (n_z_alpha == 0 & n_z_tau == 0) {


      out <- statmod::gauss.quad(12, "legendre")

      adtneh_ll_fn_deriv_i1 <-
          function(alpha0,
                   alpha_k=NULL,
                   beta,
                   tau0,
                   x,
                   d,
                   pophaz,
                   cumpophaz,
                   alpha_pophaz,
                   z_alpha=NULL,
                   n_z_alpha,
                   n_z_tau,
                   out){
          .e1 <- x / 2
          .e2 <- x / tau0
          .e3 <- alpha0 - 1
          .e4 <- beta - 1
          .e6 <- outer(.e1, out$nodes, "^") + .e1
          .e7 <- .e6 / tau0
          .e8 <- 1 - .e2
          .e9 <- .e8 ^ .e4
          .e10 <- .e2 ^ .e3
          .e11 <- 1 - .e7
          .e12 <- .e7 ^ .e3
          .e13 <- .e11 ^ .e4
          .e15 <- .e9 * .e10 + pophaz * exp(alpha_pophaz)
          .e16 <- .e12 * .e13
          .e17 <- alpha0 - 2
          .e18 <- beta - 2
          .e19 <- d * .e9
          .e20 <- log(tau0)
          .e21 <- tau0 ^ 2
          cbind(
            alpha0 = -(
              .e19 * (log(x) - .e20) * .e10 / .e15 - x *
                rowSums(x = .e16 * (log(.e6) - .e20) *
                          out$weights) / 2),
            beta = -(.e19 * log(.e8) * .e10 / .e15 -
                       x * rowSums(x = .e16 * log(.e11) * out$weights) / 2),
            tau0 = -(x * (d * (.e8 ^ .e18 * .e4 * .e10 -
                                 .e9 * .e3 * .e2 ^ .e17) / (.e21 * .e15) -
                            rowSums(x = (.e12 * .e11 ^ .e18 * .e4 -
                                           .e7 ^ .e17 * .e13 * .e3) *
                                      out$weights * .e6 / .e21) / 2
            ))
          )
        }

      adtneh_ll_fn_deriv_i2 <-
        function(alpha0,
                 alpha_k=NULL,
                 beta,
                 tau0,
                 x,
                 d,
                 pophaz,
                 cumpophaz,
                 alpha_pophaz,
                 z_alpha=NULL,
                 n_z_alpha,
                 n_z_tau){
          .e1 <- beta(alpha0, beta)
          .e3 <- digamma(alpha0 + beta)
          .e4 <- tau0 * .e1
          cbind(alpha0 = .e4 * (digamma(alpha0) - .e3),
                beta = .e4 *(digamma(beta) - .e3),
                tau0 = .e1)

        }

      adtneh_ll_fn_deriv_i <-
        function(alpha0,
                 alpha_k=NULL,
                 beta,
                 tau0,
                 x,
                 d,
                 pophaz,
                 cumpophaz,
                 alpha_pophaz,
                 z_alpha=NULL,
                 n_z_alpha,
                 n_z_tau,
                 out){
          gradi <- rbind(
            adtneh_ll_fn_deriv_i1(
              alpha0,
              beta,
              tau0,
              x[which((x <= tau0))],
              d[which((x <= tau0))],
              pophaz[which((x <= tau0))],
              cumpophaz[which((x <= tau0))],
              alpha_pophaz[which((x <= tau0))],
              n_z_alpha,
              n_z_tau,
              out
            ),
            t(sapply(1:length(which((x > tau0))),
                     function(i)adtneh_ll_fn_deriv_i2(
                       alpha0,
                       beta,
                       tau0,
                       x[which((x > tau0))],
                       d[which((x > tau0))],
                       pophaz[which((x > tau0))],
                       cumpophaz[which((x > tau0))],
                       alpha_pophaz[which((x > tau0))],
                       n_z_alpha,
                       n_z_tau
                     )))

          )
          return(gradi)

        }




      adtneh_ll_fn_deriv_all <- function(theta
      ) {
        alpha0 <- theta[1]
        beta <- (theta[2])
        tau0 <- theta[3]
        alpha <- (alpha0)
        tau <- (tau0)
        u <- x / (tau)
        alpha_pophaz <- c(matrix(1, ncol = 1, nrow = length(pophaz)) %*% 0)
        pophaz <- c(pophaz)
        cumpophaz <- c(cumpophaz)

        gradval1 <- adtneh_ll_fn_deriv_i(alpha0, beta, tau0, x, d,
                                         pophaz, cumpophaz,
                                         alpha_pophaz,
                                         n_z_alpha, n_z_tau, out)


        gradval <- colSums(gradval1, na.rm = TRUE)
        return(gradval)
      }






    }

  }

  adtneh_ll_fn <- function(alpha0 = NULL,
                           alpha_k = NULL,
                           beta = NULL,
                           tau0 = NULL,
                           tau_z = NULL,
                           alpha = NULL) {

    if (pophaz.alpha) {
      if (n_z_alpha == 0 & n_z_tau == 0) {
        theta <- c(alpha0, beta, tau0, alpha)
      } else if (n_z_alpha > 0 & n_z_tau == 0) {
        theta <- c(alpha0, alpha_k, beta, tau0, alpha)
      }else if (n_z_alpha == 0 & n_z_tau > 0) {
        theta <- c(alpha0, beta, tau0, tau_z, alpha)
      }else if (n_z_alpha > 0 & n_z_tau > 0) {
        theta <- c(alpha0, alpha_k, beta, tau0, tau_z, alpha)
      }
    }else{
      if (n_z_alpha == 0 & n_z_tau == 0) {
        theta <- c(alpha0, beta, tau0)
      } else if (n_z_alpha > 0 & n_z_tau == 0) {
        theta <- c(alpha0, alpha_k, beta, tau0)
      }else if (n_z_alpha == 0 & n_z_tau > 0) {
        theta <- c(alpha0, beta, tau0, tau_z)
      }else if (n_z_alpha > 0 & n_z_tau > 0) {
        theta <- c(alpha0, alpha_k, beta, tau0, tau_z)
      }
    }


    if (anyNA(theta))
      return(1e10)
    if (pophaz.alpha) {
      alpha_pophaz <- matrix(1,
                             ncol = 1,
                             nrow = length(pophaz)) %*% theta[n_z_alpha + 2 + n_z_tau + 1 + 1]
    }else{

      alpha_pophaz <- matrix(1, ncol = 1,
                             nrow = length(pophaz)) %*% 0
    }

    n_z_tau <- ncol(z_tau)
    n_z_alpha <- ncol(z_alpha)
    n_z_tau_ad <- n_z_tau - 1
    n_z_alpha_ad <- n_z_alpha - 1
    alpha0 <- theta[1]

    if (n_z_tau > 0 & n_z_alpha > 0) {
      alpha_k <- theta[2:(n_z_alpha + 1)]
      beta <- (theta[n_z_alpha + 2])
      tau0 <- theta[n_z_alpha + 2 + 1]
      tau_z <- theta[(n_z_alpha + 2 + 1 + 1):(n_z_alpha + 2 + n_z_tau + 1)]

      if (any(alpha0 + z_alpha %*% alpha_k <= 0) | any(tau0 + z_tau %*% tau_z <= 0)) {
        sum_minusloklik <- 1e10
      }else {
        sum_minusloklik <- -sum(d * log(c(exp(alpha_pophaz)) * pophaz +
                                          lexc_ad2(z_tau, z_alpha, x, theta)) -
                                  cumLexc_ad2(z_tau, z_alpha, x, theta) -
                                  cumpophaz * c(exp(alpha_pophaz)), na.rm = TRUE)
      }
    }else if (n_z_tau > 0 & n_z_alpha == 0) {
      beta <- (theta[2])
      tau0 <- theta[2 + 1]
      tau_z <- theta[(2 + 1 + 1):(2 + n_z_tau + 1)]

      if (any(alpha0  <= 0) | any(tau0 + z_tau %*% tau_z <= 0)) {
        sum_minusloklik <- 1e10
      }else {
        sum_minusloklik <- -sum(d * log(c(exp(alpha_pophaz)) * pophaz +
                                          lexc_ad2(z_tau, z_alpha, x, theta)) -
                                  cumLexc_ad2(z_tau, z_alpha, x, theta) -
                                  cumpophaz * c(exp(alpha_pophaz)), na.rm = TRUE)
      }
    }else if (n_z_tau == 0 & n_z_alpha > 0) {
      alpha_k <- theta[2:(n_z_alpha + 1)]
      beta <- (theta[n_z_alpha + 2])
      tau0 <- theta[n_z_alpha + 2 + 1]
      if (any(alpha0 + z_alpha %*% alpha_k <= 0) | tau0  <= 0) {
        sum_minusloklik <- 1e10
      }else {
        sum_minusloklik <- -sum(d * log(c(exp(alpha_pophaz)) * pophaz +
                                          lexc_ad2(z_tau, z_alpha, x, theta)) -
                                  cumLexc_ad2(z_tau, z_alpha, x, theta) -
                                  cumpophaz * c(exp(alpha_pophaz)), na.rm = TRUE)
      }
    }else if (n_z_tau == 0 & n_z_alpha == 0) {
      tau0 <- theta[n_z_alpha + 2 + 1]

      if (alpha0  <= 0 | tau0 <= 0) {
        sum_minusloklik <- 1e10
      }else {
        sum_minusloklik <- -sum(d * log(c(exp(alpha_pophaz)) * pophaz +
                                          lexc_ad2(z_tau, z_alpha, x, theta)) -
                                  cumLexc_ad2(z_tau, z_alpha, x, theta) -
                                  cumpophaz * c(exp(alpha_pophaz)), na.rm = TRUE)
      }
    }




    if (anyNA(sum_minusloklik)) {
      return(1e10)
    } else {
      return(sum_minusloklik)
    }

  }








  adtneh_ll_fn3 <- function(alpha0, alpha_k, beta,
                            tau0, tau_z, alpha = NULL) {

    if (pophaz.alpha) {
      if (n_z_alpha == 0 & n_z_tau == 0) {
        theta <- c(alpha0, beta, tau0, alpha)
      } else if (n_z_alpha > 0 & n_z_tau == 0) {

        if (optim_func == "bbmle") {

          alpha_k <- eval(
            parse(text = paste("c(",
                               paste(paste("tau", 1:n_z_tau,
                                           sep = "_"),
                                     collapse = ","),
                               ")", sep = "")))
        }
        theta <- c(alpha0, alpha_k, beta, tau0, alpha)
      }else if (n_z_alpha == 0 & n_z_tau > 0) {
        if (optim_func == "bbmle") {
          tau_z <- eval(
            parse(text = paste("c(",
                               paste(paste("tau", 1:n_z_tau,
                                           sep = "_"),
                                     collapse = ","),
                               ")", sep = "")))
        }
        theta <- c(alpha0, beta, tau0, tau_z, alpha)
      }else if (n_z_alpha > 0 & n_z_tau > 0) {
        if (optim_func == "bbmle") {
          alpha_k <- eval(
            parse(text = paste("c(",
                               paste(paste("tau", 1:n_z_tau,
                                           sep = "_"),
                                     collapse = ","),
                               ")", sep = "")))
          tau_z <- eval(
            parse(text = paste("c(",
                               paste(paste("tau", 1:n_z_tau,
                                           sep = "_"),
                                     collapse = ","),
                               ")", sep = "")))
        }
        theta <- c(alpha0, alpha_k, beta, tau0, tau_z, alpha)
      }
    }else{
      if (n_z_alpha == 0 & n_z_tau == 0) {
        theta <- c(alpha0, beta, tau0)
      } else if (n_z_alpha > 0 & n_z_tau == 0) {

        if (optim_func == "bbmle") {


          alpha_k <- eval(parse(text = paste("c(",
                                             paste(paste("alpha", 1:n_z_alpha,
                                                         sep = "_"),
                                                   collapse = ","),
                                             ")", sep = "")))
        }
        theta <- c(alpha0, alpha_k, beta, tau0)
      } else if (n_z_alpha == 0 & n_z_tau > 0) {
        if (optim_func == "bbmle") {


          tau_z <- eval(
            parse(text = paste("c(",
                               paste(paste("tau", 1:n_z_tau,
                                           sep = "_"),
                                     collapse = ","),
                               ")", sep = "")))
        }

        theta <- c(alpha0, beta, tau0, tau_z)
      } else if (n_z_alpha > 0 & n_z_tau > 0) {
        if (optim_func == "bbmle") {
          alpha_k <- eval(
            parse(text = paste("c(",
                               paste(paste("tau", 1:n_z_tau,
                                           sep = "_"),
                                     collapse = ","),
                               ")", sep = "")))
          tau_z <- eval(
            parse(text = paste("c(",
                               paste(paste("tau", 1:n_z_tau,
                                           sep = "_"),
                                     collapse = ","),
                               ")", sep = "")))
        }
        theta <- c(alpha0, alpha_k, beta, tau0, tau_z)
      }
    }

    n_z_tau <- ncol(z_tau)
    n_z_alpha <- ncol(z_alpha)
    n_z_tau_ad <- n_z_tau - 1
    n_z_alpha_ad <- n_z_alpha - 1
    if (anyNA(theta))
      return(1e10)
    if (pophaz.alpha) {
      alpha_pophaz <- matrix(1,
                             ncol = 1,
                             nrow = length(pophaz)) %*% theta[n_z_alpha + 2 + n_z_tau + 1 + 1]
    }else{

      alpha_pophaz <- matrix(1, ncol = 1,
                             nrow = length(pophaz)) %*% 0
    }


    alpha0 <- theta[1]

    if (n_z_tau > 0 & n_z_alpha > 0) {
      alpha_k <- theta[2:(n_z_alpha + 1)]
      beta <- (theta[n_z_alpha + 2])
      tau0 <- theta[n_z_alpha + 2 + 1]
      tau_z <- theta[(n_z_alpha + 2 + 1 + 1):(n_z_alpha + 2 + n_z_tau + 1)]

      if (any(alpha0 + z_alpha %*% alpha_k <= 0) | any(tau0 + z_tau %*% tau_z <= 0)) {
        sum_minusloklik <- 1e10
      }else {
        sum_minusloklik <- -sum(d * log(c(exp(alpha_pophaz)) * pophaz +
                                          lexc_ad2(z_tau, z_alpha, x, theta)) -
                                  cumLexc_ad2(z_tau, z_alpha, x, theta) -
                                  cumpophaz * c(exp(alpha_pophaz)), na.rm = TRUE)
      }
    }else if (n_z_tau > 0 & n_z_alpha == 0) {

      beta <- (theta[2])
      tau0 <- theta[2 + 1]
      tau_z <- theta[(2 + 1 + 1):(2 + n_z_tau + 1)]

      if (any(alpha0  <= 0) | any(tau0 + z_tau %*% tau_z <= 0)) {
        sum_minusloklik <- 1e10
      }else {
        sum_minusloklik <- -sum(d * log(c(exp(alpha_pophaz)) * pophaz +
                                          lexc_ad2(z_tau, z_alpha, x, theta)) -
                                  cumLexc_ad2(z_tau, z_alpha, x, theta) -
                                  cumpophaz * c(exp(alpha_pophaz)), na.rm = TRUE)
      }
    }else if (n_z_tau == 0 & n_z_alpha > 0) {

      alpha_k <- theta[2:(n_z_alpha + 1)]
      beta <- (theta[n_z_alpha + 2])
      tau0 <- theta[n_z_alpha + 2 + 1]
      if (any(alpha0 + z_alpha %*% alpha_k <= 0) | tau0  <= 0) {
        sum_minusloklik <- 1e10
      }else {
        sum_minusloklik <- -sum(d * log(c(exp(alpha_pophaz)) * pophaz +
                                          lexc_ad2(z_tau, z_alpha, x, theta)) -
                                  cumLexc_ad2(z_tau, z_alpha, x, theta) -
                                  cumpophaz * c(exp(alpha_pophaz)), na.rm = TRUE)
      }
    }else if (n_z_tau == 0 & n_z_alpha == 0) {
      beta <- (theta[n_z_alpha + 2])
      tau0 <- theta[n_z_alpha + 2 + 1]

      if (alpha0  <= 0 | tau0 <= 0) {
        sum_minusloklik <- 1e10
      }else {
        sum_minusloklik <- -sum(d * log(c(exp(alpha_pophaz)) * pophaz +
                                          lexc_ad2(z_tau, z_alpha, x, theta)) -
                                  cumLexc_ad2(z_tau, z_alpha, x, theta) -
                                  cumpophaz * c(exp(alpha_pophaz)), na.rm = TRUE)
      }
    }




    if (anyNA(sum_minusloklik)) {
      return(1e10)
    } else {
      return(sum_minusloklik)
    }

  }




  if (!is.null(theta_lower) && !is.null(theta_upper)) {
    start_val <- gen_start_values(nvalues = 10000,
                                  min = theta_lower,
                                  max = theta_upper,
                                  npar = length(theta_init))
  } else if (NCD == "iter") {
    stop("add bounds for initial values")
  }


  if (optim_func == "optim") {
    if (method_opt == "L-BFGS-B") {
      if (gradient) {

        NCD <- ncoor_des

        if (is.null(NCD)) {
          optimized <- stats::optim(par = theta,
                                    fn = adtneh_ll_fn,
                                    gr = adtneh_ll_fn_deriv_all,
                                    method = method_opt,
                                    hessian = TRUE,
                                    lower = theta_lower,
                                    upper = theta_upper,
                                    control = list(maxit = maxit_opt,
                                                   trace = trace,
                                                   REPORT = 500))
        } else if (is.numeric(NCD)) {
          optimized_theta <- try(stats::optim(par = theta,
                                              fn = adtneh_ll_fn,
                                              gr = adtneh_ll_fn_deriv_all,
                                              method = method_opt,
                                              hessian = TRUE,
                                              lower = theta_lower,
                                              upper = theta_upper,
                                              control = list(maxit = maxit_opt,
                                                             trace = trace,
                                                             REPORT = 500))$par, TRUE)

          if (any(is.na(optimized_theta) || (inherits(optimized_theta, "try-error")))) {
            init.cd <- theta
          }  else {
            init.cd <- optimized_theta
          }


          for (j in 1:NCD) {
            message("cycle = ", j, "\n")

            for (i in 1:length(init.cd)) {
              tempf <- Vectorize(function(par){
                val <- replace(init.cd, i, par)
                return(adtneh_ll_fn(val))
              })

              tempgr <- Vectorize(function(par){
                val <- replace(init.cd, i, par)
                return(adtneh_ll_fn_deriv_all(val)[i])
              })


              new.val <- try(stats::optim(par = init.cd[i],
                                          fn = tempf,
                                          gr = tempgr,
                                          method = method_opt,
                                          hessian = TRUE,
                                          lower = theta_lower[i],
                                          upper = theta_upper[i],
                                          control = list(
                                            maxit = maxit_opt, trace = trace,
                                            REPORT = 500))$par, TRUE)
              message("new initial value of parameter number ",i," is equal to: ",new.val,"\n")
              if (anyNA(new.val) | inherits(new.val, "try-error")) {
                new.val <- (theta_lower[i] + theta_upper[i])/2
              }
              init.cd <- replace(init.cd, i,  new.val )
            }
          }

          if (any(is.na(init.cd))) {
            initG <- theta
          }  else {
            initG <- init.cd
          }

          optimized <- try(stats::optim(par = initG ,
                                        fn = adtneh_ll_fn,
                                        gr = adtneh_SDtheta,
                                        method = method_opt,
                                        hessian = TRUE,
                                        lower = theta_lower,
                                        upper = theta_upper,
                                        control = list(maxit = maxit_opt,
                                                       trace = trace,
                                                       REPORT = 500)), TRUE)



        }else if ( NCD == "iter") {
          optimized_theta <- try(stats::optim(par = theta,
                                              fn = adtneh_ll_fn,
                                              gr = adtneh_SDtheta,
                                              method = method_opt,
                                              hessian = TRUE,
                                              lower = theta_lower,
                                              upper = theta_upper,
                                              control = list(maxit = maxit_opt,
                                                             trace = trace,
                                                             REPORT = 500))$par, TRUE)

          if (any(is.na(optimized_theta) || (inherits(optimized_theta, "try-error")))) {
            init.cd <- theta
          }  else {
            init.cd <- optimized_theta
          }

          j = 0
          mydiff <- matrix(data = 1,
                           nrow = 10000,
                           ncol = 1)
          loglik_val <- matrix(data = 0,
                               nrow = 10000,
                               ncol = 1)
          theta_val <- matrix(data = 0,
                              nrow = 10000,
                              ncol = length(theta))
          while (abs(mydiff[j + 1, ]) > iter_eps && (j = j + 1) < 10000) {
            message("cycle = ", j, "\n")

            for (i in 1:length(init.cd)) {
              tempf <- Vectorize(function(par){
                val <- replace(init.cd, i, par)
                return(adtneh_ll_fn(val))
              })


              new.valpar <- try(stats::optim(par = init.cd[i],
                                             fn = tempf,
                                             method = method_opt,
                                             hessian = TRUE,
                                             lower = theta_lower[i],
                                             upper = theta_upper[i],
                                             control = list(
                                               maxit = maxit_opt, trace = trace,
                                               REPORT = 500)), TRUE)

              new.val <- try(new.valpar$par, TRUE)

              message("new initial value of parameter number ", i," is equal to: ", new.val,"\n")
              if (anyNA(new.val) | inherits(new.val, "try-error")) {
                new.val <- (theta_lower[i] + theta_upper[i])/2
              }



              init.cd <- replace(init.cd, i,  new.val )

            }
            theta_val[j, ] <-  theta
            theta_val[j + 1, ] <- try(init.cd, TRUE)
            loglik_val[j + 1,] <- try(c(new.valpar$value), TRUE)
            mydiff[j + 1, ] <- try(c(loglik_val[j + 1, ] -
                                       loglik_val[j, ]), TRUE)
            message("difference between loglikelihood at iter = ",
                j + 1,
                " equal ",
                mydiff[j + 1, ],
                ". The j+1 th logliklihood equal ",
                loglik_val[j + 1, ],
                ".\n")
            if (abs(loglik_val[(j + 1), ]) == Inf |
                abs(loglik_val[(j + 1), ]) == .Machine$double.xmax &
                loglik_val[(j + 1), ] > 0) {
              loglik_val[(j + 1), ] = 10000
              mydiff[(j + 1), ]     = 10000
            } else if (abs(loglik_val[(j + 1), ]) == Inf |
                       abs(loglik_val[(j + 1), ]) == .Machine$double.xmax &
                       loglik_val[(j + 1), ] < 0) {
              loglik_val[(j + 1), ] = -10000
              mydiff[(j + 1), ]     = -10000
            }

            if (j >= 2) {

              par(mfrow = c(1, 2))
              plot(2:(j + 1),
                   loglik_val[2:(j + 1),],
                   type = "l",
                   col = 2,
                   ylab = "loglik",
                   xlab = "iter")
              graphics::grid()
              plot(2:(j + 1),
                   abs(mydiff[2:(j + 1),]),
                   type = "l",
                   col = "green",
                   ylab = "loglik_abs_diff",
                   xlab = "iter")
              graphics::grid()
              par(mfrow = c(1, 1))
            }

          }


          if (any(is.na(init.cd))) {
            initG <- theta
          }  else {
            initG <- init.cd
          }

          optimized <- try(stats::optim(par = initG ,
                                        fn = adtneh_ll_fn,
                                        gr = adtneh_SDtheta,
                                        method = method_opt,
                                        hessian = TRUE,
                                        lower = theta_lower,
                                        upper = theta_upper,
                                        control = list(maxit = maxit_opt,
                                                       trace = trace,
                                                       REPORT = 500)), TRUE)



        }

      }else {
        NCD <- ncoor_des

        if (is.null(NCD)) {
          optimized <- stats::optim(par = theta ,
                                    fn = adtneh_ll_fn,
                                    method = method_opt,
                                    hessian = TRUE,
                                    lower = theta_lower,
                                    upper = theta_upper,
                                    control = list(maxit = maxit_opt,
                                                   trace = trace,
                                                   REPORT = 500))

        } else if (is.numeric(NCD)) {

          optimized_theta <- try(stats::optim(par = theta,
                                              fn = adtneh_ll_fn,
                                              method = method_opt,
                                              hessian = TRUE,
                                              lower = theta_lower,
                                              upper = theta_upper,
                                              control = list(maxit = maxit_opt,
                                                             trace = trace,
                                                             REPORT = 500))$par, TRUE)

          if (any(is.na(optimized_theta) || (inherits(optimized_theta, "try-error")))) {
            init.cd <- theta
          }  else {
            init.cd <- optimized_theta
          }
          for (j in 1:NCD) {
            message("cycle = ", j, "\n")

            for (i in 1:length(init.cd)) {
              tempf <- Vectorize(function(par){
                val <- replace(init.cd, i, par)
                return(adtneh_ll_fn(val))
              })

              new.val <- try(stats::optim(par = init.cd[i],
                                          fn = tempf,
                                          method = method_opt,
                                          hessian = TRUE,
                                          lower = theta_lower[i],
                                          upper = theta_upper[i],
                                          control = list(
                                            maxit = maxit_opt, trace = trace,
                                            REPORT = 500))$par, TRUE)
              message("new initial value of parameter number ",i," is equal to: ",new.val,"\n")
              if (anyNA(new.val) | inherits(new.val, "try-error")) {
                new.val <- (theta_lower[i] + theta_upper[i])/2
              }

              init.cd <- replace(init.cd, i,  new.val )
            }
          }

          if (any(is.na(init.cd))) {
            initG <- theta
          }  else {
            initG <- init.cd
          }

          optimized <- try(stats::optim(par = initG ,
                                        fn = adtneh_ll_fn,
                                        method = method_opt,
                                        hessian = TRUE,
                                        lower = theta_lower,
                                        upper = theta_upper,
                                        control = list(maxit = maxit_opt,
                                                       trace = trace,
                                                       REPORT = 500)), TRUE)

        } else if (NCD == "iter") {
          optimized_theta <- try(stats::optim(par = theta,
                                              fn = adtneh_ll_fn,
                                              gr = adtneh_SDtheta,
                                              method = method_opt,
                                              hessian = TRUE,
                                              lower = theta_lower,
                                              upper = theta_upper,
                                              control = list(maxit = maxit_opt,
                                                             trace = trace,
                                                             REPORT = 500))$par,
                                 TRUE)

          if (any(is.na(optimized_theta) || (inherits(optimized_theta, "try-error")))) {
            init.cd <- theta
          }  else {
            init.cd <- optimized_theta
          }

          j = 0
          mydiff <- matrix(data = 1,
                           nrow = 10000,
                           ncol = 1)
          loglik_val <- matrix(data = 0,
                               nrow = 10000,
                               ncol = 1)
          theta_val <- matrix(data = 0,
                              nrow = 10000,
                              ncol = length(theta))
          while (abs(mydiff[j + 1, ]) > iter_eps && (j = j + 1) < 10000) {
            message("cycle = ", j, "\n")

            for (i in 1:length(init.cd)) {
              tempf <- Vectorize(function(par){
                val <- replace(init.cd, i, par)
                return(adtneh_ll_fn(val))
              })


              new.valpar <- try(stats::optim(par = init.cd[i],
                                             fn = tempf,
                                             method = method_opt,
                                             hessian = TRUE,
                                             lower = theta_lower[i],
                                             upper = theta_upper[i],
                                             control = list(
                                               maxit = maxit_opt, trace = trace,
                                               REPORT = 500)), TRUE)

              new.val <- try(new.valpar$par, TRUE)

              message("new initial value of parameter number ",i," is equal to: ",new.val,"\n")
              if (anyNA(new.val) | inherits(new.val, "try-error")) {
                new.val <- (theta_lower[i] + theta_upper[i])/2
              }



              init.cd <- replace(init.cd, i,  new.val )

            }
            theta_val[j, ] <-  theta
            theta_val[j + 1, ] <- try(init.cd, TRUE)
            loglik_val[j + 1,] <- try(c(new.valpar$value), TRUE)
            mydiff[j + 1, ] <- try(c(loglik_val[j + 1, ] -
                                       loglik_val[j, ]), TRUE)
            message("difference between loglikelihood at iter = ",
                j + 1,
                " equal ",
                mydiff[j + 1, ],
                ". The j+1 th logliklihood equal ",
                loglik_val[j + 1, ],
                ".\n")
            if (abs(loglik_val[(j + 1), ]) == Inf |
                abs(loglik_val[(j + 1), ]) == .Machine$double.xmax &
                loglik_val[(j + 1), ] > 0) {
              loglik_val[(j + 1), ] = 10000
              mydiff[(j + 1), ]     = 10000
            } else if (abs(loglik_val[(j + 1), ]) == Inf |
                       abs(loglik_val[(j + 1), ]) == .Machine$double.xmax  &
                       loglik_val[(j + 1), ] < 0) {
              loglik_val[(j + 1), ] = -10000
              mydiff[(j + 1), ]     = -10000
            }

            if (j >= 2) {
              par(mfrow = c(1, 2))
              plot(2:(j + 1),
                   loglik_val[2:(j + 1),],
                   type = "l",
                   col = 2,
                   ylab = "loglik",
                   xlab = "iter")
              graphics::grid()
              plot(2:(j + 1),
                   abs(mydiff[2:(j + 1),]),
                   type = "l",
                   col = "green",
                   ylab = "loglik_abs_diff",
                   xlab = "iter")
              graphics::grid()
              par(mfrow = c(1, 1))
            }

          }


          if (any(is.na(init.cd))) {
            initG <- theta
          }  else {
            initG <- init.cd
          }

          optimized <- try(stats::optim(par = initG ,
                                        fn = adtneh_ll_fn,
                                        method = method_opt,
                                        hessian = TRUE,
                                        lower = theta_lower,
                                        upper = theta_upper,
                                        control = list(maxit = maxit_opt,
                                                       trace = trace,
                                                       REPORT = 500)), TRUE)



        }
      }




    } else {
      if (gradient) {
        NCD <- ncoor_des
        if (is.null(NCD)) {
          optimized <- stats::optim(par = theta,
                                    fn = adtneh_ll_fn,
                                    gr = adtneh_SDtheta,
                                    method = method_opt,
                                    control = list(pgtol = 1e-15,
                                                   maxit = maxit_opt,
                                                   trace = trace,
                                                   REPORT = 500))
        } else if (is.numeric(NCD)) {
          optimized_theta <- try(stats::optim(par = theta,
                                              fn = adtneh_ll_fn,
                                              gr = adtneh_SDtheta,
                                              method = method_opt,
                                              hessian = TRUE,
                                              lower = theta_lower,
                                              upper = theta_upper,
                                              control = list(maxit = maxit_opt,
                                                             trace = trace,
                                                             REPORT = 500))$par, TRUE)

          if (any(is.na(optimized_theta) || (inherits(optimized_theta, "try-error")))) {
            init.cd <- theta
          }  else {
            init.cd <- optimized_theta
          }


          for (j in 1:NCD) {
            message("cycle = ", j, "\n")
            for (i in 1:length(init.cd)) {
              tempf <- Vectorize(function(par){
                val <- replace(init.cd, i, par)
                return(adtneh_ll_fn(val))
              })

              new.val <- try(stats::optim(par = init.cd[i],
                                          fn = tempf,
                                          method = method_opt,
                                          control = list(pgtol = 1e-15,
                                                         maxit = maxit_opt,
                                                         trace = trace,
                                                         REPORT = 500))$par, TRUE)
              message("new initial value of parameter number ",i," is equal to: ",new.val,"\n")
              if (anyNA(new.val) | inherits(new.val, "try-error")) {
                new.val <- init.cd[i] + 1
              }
              init.cd <- replace(init.cd, i,  new.val )
            }
          }

          if (any(is.na(init.cd))) {
            initG <- theta
          }  else {
            initG <- init.cd
          }

          optimized <- try(stats::optim(par = initG ,
                                        fn = adtneh_ll_fn,
                                        gr = adtneh_SDtheta,
                                        method = method_opt,
                                        control = list(pgtol = 1e-15,
                                                       maxit = maxit_opt,
                                                       trace = trace,
                                                       REPORT = 500)), TRUE)




        } else if (NCD == "iter") {
          optimized_theta <- try(stats::optim(par = theta,
                                              fn = adtneh_ll_fn,
                                              gr = adtneh_SDtheta,
                                              method = method_opt,
                                              hessian = TRUE,
                                              lower = theta_lower,
                                              upper = theta_upper,
                                              control = list(maxit = maxit_opt,
                                                             trace = trace,
                                                             REPORT = 500))$par, TRUE)

          if (any(is.na(optimized_theta) || (inherits(optimized_theta, "try-error")))) {
            init.cd <- theta
          }  else {
            init.cd <- optimized_theta
          }

          j = 0
          mydiff <- matrix(data = 1,
                           nrow = 10000,
                           ncol = 1)
          loglik_val <- matrix(data = 0,
                               nrow = 10000,
                               ncol = 1)
          theta_val <- matrix(data = 0,
                              nrow = 10000,
                              ncol = length(theta))
          while (abs(mydiff[j + 1, ]) > iter_eps && (j = j + 1) < 10000) {
            message("cycle = ", j, "\n")

            for (i in 1:length(init.cd)) {
              tempf <- Vectorize(function(par){
                val <- replace(init.cd, i, par)
                return(adtneh_ll_fn(val))
              })


              new.valpar <- try(stats::optim(par = init.cd[i],
                                             fn = tempf,
                                             method = method_opt,
                                             hessian = TRUE,
                                             control = list(
                                               maxit = maxit_opt, trace = trace,
                                               REPORT = 500)), TRUE)

              new.val <- try(new.valpar$par, TRUE)

              message("new initial value of parameter number ",i," is equal to: ",new.val,"\n")
              if (anyNA(new.val) | inherits(new.val, "try-error")) {
                new.val <- init.cd[i] + 1
              }



              init.cd <- replace(init.cd, i,  new.val )

            }
            theta_val[j, ] <-  theta
            theta_val[j + 1, ] <- try(init.cd, TRUE)
            loglik_val[j + 1,] <- try(c(new.valpar$value), TRUE)
            mydiff[j + 1, ] <- try(c(loglik_val[j + 1, ] -
                                       loglik_val[j, ]), TRUE)
            message("difference between loglikelihood at iter = ",
                j + 1,
                " equal ",
                mydiff[j + 1, ],
                ". The j+1 th logliklihood equal ",
                loglik_val[j + 1, ],
                ".\n")

            if (abs(loglik_val[(j + 1), ]) == Inf  |
                abs(loglik_val[(j + 1), ]) == .Machine$double.xmax &
                loglik_val[(j + 1), ] > 0) {
              loglik_val[(j + 1), ] = 10000
              mydiff[(j + 1), ]     = 10000
            } else if (abs(loglik_val[(j + 1), ]) == Inf |
                       abs(loglik_val[(j + 1), ]) == .Machine$double.xmax &
                       loglik_val[(j + 1), ] < 0) {
              loglik_val[(j + 1), ] = -10000
              mydiff[(j + 1), ]     = -10000
            }

            if (j >= 2) {
              par(mfrow = c(1, 2))
              plot(2:(j + 1),
                   loglik_val[2:(j + 1),],
                   type = "l",
                   col = 2,
                   ylab = "loglik",
                   xlab = "iter")
              graphics::grid()
              plot(2:(j + 1),
                   abs(mydiff[2:(j + 1),]),
                   type = "l",
                   col = "green",
                   ylab = "loglik_abs_diff",
                   xlab = "iter")
              graphics::grid()
              par(mfrow = c(1, 1))
            }
          }


          if (any(is.na(init.cd))) {
            initG <- theta
          }  else {
            initG <- init.cd
          }

          optimized <- try(stats::optim(par = initG ,
                                        fn = adtneh_ll_fn,
                                        gr = adtneh_SDtheta,
                                        method = method_opt,
                                        hessian = TRUE,
                                        control = list(maxit = maxit_opt,
                                                       trace = trace,
                                                       REPORT = 500)), TRUE)


        }
      }
      else{

        NCD <- ncoor_des
        if (is.null(NCD)) {
          optimized <- stats::optim(par = theta,
                                    fn = adtneh_ll_fn,
                                    method = method_opt,
                                    control = list(pgtol = 1e-15,
                                                   maxit = maxit_opt,
                                                   trace = trace,
                                                   REPORT = 500))
        } else if (is.numeric(NCD)) {
          optimized_theta <- try(stats::optim(par = theta,
                                              fn = adtneh_ll_fn,
                                              method = method_opt,
                                              hessian = TRUE,
                                              control = list(maxit = maxit_opt,
                                                             trace = trace,
                                                             REPORT = 500))$par, TRUE)

          if (any(is.na(optimized_theta) || (inherits(optimized_theta, "try-error")))) {
            init.cd <- theta
          }  else {
            init.cd <- optimized_theta
          }

          for (j in 1:NCD) {
            message("cycle = ", j, "\n")

            for (i in 1:length(init.cd)) {
              tempf <- Vectorize(function(par){
                val <- replace(init.cd, i, par)
                return(adtneh_ll_fn(val))
              })
              new.val <- try(stats::optim(par = init.cd[i],
                                          fn = tempf,
                                          method = method_opt,
                                          control = list(pgtol = 1e-15,
                                                         maxit = maxit_opt,
                                                         trace = trace,
                                                         REPORT = 500))$par, TRUE)
              message("new initial value of parameter number ",i," is equal to: ",new.val,"\n")
              if (anyNA(new.val) | inherits(new.val, "try-error")) {
                new.val <- init.cd[i] + 1
              }
              init.cd <- replace(init.cd, i,  new.val )
            }
          }

          if (any(is.na(init.cd))) {
            initG <- theta
          }  else {
            initG <- init.cd
          }

          optimized <- try(stats::optim(par = initG ,
                                        fn = adtneh_ll_fn,
                                        method = method_opt,
                                        control = list(pgtol = 1e-15,
                                                       maxit = maxit_opt,
                                                       trace = trace,
                                                       REPORT = 500)), TRUE)




        } else if ( NCD == "iter") {
          optimized_theta <- try(stats::optim(par = theta,
                                              fn = adtneh_ll_fn,
                                              method = method_opt,
                                              hessian = TRUE,
                                              control = list(maxit = maxit_opt,
                                                             trace = trace,
                                                             REPORT = 500))$par, TRUE)

          if (any(is.na(optimized_theta) || (inherits(optimized_theta, "try-error")))) {
            init.cd <- theta
          }  else {
            init.cd <- optimized_theta
          }

          j = 0
          mydiff <- matrix(data = 1,
                           nrow = 10000,
                           ncol = 1)
          loglik_val <- matrix(data = 0,
                               nrow = 10000,
                               ncol = 1)
          theta_val <- matrix(data = 0,
                              nrow = 10000,
                              ncol = length(theta))
          while (abs(mydiff[j + 1, ]) > iter_eps && (j = j + 1) < 10000) {
            message("cycle = ", j, "\n")

            for (i in 1:length(init.cd)) {
              tempf <- Vectorize(function(par){
                val <- replace(init.cd, i, par)
                return(adtneh_ll_fn(val))
              })


              new.valpar <- try(stats::optim(par = init.cd[i],
                                             fn = tempf,
                                             method = method_opt,
                                             hessian = TRUE,
                                             control = list(
                                               maxit = maxit_opt, trace = trace,
                                               REPORT = 500)), TRUE)

              new.val <- try(new.valpar$par, TRUE)

              message("new initial value of parameter number ",i," is equal to: ",new.val,"\n")
              if (anyNA(new.val) | inherits(new.val, "try-error")) {
                new.val <- init.cd[i] + 1
              }



              init.cd <- replace(init.cd, i,  new.val )

            }
            theta_val[j, ] <-  theta
            theta_val[j + 1, ] <- try(init.cd, TRUE)
            loglik_val[j + 1,] <- try(c(new.valpar$value), TRUE)
            mydiff[j + 1, ] <- try(c(loglik_val[j + 1, ] -
                                       loglik_val[j, ]), TRUE)
            message("difference between loglikelihood at iter = ",
                j + 1,
                " equal ",
                mydiff[j + 1, ],
                ". The j+1 th logliklihood equal ",
                loglik_val[j + 1, ],
                ".\n")

            if (abs(loglik_val[(j + 1), ]) == Inf |
                abs(loglik_val[(j + 1), ]) == .Machine$double.xmax &
                loglik_val[(j + 1), ] > 0) {
              loglik_val[(j + 1), ] = 10000
              mydiff[(j + 1), ]     = 10000
            } else if (abs(loglik_val[(j + 1), ]) == Inf |
                       abs(loglik_val[(j + 1), ]) == .Machine$double.xmax &
                       loglik_val[(j + 1), ] < 0) {
              loglik_val[(j + 1), ] = -10000
              mydiff[(j + 1), ]     = -10000
            }

            if (j >= 2) {
              par(mfrow = c(1, 2))
              plot(2:(j + 1),
                   loglik_val[2:(j + 1),],
                   type = "l",
                   col = 2,
                   ylab = "loglik",
                   xlab = "iter")
              graphics::grid()
              plot(2:(j + 1),
                   abs(mydiff[2:(j + 1),]),
                   type = "l",
                   col = "green",
                   ylab = "loglik_abs_diff",
                   xlab = "iter")
              graphics::grid()
              par(mfrow = c(1, 1))
            }
          }


          if (any(is.na(init.cd))) {
            initG <- theta
          }  else {
            initG <- init.cd
          }

          optimized <- try(stats::optim(par = initG ,
                                        fn = adtneh_ll_fn,
                                        method = method_opt,
                                        hessian = TRUE,
                                        control = list(maxit = maxit_opt,
                                                       trace = trace,
                                                       REPORT = 500)), TRUE)


        }
      }}

  } else if (optim_func == "bbmle") {

    if (method_opt == "L-BFGS-B" |
        method_opt == "Rvmmin" |
        method_opt == "Rcgmin" |
        method_opt == "nlminb" |
        method_opt == "bobyqa" ) {
      if (gradient) {

        NCD <- ncoor_des

        if (is.null(NCD)) {
          optimized <- stats::optim(par = theta,
                                    fn = adtneh_ll_fn,
                                    gr = adtneh_SDtheta,
                                    method = method_opt,
                                    hessian = TRUE,
                                    lower = theta_lower,
                                    upper = theta_upper,
                                    control = list(maxit = maxit_opt,
                                                   trace = trace,
                                                   REPORT = 500))
        } else {

          init.cd <- theta
          for (j in 1:NCD) {
            message("cycle = ", j, "\n")

            for (i in 1:length(init.cd)) {
              tempf <- Vectorize(function(par){
                val <- replace(init.cd, i, par)
                return(adtneh_ll_fn(val))
              })

              new.val <- try(stats::optim(par = init.cd[i],
                                          fn = tempf,
                                          method = method_opt,
                                          hessian = TRUE,
                                          lower = theta_lower[i],
                                          upper = theta_upper[i],
                                          control = list(
                                            maxit = maxit_opt, trace = trace,
                                            REPORT = 500))$par, TRUE)
              message("new initial value of parameter number ",i," is equal to: ",new.val,"\n")
              if (anyNA(new.val) | inherits(new.val, "try-error")) {
                new.val <- (theta_lower[i] + theta_upper[i])/2
              }
              init.cd <- replace(init.cd, i,  new.val )
            }
          }

          if (any(is.na(init.cd))) {
            initG <- theta
          }  else {
            initG <- init.cd
          }

          optimized <- try(stats::optim(par = initG ,
                                        fn = adtneh_ll_fn,
                                        gr = adtneh_SDtheta,
                                        method = method_opt,
                                        hessian = TRUE,
                                        lower = theta_lower,
                                        upper = theta_upper,
                                        control = list(maxit = maxit_opt,
                                                       trace = trace,
                                                       REPORT = 500)), TRUE)



        }



      }else {
        NCD <- ncoor_des
        if (pophaz.alpha) {
          alpha <- theta[length(theta)]

        }else {
          alpha <- 1
        }

        if (is.null(NCD)) {
          par_adtneh2 <- genpar_adtneh2(theta, z_tau, z_alpha)

          if (n_z_tau > 0 & n_z_alpha > 0) {
            alpha0 <- par_adtneh2["alpha0"]
            alpha_k <- par_adtneh2[
              c(which(stringr::str_detect(names(par_adtneh2),"alpha_k")))]
            beta <- par_adtneh2["beta"]
            tau0 <- par_adtneh2["tau0"]
            tau_z <- par_adtneh2[
              c(which(stringr::str_detect(names(par_adtneh2),"tau_z")))]
            if (pophaz.alpha) {

              start <- c(alpha0 = as.numeric(alpha0),
                         alpha_k = as.numeric(alpha_k),
                         beta = as.numeric(beta),
                         tau0 = as.numeric(tau0),
                         tau_z = as.numeric(tau_z),
                         alpha = as.numeric(alpha))
            }else{

              start <- c(alpha0 = as.numeric(alpha0),
                         alpha_k = as.numeric(alpha_k),
                         beta = as.numeric(beta),
                         tau0 = as.numeric(tau0),
                         tau_z = as.numeric(tau_z))
            }
            names(start) <- names(theta2)



          }else if (n_z_tau > 0 & n_z_alpha == 0) {
            alpha0 <- par_adtneh2["alpha0"]
            beta <- par_adtneh2["beta"]
            tau0 <- par_adtneh2["tau0"]
            tau_z <- par_adtneh2[
              c(which(stringr::str_detect(names(par_adtneh2),"tau_z")))]

            if (pophaz.alpha) {
              start <- c(alpha0 = as.numeric(alpha0),
                         beta = as.numeric(beta),
                         tau0 = as.numeric(tau0),
                         tau_z = as.numeric(tau_z),
                         alpha = as.numeric(alpha))
            }else {
              start <- c(alpha0 = as.numeric(alpha0),
                         beta = as.numeric(beta),
                         tau0 = as.numeric(tau0),
                         tau_z = as.numeric(tau_z))
            }

            names(start) <- names(theta2)

          }else if (n_z_tau == 0 & n_z_alpha > 0) {

            alpha0 <- par_adtneh2["alpha0"]
            alpha_k <- par_adtneh2[
              c(which(stringr::str_detect(names(par_adtneh2),"alpha_k")))]
            beta <- par_adtneh2["beta"]
            tau0 <- par_adtneh2["tau0"]
            if (pophaz.alpha) {
              start <- c(alpha0 = as.numeric(alpha0),
                         alpha_k = as.numeric(alpha_k),
                         beta = as.numeric(beta),
                         tau0 = as.numeric(tau0),
                         alpha = as.numeric(alpha))
            }else {
              start <- c(alpha0 = as.numeric(alpha0),
                         alpha_k = as.numeric(alpha_k),
                         beta = as.numeric(beta),
                         tau0 = as.numeric(tau0))
            }

            names(start) <- names(theta2)

          }else if (n_z_tau == 0 & n_z_alpha == 0) {

            alpha0 <- par_adtneh2["alpha0"]
            beta <- par_adtneh2["beta"]
            tau0 <- par_adtneh2["tau0"]

            if (pophaz.alpha) {
              start <- c(alpha0 = as.numeric(alpha0),
                         beta = as.numeric(beta),
                         tau0 = as.numeric(tau0),
                         alpha = as.numeric(alpha))
            }else {
              start <- c(alpha0 = as.numeric(alpha0),
                         beta = as.numeric(beta),
                         tau0 = as.numeric(tau0))
            }

            names(start) <- names(theta2)

          }


          add.arguments <- function(f, params) {
            t <- paste("arg <- alist(",
                       paste(sapply(1:length(params),
                                    function(i)
                                      paste(names(params)[i], "=", sep = "")),
                             collapse = ","),
                       ")", sep = "")
            formals(f) <- eval(parse(text = t))
            f
          }

          adtneh_ll_fn2noalpha2 <- add.arguments(adtneh_ll_fn3, start)
          bbmle::parnames(adtneh_ll_fn2noalpha2) <- names(theta2)


          start2 <- paste("list(",
                          paste(sapply(1:length(start),
                                       function(i)
                                         paste(names(start)[i], "=",(start)[i], sep = "")),
                                collapse = ","),
                          ")", sep = "")
          start3 <- eval(parse(text = start2))

          theta_lower2 <- paste("list(",
                                paste(sapply(1:length(start),
                                             function(i)
                                               paste(names(start)[i], "=",
                                                     (theta_lower)[i], sep = "")),
                                      collapse = ","),
                                ")", sep = "")

          theta_upper2 <- paste("list(",
                                paste(sapply(1:length(start),
                                             function(i)
                                               paste(names(start)[i], "=",
                                                     (theta_upper)[i], sep = "")),
                                      collapse = ","),
                                ")", sep = "")

          theta_lower3 <- eval(parse(text = theta_lower2))
          theta_upper3 <- eval(parse(text = theta_upper2))

          optimized <- bbmle::mle2(
            start = start3,
            minuslogl = adtneh_ll_fn2noalpha2,
            fixed = fixed,
            optimizer = optimizer,
            method = method_opt,
            lower = unlist(theta_lower3),
            upper = unlist(theta_upper3),
            control = list(maxit = maxit_opt,
                           trace = trace))

        } else {
          init.cd <- theta
          for (j in 1:NCD) {
            message("cycle = ", j, "\n")

            for (i in 1:length(init.cd)) {
              tempf <- Vectorize(function(par){
                val <- replace(init.cd, i, par)
                return(adtneh_ll_fn(val))
              })

              new.val <- try(stats::optim(par = init.cd[i],
                                          fn = tempf,
                                          method = method_opt,
                                          hessian = TRUE,
                                          lower = theta_lower[i],
                                          upper = theta_upper[i],
                                          control = list(
                                            maxit = maxit_opt, trace = trace,
                                            REPORT = 500))$par, TRUE)
              message("new initial value of parameter number ",i," is equal to: ",new.val,"\n")
              if (anyNA(new.val) | inherits(new.val, "try-error")) {
                new.val <- (theta_lower[i] + theta_upper[i])/2
              }

              init.cd <- replace(init.cd, i,  new.val )
            }
          }

          if (any(is.na(init.cd))) {
            initG <- theta
          }  else {
            initG <- init.cd
          }

          optimized <- try(stats::optim(par = initG ,
                                        fn = adtneh_ll_fn,
                                        method = method_opt,
                                        hessian = TRUE,
                                        lower = theta_lower,
                                        upper = theta_upper,
                                        control = list(maxit = maxit_opt,
                                                       trace = trace,
                                                       REPORT = 500)), TRUE)

        }
      }




    }
    else
    { if (pophaz.alpha) {
      alpha <- theta[length(theta)]

    }else {
      alpha <- 1
    }

      if (method_opt == "all.methods") {

        NCD <- ncoor_des

        if (is.null(NCD)) {
          par_adtneh2 <- genpar_adtneh2(theta, z_tau, z_alpha)

          if (n_z_tau > 0 & n_z_alpha > 0) {
            alpha0 <- par_adtneh2["alpha0"]
            alpha_k <- par_adtneh2[
              c(which(stringr::str_detect(names(par_adtneh2),"alpha_k")))]
            beta <- par_adtneh2["beta"]
            tau0 <- par_adtneh2["tau0"]
            tau_z <- par_adtneh2[
              c(which(stringr::str_detect(names(par_adtneh2),"tau_z")))]
            if (pophaz.alpha) {
              start <- c(alpha0 = as.numeric(alpha0),
                         alpha_k = as.numeric(alpha_k),
                         beta = as.numeric(beta),
                         tau0 = as.numeric(tau0),
                         tau_z = as.numeric(tau_z),
                         alpha = as.numeric(alpha))
            }else{
              start <- c(alpha0 = as.numeric(alpha0),
                         alpha_k = as.numeric(alpha_k),
                         beta = as.numeric(beta),
                         tau0 = as.numeric(tau0),
                         tau_z = as.numeric(tau_z))
            }
            names(start) <- names(theta2)



          }else if (n_z_tau > 0 & n_z_alpha == 0) {
            alpha0 <- par_adtneh2["alpha0"]
            beta <- par_adtneh2["beta"]
            tau0 <- par_adtneh2["tau0"]
            tau_z <- par_adtneh2[
              c(which(stringr::str_detect(names(par_adtneh2),"tau_z")))]

            if (pophaz.alpha) {
              start <- c(alpha0 = as.numeric(alpha0),
                         beta = as.numeric(beta),
                         tau0 = as.numeric(tau0),
                         tau_z = as.numeric(tau_z),
                         alpha = as.numeric(alpha))
            }else {
              start <- c(alpha0 = as.numeric(alpha0),
                         beta = as.numeric(beta),
                         tau0 = as.numeric(tau0),
                         tau_z = as.numeric(tau_z))
            }

            names(start) <- names(theta2)

          }else if (n_z_tau == 0 & n_z_alpha > 0) {

            alpha0 <- par_adtneh2["alpha0"]
            alpha_k <- par_adtneh2[
              c(which(stringr::str_detect(names(par_adtneh2),"alpha_k")))]
            beta <- par_adtneh2["beta"]
            tau0 <- par_adtneh2["tau0"]
            if (pophaz.alpha) {
              start <- c(alpha0 = as.numeric(alpha0),
                         alpha_k = as.numeric(alpha_k),
                         beta = as.numeric(beta),
                         tau0 = as.numeric(tau0),
                         alpha = as.numeric(alpha))
            }else {
              start <- c(alpha0 = as.numeric(alpha0),
                         alpha_k = as.numeric(alpha_k),
                         beta = as.numeric(beta),
                         tau0 = as.numeric(tau0))
            }

            names(start) <- names(theta2)

          }else if (n_z_tau == 0 & n_z_alpha == 0) {

            alpha0 <- par_adtneh2["alpha0"]
            beta <- par_adtneh2["beta"]
            tau0 <- par_adtneh2["tau0"]

            if (pophaz.alpha) {
              start <- c(alpha0 = as.numeric(alpha0),
                         beta = as.numeric(beta),
                         tau0 = as.numeric(tau0),
                         alpha = as.numeric(alpha))
            }else {
              start <- c(alpha0 = as.numeric(alpha0),
                         beta = as.numeric(beta),
                         tau0 = as.numeric(tau0))
            }

            names(start) <- names(theta2)

          }


          add.arguments <- function(f, params) {
            t <- paste("arg <- alist(",
                       paste(sapply(1:length(params),
                                    function(i)
                                      paste(names(params)[i], "=", sep = "")),
                             collapse = ","),
                       ")", sep = "")
            formals(f) <- eval(parse(text = t))
            f
          }

          adtneh_ll_fn2noalpha2 <- add.arguments(adtneh_ll_fn3, start)
          bbmle::parnames(adtneh_ll_fn2noalpha2) <- names(theta2)

          start2 <- paste("list(",
                          paste(sapply(1:length(start),
                                       function(i)
                                         paste(names(start)[i], "=",(start)[i], sep = "")),
                                collapse = ","),
                          ")", sep = "")
          start3 <- eval(parse(text = start2))

          theta_lower2 <- paste("list(",
                                paste(sapply(1:length(start),
                                             function(i)
                                               paste(names(start)[i], "=",
                                                     (theta_lower)[i], sep = "")),
                                      collapse = ","),
                                ")", sep = "")

          theta_upper2 <- paste("list(",
                                paste(sapply(1:length(start),
                                             function(i)
                                               paste(names(start)[i], "=",
                                                     (theta_upper)[i], sep = "")),
                                      collapse = ","),
                                ")", sep = "")

          theta_lower3 <- eval(parse(text = theta_lower2))
          theta_upper3 <- eval(parse(text = theta_upper2))




          optimized2 <- optimx::optimx(
            par = unlist(start3),
            fn = adtneh_ll_fn,
            lower = unlist(theta_lower3),
            upper = unlist(theta_upper3),
            control = list(all.methods = TRUE,
                           maximize = FALSE,
                           trace = trace,
                           save.failures = TRUE))

          optimized <- bbmle::mle2(
            start = start3,
            minuslogl = adtneh_ll_fn2noalpha2,
            fixed = fixed,
            lower = unlist(theta_lower3),
            upper = unlist(theta_upper3),
            control = list(all.methods = TRUE,
                           maximize = FALSE,
                           trace = trace,
                           save.failures = TRUE))


        } else {
          init.cd <- theta
          for (j in 1:NCD) {
            message("cycle = ", j, "\n")

            for (i in 1:length(init.cd)) {
              tempf <- Vectorize(function(par){
                val <- replace(init.cd, i, par)
                return(adtneh_ll_fn(val))
              })

              new.val <- try(stats::optim(par = init.cd[i],
                                          fn = tempf,
                                          method = method_opt,
                                          hessian = TRUE,
                                          lower = theta_lower[i],
                                          upper = theta_upper[i],
                                          control = list(
                                            maxit = maxit_opt, trace = trace,
                                            REPORT = 500))$par, TRUE)
              message("new initial value of parameter number ",i," is equal to: ",new.val,"\n")
              if (anyNA(new.val) | inherits(new.val, "try-error")) {
                new.val <- (theta_lower[i] + theta_upper[i])/2
              }

              init.cd <- replace(init.cd, i,  new.val )
            }
          }

          if (any(is.na(init.cd))) {
            initG <- theta
          }  else {
            initG <- init.cd
          }

          optimized <- try(stats::optim(par = initG ,
                                        fn = adtneh_ll_fn,
                                        method = method_opt,
                                        hessian = TRUE,
                                        lower = theta_lower,
                                        upper = theta_upper,
                                        control = list(maxit = maxit_opt,
                                                       trace = trace,
                                                       REPORT = 500)), TRUE)

        }

      }else{
        if (gradient) {


          NCD <- ncoor_des
          if (is.null(NCD)) {
            optimized <- stats::optim(par = theta,
                                      fn = adtneh_ll_fn,
                                      gr = adtneh_SDtheta,
                                      method = method_opt,
                                      control = list(pgtol = 1e-15,
                                                     maxit = maxit_opt,
                                                     trace = trace,
                                                     REPORT = 500))
          } else{
            init.cd <- theta
            for (j in 1:NCD) {
              message("cycle = ", j, "\n")
              for (i in 1:length(init.cd)) {
                tempf <- Vectorize(function(par){
                  val <- replace(init.cd, i, par)
                  return(adtneh_ll_fn(val))
                })

                new.val <- try(stats::optim(par = init.cd[i],
                                            fn = tempf,
                                            method = method_opt,
                                            control = list(pgtol = 1e-15,
                                                           maxit = maxit_opt,
                                                           trace = trace,
                                                           REPORT = 500))$par, TRUE)
                message("new initial value of parameter number ",i," is equal to: ",new.val,"\n")
                if (anyNA(new.val) | inherits(new.val, "try-error")) {
                  new.val <- (theta_lower[i] + theta_upper[i])/2
                }
                init.cd <- replace(init.cd, i,  new.val )
              }
            }

            if (any(is.na(init.cd))) {
              initG <- theta
            }  else {
              initG <- init.cd
            }

            optimized <- try(stats::optim(par = initG ,
                                          fn = adtneh_ll_fn,
                                          gr = adtneh_SDtheta,
                                          method = method_opt,
                                          control = list(pgtol = 1e-15,
                                                         maxit = maxit_opt,
                                                         trace = trace,
                                                         REPORT = 500)), TRUE)




          }

        }else{

          par_adtneh2 <- genpar_adtneh2(theta, z_tau, z_alpha)

          if (n_z_tau > 0 & n_z_alpha > 0) {
            alpha0 <- par_adtneh2["alpha0"]
            alpha_k <- par_adtneh2[
              c(which(stringr::str_detect(names(par_adtneh2),"alpha_k")))]
            beta <- par_adtneh2["beta"]
            tau0 <- par_adtneh2["tau0"]
            tau_z <- par_adtneh2[
              c(which(stringr::str_detect(names(par_adtneh2),"tau_z")))]
            if (pophaz.alpha) {
              start <- c(alpha0 = as.numeric(alpha0),
                         alpha_k = as.numeric(alpha_k),
                         beta = as.numeric(beta),
                         tau0 = as.numeric(tau0),
                         tau_z = as.numeric(tau_z),
                         alpha = as.numeric(alpha))
            }else{
              start <- c(alpha0 = as.numeric(alpha0),
                         alpha_k = as.numeric(alpha_k),
                         beta = as.numeric(beta),
                         tau0 = as.numeric(tau0),
                         tau_z = as.numeric(tau_z))
            }
            names(start) <- names(theta2)



          }else if (n_z_tau > 0 & n_z_alpha == 0) {
            alpha0 <- par_adtneh2["alpha0"]
            beta <- par_adtneh2["beta"]
            tau0 <- par_adtneh2["tau0"]
            tau_z <- par_adtneh2[
              c(which(stringr::str_detect(names(par_adtneh2),"tau_z")))]

            if (pophaz.alpha) {
              start <- c(alpha0 = as.numeric(alpha0),
                         beta = as.numeric(beta),
                         tau0 = as.numeric(tau0),
                         tau_z = as.numeric(tau_z),
                         alpha = as.numeric(alpha))
            }else {
              start <- c(alpha0 = as.numeric(alpha0),
                         beta = as.numeric(beta),
                         tau0 = as.numeric(tau0),
                         tau_z = as.numeric(tau_z))
            }

            names(start) <- names(theta2)

          }else if (n_z_tau == 0 & n_z_alpha > 0) {

            alpha0 <- par_adtneh2["alpha0"]
            alpha_k <- par_adtneh2[
              c(which(stringr::str_detect(names(par_adtneh2),"alpha_k")))]
            beta <- par_adtneh2["beta"]
            tau0 <- par_adtneh2["tau0"]
            if (pophaz.alpha) {
              start <- c(alpha0 = as.numeric(alpha0),
                         alpha_k = as.numeric(alpha_k),
                         beta = as.numeric(beta),
                         tau0 = as.numeric(tau0),
                         alpha = as.numeric(alpha))
            }else {
              start <- c(alpha0 = as.numeric(alpha0),
                         alpha_k = as.numeric(alpha_k),
                         beta = as.numeric(beta),
                         tau0 = as.numeric(tau0))
            }

            names(start) <- names(theta2)

          }else if (n_z_tau == 0 & n_z_alpha == 0) {

            alpha0 <- par_adtneh2["alpha0"]
            beta <- par_adtneh2["beta"]
            tau0 <- par_adtneh2["tau0"]

            if (pophaz.alpha) {
              start <- c(alpha0 = as.numeric(alpha0),
                         beta = as.numeric(beta),
                         tau0 = as.numeric(tau0),
                         alpha = as.numeric(alpha))
            }else {
              start <- c(alpha0 = as.numeric(alpha0),
                         beta = as.numeric(beta),
                         tau0 = as.numeric(tau0))
            }

            names(start) <- names(theta2)

          }


          add.arguments <- function(f, params) {
            t <- paste("arg <- alist(",
                       paste(sapply(1:length(params),
                                    function(i)
                                      paste(names(params)[i], "=", sep = "")),
                             collapse = ","),
                       ")", sep = "")
            formals(f) <- eval(parse(text = t))
            f
          }

          adtneh_ll_fn2noalpha2 <- add.arguments(adtneh_ll_fn3, start)
          bbmle::parnames(adtneh_ll_fn2noalpha2) <- names(theta2)

          if (pophaz.alpha) {
            alpha <- theta[length(theta)]

          }else {
            alpha <- 1
          }

          start2 <- paste("list(",
                          paste(sapply(1:length(start),
                                       function(i)
                                         paste(names(start)[i], "=",(start)[i], sep = "")),
                                collapse = ","),
                          ")", sep = "")
          start3 <- eval(parse(text = start2))

          theta_lower2 <- paste("list(",
                                paste(sapply(1:length(start),
                                             function(i)
                                               paste(names(start)[i], "=",
                                                     (theta_lower)[i], sep = "")),
                                      collapse = ","),
                                ")", sep = "")

          theta_upper2 <- paste("list(",
                                paste(sapply(1:length(start),
                                             function(i)
                                               paste(names(start)[i], "=",
                                                     (theta_upper)[i], sep = "")),
                                      collapse = ","),
                                ")", sep = "")

          theta_lower3 <- eval(parse(text = theta_lower2))
          theta_upper3 <- eval(parse(text = theta_upper2))

          optimized <- bbmle::mle2(
            start = start3,
            minuslogl = adtneh_ll_fn2noalpha2,
            fixed = fixed,
            optimizer = optimizer,
            method = method_opt,
            control = list(maxit = maxit_opt,
                           trace = trace))






        }
      }
    }

  }


  if ( inherits(optimized, "try-error")) {
    stop("The function cannot be evaluated using these initial parameters!")
  }

  if (optim_func == "bbmle") {


    par <- try(optimized@details$par, TRUE)


    AIC <- try(2* optimized@details$value + 2*length(par), TRUE)
    loglik <- try(-as.numeric(bbmle::logLik(optimized)) , TRUE)

    evaluations <- try(optimized@details$counts, TRUE)
    iterations <- try(evaluations, TRUE)
    convergence <- try(optimized@details$convergence, TRUE)
    message <- try(optimized@details$message, TRUE)
    par_adtneh2 <- genpar_adtneh2(par, z_tau, z_alpha)
    alpha0 <- par_adtneh2["alpha0"]
    alpha_k <- par_adtneh2[c(which(stringr::str_detect(names(par_adtneh2),
                                                       "alpha_k")))]
    beta <- par_adtneh2["beta"]
    tau0 <- par_adtneh2["tau0"]
    tau_z <- par_adtneh2[c(which(stringr::str_detect(names(par_adtneh2),
                                                     "tau_z")))]
    varcov_star <- optimized@vcov
    if (any(inherits(optimized@vcov, "try-error"))) {
      varcov_star <- try(solve(optimized@details$hessian), TRUE)
    }
  }else{
    par <- try(optimized$par, TRUE)
    AIC <- try(2*optimized$value + 2*length(par), TRUE)
    loglik <- try(-optimized$value, TRUE)

    evaluations <- try(optimized$counts, TRUE)
    iterations <- try(evaluations, TRUE)
    convergence <- try(optimized$convergence, TRUE)
    message <- try(optimized$message, TRUE)
    SD <- try(numDeriv::hessian(adtneh_ll_fn, par), TRUE)

    if (hessian_varcov) {
      varcov_star <- try(solve(SD), TRUE)
    } else{
      varcov_star <- try(adtneh_vartheta(z_tau,
                                         z_alpha,
                                         time = x,
                                         theta = par,
                                         d,
                                         pophaz), TRUE)
    }

    if (any(inherits(varcov_star, "try-error"))) {
      varcov_star   <- try(solve(optimized$hessian), TRUE)
    }
  }






  std_err_star <- try(sqrt(diag(varcov_star)), TRUE)

  thetaest <- matrix(par, nrow = 1)
  n_z_tau <- ncol(z_tau)
  n_z_tau_ad <- n_z_tau - 1

  if (pophaz.alpha) {
    if (n_z_tau > 0 & n_z_alpha > 0) {
      colnames(thetaest) <- c('alpha0',
                              paste('alpha', 1:n_z_alpha, colnames(z_alpha) ,sep = "_"),
                              'beta',
                              'tau0',
                              paste('tau', 1:n_z_tau, colnames(z_tau), sep = "_"),
                              'pophaz.alpha')


    }else if (n_z_tau > 0 & n_z_alpha == 0) {
      colnames(thetaest) <- c('alpha0',
                              'beta',
                              'tau0',
                              paste('tau', 1:n_z_tau, colnames(z_tau), sep = "_"),
                              'pophaz.alpha')
    }else if (n_z_tau == 0 & n_z_alpha > 0) {
      colnames(thetaest) <- c('alpha0',
                              paste('alpha', 1:n_z_alpha, colnames(z_alpha) ,sep = "_"),
                              'beta',
                              'tau0',
                              'pophaz.alpha')}else
                                if (n_z_tau == 0 & n_z_alpha == 0) {

                                  colnames(thetaest) <- c('alpha0',
                                                          'beta',
                                                          'tau0',
                                                          'pophaz.alpha')

                                }
  }else
  {
    if (n_z_tau > 0 & n_z_alpha > 0) {
      colnames(thetaest) <- c('alpha0',
                              paste('alpha', 1:n_z_alpha, colnames(z_alpha) ,sep = "_"),
                              'beta',
                              'tau0',
                              paste('tau', 1:n_z_tau, colnames(z_tau), sep = "_"))


    }else if (n_z_tau > 0 & n_z_alpha == 0) {
      colnames(thetaest) <- c('alpha0',
                              'beta',
                              'tau0',
                              paste('tau', 1:n_z_tau, colnames(z_tau), sep = "_"))
    }else if (n_z_tau == 0 & n_z_alpha > 0) {
      colnames(thetaest) <- c('alpha0',
                              paste('alpha', 1:n_z_alpha, colnames(z_alpha) ,sep = "_"),
                              'beta',
                              'tau0')
    }else if (n_z_tau == 0 & n_z_alpha == 0) {

      colnames(thetaest) <- c('alpha0',
                              'beta',
                              'tau0')
    }
  }


  id_beta <- which(colnames(thetaest) == 'beta')
  id_alpha <- which(stringr::str_detect(colnames(thetaest),
                                        "alpha"))

  id_tau <- which(stringr::str_detect(colnames(thetaest),
                                      "tau"))



  thetaest2 <- thetaest
  thetaest2[-id_beta] <- (thetaest[-id_beta])
  thetaest2[id_beta] <- thetaest[id_beta]
  thetaest2[id_alpha] <- cumsum(thetaest[id_alpha])
  thetaest2[id_tau] <- cumsum(thetaest[id_tau])


  if (pophaz.alpha) {
    id_pophaz.alpha <- which(colnames(thetaest) == 'pophaz.alpha')
    thetaest2[id_pophaz.alpha] <- exp(thetaest[id_pophaz.alpha])
  }
  Dthetaest2 <- thetaest2
  Dthetaest2[-id_beta] <- 1
  Dthetaest2[id_beta] <- 1
  if (pophaz.alpha) {
    Dthetaest2[id_pophaz.alpha] <- exp(thetaest[id_pophaz.alpha])

  }
  varcov <- try(c(Dthetaest2^2)*varcov_star,TRUE)
  std_err <- try((sqrt(diag(varcov))), TRUE)

  if (inherits(std_err, "try-error")) {
    #error print to be implemented
  }else{
    std_err <- try(matrix(std_err, ncol = length(std_err)), TRUE)

    if (pophaz.alpha) {

      if (n_z_tau > 0 & n_z_alpha > 0) {
        colnames(std_err) <- c('alpha0',
                               paste('alpha', 1:n_z_alpha, colnames(z_alpha) ,sep = "_"),
                               'beta',
                               'tau0',
                               paste('tau', 1:n_z_tau, colnames(z_tau), sep = "_"),
                               'pophaz.alpha')


      }else if (n_z_tau > 0 & n_z_alpha == 0) {
        colnames(std_err) <- c('alpha0',
                               'beta',
                               'tau0',
                               paste('tau', 1:n_z_tau, colnames(z_tau), sep = "_"),
                               'pophaz.alpha')
      }else if (n_z_tau == 0 & n_z_alpha > 0) {
        colnames(std_err) <- c('alpha0',
                               paste('alpha', 1:n_z_alpha, colnames(z_alpha) ,sep = "_"),
                               'beta',
                               'tau0',
                               'pophaz.alpha')
      }else
        if (n_z_tau == 0 & n_z_alpha == 0) {

          colnames(std_err) <- c('alpha0',
                                 'beta',
                                 'tau0',
                                 'pophaz.alpha')

        }

    }else{

      if (n_z_tau > 0 & n_z_alpha > 0) {
        colnames(std_err) <- c('alpha0',
                               paste('alpha', 1:n_z_alpha, colnames(z_alpha) ,sep = "_"),
                               'beta',
                               'tau0',
                               paste('tau', 1:n_z_tau, colnames(z_tau), sep = "_"))


      }else if (n_z_tau > 0 & n_z_alpha == 0) {
        colnames(std_err) <- c('alpha0',
                               'beta',
                               'tau0',
                               paste('tau', 1:n_z_tau, colnames(z_tau), sep = "_"))
      }else if (n_z_tau == 0 & n_z_alpha > 0) {
        colnames(std_err) <- c('alpha0',
                               paste('alpha', 1:n_z_alpha, colnames(z_alpha) ,sep = "_"),
                               'beta',
                               'tau0')
      }else if (n_z_tau == 0 & n_z_alpha == 0) {

        colnames(std_err) <- c('alpha0',
                               'beta',
                               'tau0')
      }

    }

    colnames(thetaest) <- try(paste0(colnames(thetaest),"*"),TRUE)
    std_err_star <- matrix(std_err_star, ncol = length(std_err_star))
    colnames(std_err_star) <- try(paste0(colnames(thetaest),"*"), TRUE)
  }

  if (is.null(ncoor_des)) {
    iter_coords <- 0

  }else if (is.numeric(ncoor_des)) {
    iter_coords <- ncoor_des
  }else if (ncoor_des == "iter") {
    iter_coords <- j
  }

  if (method_opt == "all.methods") {
    return(
      list(iter_coords = iter_coords,
           coefficients = thetaest,
           estimates = thetaest2,
           loglik = loglik,
           iterations = iterations,
           evaluations = evaluations,
           convergence = convergence,
           message = message,
           varcov = varcov,
           varcov_star = varcov_star,
           std_err = std_err,
           std_err_star = std_err_star,
           AIC = AIC,
           optimized2 = optimized2
      )
    )
  }else{
    return(
      list(iter_coords = iter_coords,
           coefficients = thetaest,
           estimates = thetaest2,
           loglik = loglik,
           iterations = iterations,
           evaluations = evaluations,
           convergence = convergence,
           message = message,
           varcov = varcov,
           varcov_star = varcov_star,
           std_err = std_err,
           std_err_star = std_err_star,
           AIC = AIC
      )
    )
  }


}
