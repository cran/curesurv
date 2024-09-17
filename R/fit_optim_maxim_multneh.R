fit.opt.maxim.multneh <- function(x = x,
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
                      trace = trace) {
  theta <- theta_init
  n_z_tau <- ncol(z_tau)
  n_z_alpha <- ncol(z_alpha)
  n_z_tau_ad <- n_z_tau - 1
  n_z_alpha_ad <- n_z_alpha - 1

  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))

  if (pophaz.alpha) {
    if (n_z_tau >= 1 & n_z_alpha >= 1) {
      names(theta) <- c('alpha0',
                              paste('alpha', 1:n_z_alpha, colnames(z_alpha),
                                    sep = "_"),
                        'beta',
                        'tau0',
                        paste('tau',
                              1:n_z_tau, colnames(z_tau), sep = "_"),
                              'pophaz.alpha')


    } else if (n_z_tau >= 1 & n_z_alpha == 0) {
      names(theta) <- c('alpha0',
                        'beta',
                        'tau0',
                        paste('tau', 1:n_z_tau, colnames(z_tau), sep = "_"),
                        'pophaz.alpha')


    } else if (n_z_tau == 0 & n_z_alpha >= 1) {
      names(theta) <- c('alpha0',
                        paste('alpha', 1:n_z_alpha, colnames(z_alpha) , sep = "_"),
                        'beta',
                        'tau0',
                        'pophaz.alpha')


    }
    else if (n_z_tau == 0 & n_z_alpha == 0) {
      names(theta) <- c('alpha0','beta','tau0', 'pophaz.alpha')
    }
  }else{
    if (n_z_tau >= 1 & n_z_alpha >= 1) {
      names(theta) <- c('alpha0',
                        paste('alpha', 1:n_z_alpha, colnames(z_alpha) , sep = "_"),
                        'beta',
                        'tau0',
                        paste('tau', 1:n_z_tau, colnames(z_tau), sep = "_"))


    }else  if (n_z_tau == 0 & n_z_alpha >= 1) {
      names(theta) <- c('alpha0',
                        paste('alpha', 1:n_z_alpha, colnames(z_alpha) , sep = "_"),
                        'beta',
                        'tau0')


    } else if (n_z_tau >= 1 & n_z_alpha == 0) {
      names(theta) <- c('alpha0',
                        'beta',
                        'tau0',
                        paste('tau', 1:n_z_tau, colnames(z_tau), sep = "_"))


    } else if (n_z_tau == 0 & n_z_alpha == 0) {
      names(theta) <- c('alpha0','beta','tau0')
    }
  }




  multneh_ll_fn <- function(theta) {
    if (anyNA(theta))
      return(Inf)
    if (pophaz.alpha) {
      alpha <- matrix(1, ncol = 1, nrow = length(pophaz)) %*%
        theta[n_z_alpha + 2 + n_z_tau + 1 + 1]
    }else{
      alpha <- matrix(1, ncol = 1, nrow = length(pophaz)) %*% 0
    }
    sum_minusloklik <- -sum(d * log(exp(alpha) * pophaz +
                                      lexc_mul(z_tau, z_alpha, x, theta)) -
                              cumLexc_mul(z_tau, z_alpha, x, theta) -
                              cumpophaz * exp(alpha), na.rm = TRUE)

    if (anyNA(sum_minusloklik)) {
      return(1e10)
    } else {
      return(sum_minusloklik)
    }
  }








  gradient_fn <- function(theta) {
      FD_alpha0 <- ((d * FD_lexc_d_alpha0(z_tau, z_alpha, x, theta)) /
                    (pophaz + lexc_mul(z_tau, z_alpha, x, theta))) -
      FD_LambdaE_d_alpha0(z_tau, z_alpha, x, theta)
    FD_alphak <- ((d * FD_lexc_d_alphak(z_tau, z_alpha, x, theta)) /
                    c(pophaz + lexc_mul(z_tau, z_alpha, x, theta))) -
      FD_LambdaE_d_alphak(z_tau, z_alpha, x, theta)


    FD_beta <- ((d * FD_lexc_d_beta(z_tau, z_alpha, x, theta)) /
                  (pophaz + lexc_mul(z_tau, z_alpha, x, theta))) -
      FD_LambdaE_d_beta(z_tau, z_alpha, x, theta)

    FD_tau0 <- ((d * FD_lexc_d_tau0(z_tau, z_alpha, x, theta)) /
                  (pophaz + lexc_mul(z_tau, z_alpha, x, theta))) -
      FD_LambdaE_d_tau0(z_tau, z_alpha, x, theta)

    FD_tauk <- ((d * FD_lexc_d_tauk(z_tau, z_alpha, x, theta)) /
                  c(pophaz + lexc_mul(z_tau, z_alpha, x, theta))) -
      FD_LambdaE_d_tauk(z_tau, z_alpha, x, theta)

    FD <- c(colSums(FD_alpha0),
            colSums(FD_alphak),
            colSums(FD_beta),
            colSums(FD_tau0),
            colSums(FD_tauk))
    return(FD)
  }

  if (optim_func == "optim") {

  if (method_opt == "L-BFGS-B") {

    if (gradient) {
      NCD <- ncoor_des

      if (is.null(NCD)) {
      optimized <- stats::optim(par = theta,
                               fn = multneh_ll_fn,
                               gr = gradient_fn,
                               method = method_opt,
                               hessian = TRUE,
                               lower = theta_lower,
                               upper = theta_upper,
                               control = list(maxit = maxit_opt,
                                              trace = trace,
                                              REPORT = 500))
      } else if (is.numeric(NCD)) {
        optimized_theta <- try(stats::optim(par = theta,
                                  fn = multneh_ll_fn,
                                  gr = gradient_fn,
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
              return(multneh_ll_fn(val))
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
                                      fn = multneh_ll_fn,
                                      gr = gradient_fn,
                                      method = method_opt,
                                      hessian = TRUE,
                                      lower = theta_lower,
                                      upper = theta_upper,
                                      control = list(maxit = maxit_opt,
                                                     trace = trace,
                                                     REPORT = 500)), TRUE)



      }else if (NCD == "iter") {

        optimized_theta <- try(stats::optim(par = theta,
                                            fn  = multneh_ll_fn,
                                            gr = gradient_fn,
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
                         nrow = 1000,
                         ncol = 1)
        loglik_val <- matrix(data = 0,
                             nrow = 1000,
                             ncol = 1)
        theta_val <- matrix(data = 0,
                            nrow = 1000,
                            ncol = length(theta))
        while (abs(mydiff[j + 1, ]) > iter_eps && (j = j + 1) < 10000) {
          message("cycle = ", j, "\n")

          for (i in 1:length(init.cd)) {
            tempf <- Vectorize(function(par){
              val <- replace(init.cd, i, par)
              return(multneh_ll_fn(val))
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
              abs(loglik_val[(j + 1), ]) == .Machine$double.xmax  &
              loglik_val[(j + 1), ] > 0) {
            loglik_val[(j + 1), ] = 10000
          } else if (abs(loglik_val[(j + 1), ]) == Inf |
                     abs(loglik_val[(j + 1), ]) == .Machine$double.xmax &
                     loglik_val[(j + 1), ] < 0) {
            loglik_val[(j + 1), ] = -10000
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
                                      fn  = multneh_ll_fn,
                                      gr = gradient_fn,
                                      method = method_opt,
                                      hessian = TRUE,
                                      lower = theta_lower,
                                      upper = theta_upper,
                                      control = list(maxit = maxit_opt,
                                                     trace = trace,
                                                     REPORT = 500)), TRUE)




      }
      }
    else {
      NCD <- ncoor_des
      if (is.null(NCD)) {

        optimized <- stats::optim(par = theta,
                                 fn = multneh_ll_fn,
                                 method = method_opt,
                                 hessian = TRUE,
                                 lower = theta_lower,
                                 upper = theta_upper,
                                 control = list(maxit = maxit_opt,
                                                trace = trace,
                                                REPORT = 500))
        } else if (is.numeric(NCD)) {
          init.cd <- theta
          for (j in 1:NCD) {
            for (i in 1:length(init.cd)) {
              tempf <- Vectorize(function(par){
              val <- replace(init.cd, i, par)
              return(multneh_ll_fn(val))
            })
              new.val <- try(stats::optim(par = init.cd[i],
                               fn = tempf,
                               method = method_opt,
                               hessian = TRUE,
                               lower = theta_lower[i],
                               upper = theta_upper[i],
                               control = list(maxit = maxit_opt,
                                              trace = trace,
                                              REPORT = 500))$par, TRUE)
               message("new initial value of parameter number ",i," is equal to: ", new.val,"\n")

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

          optimized <- stats::optim(par = initG ,
                             fn = multneh_ll_fn,
                             method = method_opt,
                             hessian = TRUE,
                             lower = theta_lower,
                             upper = theta_upper,
                             control = list(maxit = maxit_opt,
                                            trace = trace,
                                            REPORT = 500))
        }else if (NCD == "iter") {

          optimized_theta <- try(stats::optim(par = theta,
                                              fn = multneh_ll_fn,
                                              gr = gradient_fn,
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
                           nrow = 1000,
                           ncol = 1)
          loglik_val <- matrix(data = 0,
                               nrow = 1000,
                               ncol = 1)
          theta_val <- matrix(data = 0,
                              nrow = 1000,
                              ncol = length(theta))
          while (abs(mydiff[j + 1, ]) > iter_eps && (j = j + 1) < 10000) {
            message("cycle = ", j, "\n")

            for (i in 1:length(init.cd)) {
              tempf <- Vectorize(function(par){
                val <- replace(init.cd, i, par)
                return(multneh_ll_fn(val))
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
            } else if (abs(loglik_val[(j + 1), ]) == Inf |
                       abs(loglik_val[(j + 1), ]) == .Machine$double.xmax &
                       loglik_val[(j + 1), ] < 0) {
              loglik_val[(j + 1), ] = -10000
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
                                        fn = multneh_ll_fn,
                                        method = method_opt,
                                        hessian = TRUE,
                                        lower = theta_lower,
                                        upper = theta_upper,
                                        control = list(maxit = maxit_opt,
                                                       trace = trace,
                                                       REPORT = 500)), TRUE)

        }
      }
    }else
    {
      if (gradient) {
        NCD <- ncoor_des
        if (is.null(NCD)) {
          optimized <- stats::optim(par = theta,
                                    fn = multneh_ll_fn,
                                    gr = gradient_fn,
                                    method = method_opt,
                                    control = list(pgtol = 1e-15,
                                                   maxit = maxit_opt,
                                                   trace = trace,
                                                   REPORT = 500))
        } else if (is.numeric(NCD)) {
          optimized_theta <- try(stats::optim(par = theta,
                                              fn = multneh_ll_fn,
                                              gr = gradient_fn,
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
                return(multneh_ll_fn(val))
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
              new.val <- init.cd[i] + 1
            }
          }

          if (any(is.na(init.cd))) {
            initG <- theta
          }  else {
            initG <- init.cd
          }

          optimized <- try(stats::optim(par = initG ,
                                        fn = multneh_ll_fn,
                                        gr = gradient_fn,
                                        method = method_opt,
                                        control = list(pgtol = 1e-15,
                                                       maxit = maxit_opt,
                                                       trace = trace,
                                                       REPORT = 500)), TRUE)




        } else if (NCD == "iter") {
          optimized_theta <- try(stats::optim(par = theta,
                                              fn = multneh_ll_fn,
                                              gr = gradient_fn,
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
                           nrow = 1000,
                           ncol = 1)
          loglik_val <- matrix(data = 0,
                               nrow = 1000,
                               ncol = 1)
          theta_val <- matrix(data = 0,
                              nrow = 1000,
                              ncol = length(theta))
          while (abs(mydiff[j + 1, ]) > iter_eps && (j = j + 1) < 10000) {
            message("cycle = ", j, "\n")

            for (i in 1:length(init.cd)) {
              tempf <- Vectorize(function(par){
                val <- replace(init.cd, i, par)
                return(multneh_ll_fn(val))
              })


              new.valpar <- try(optimx::optimr(par = init.cd[i],
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
            } else if (abs(loglik_val[(j + 1), ]) == Inf |
                       abs(loglik_val[(j + 1), ]) == .Machine$double.xmax &
                       loglik_val[(j + 1), ] < 0) {
              loglik_val[(j + 1), ] = -10000
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

          optimized <- try(optimx::optimr(par = initG ,
                                          fn = multneh_ll_fn,
                                          gr = gradient_fn,
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
                                    fn = multneh_ll_fn,
                                    method = method_opt,
                                    control = list(pgtol = 1e-15,
                                                   maxit = maxit_opt,
                                                   trace = trace,
                                                   REPORT = 500))
        } else if (is.numeric(NCD)) {
          optimized_theta <- try(stats::optim(par = theta,
                                              fn = multneh_ll_fn,
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
                return(multneh_ll_fn(val))
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
                                        fn = multneh_ll_fn,
                                        method = method_opt,
                                        control = list(pgtol = 1e-15,
                                                       maxit = maxit_opt,
                                                       trace = trace,
                                                       REPORT = 500)), TRUE)




        } else if ( NCD == "iter") {
          optimized_theta <- try(stats::optim(par = theta,
                                              fn = multneh_ll_fn,
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
                           nrow = 1000,
                           ncol = 1)
          loglik_val <- matrix(data = 0,
                               nrow = 1000,
                               ncol = 1)
          theta_val <- matrix(data = 0,
                              nrow = 1000,
                              ncol = length(theta))
          while (abs(mydiff[j + 1, ]) > iter_eps && (j = j + 1) < 10000) {
            message("cycle = ", j, "\n")

            for (i in 1:length(init.cd)) {
              tempf <- Vectorize(function(par){
                val <- replace(init.cd, i, par)
                return(multneh_ll_fn(val))
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
            } else if (abs(loglik_val[(j + 1), ]) == Inf |
                       abs(loglik_val[(j + 1), ]) == .Machine$double.xmax &
                       loglik_val[(j + 1), ] < 0) {
              loglik_val[(j + 1), ] = -10000
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
                                        fn = multneh_ll_fn,
                                        method = method_opt,
                                        hessian = TRUE,
                                        control = list(maxit = maxit_opt,
                                                       trace = trace,
                                                       REPORT = 500)), TRUE)


        }
      }}


  }

  par <- optimized$par
  AIC <- 2*optimized$value + 2*length(par)
  loglik <- - optimized$value
  evaluations <- optimized$counts
  iterations <- evaluations
  convergence <- optimized$convergence
  message <- optimized$message
   SD <- numDeriv::hessian(multneh_ll_fn, par)
    if (hessian_varcov) {
      varcov_star <- try(solve(SD), TRUE)
    } else{
      stop("optimization process requires the use of  'hessian_varcov =TRUE'")
    }

    if (any(inherits(varcov_star, "try-error"))) {
    varcov_star <- try(solve(optimized$hessian), TRUE)
    }
    std_err_star <- try(sqrt(diag(varcov_star)), TRUE)

  thetaest <- matrix(par, nrow = 1)
  n_z_tau <- ncol(z_tau)
  n_z_tau_ad <- n_z_tau - 1

  if (pophaz.alpha) {
    if (n_z_tau >= 1 & n_z_alpha >= 1) {
      colnames(thetaest) <- c('alpha0',
                              paste('alpha', 1:n_z_alpha, colnames(z_alpha) ,sep = "_"),
                              'beta',
                              'tau0',
                              paste('tau', 1:n_z_tau, colnames(z_tau), sep = "_"),
                              'pophaz.alpha')



    } else if (n_z_tau >= 1 & n_z_alpha == 0) {
      colnames(thetaest) <- c('alpha0',
                              'beta',
                              'tau0',
                              paste('tau', 1:n_z_tau, colnames(z_tau), sep = "_"),
                              'pophaz.alpha')


    } else if (n_z_tau == 0 & n_z_alpha >= 1) {
      colnames(thetaest) <- c('alpha0',
                              paste('alpha', 1:n_z_alpha, colnames(z_alpha) ,sep = "_"),
                              'beta',
                              'tau0',
                              'pophaz.alpha')


    }
      else if (n_z_tau == 0 & n_z_alpha == 0) {
      colnames(thetaest) <- c('alpha0','beta','tau0', 'pophaz.alpha')
    }
  }else{
    if (n_z_tau >= 1 & n_z_alpha >= 1) {
      colnames(thetaest) <- c('alpha0',
                              paste('alpha', 1:n_z_alpha, colnames(z_alpha) ,sep = "_"),
                              'beta',
                              'tau0',
                              paste('tau', 1:n_z_tau, colnames(z_tau), sep = "_"))


    }else  if (n_z_tau == 0 & n_z_alpha >= 1) {
      colnames(thetaest) <- c('alpha0',
                              paste('alpha', 1:n_z_alpha, colnames(z_alpha) ,sep = "_"),
                              'beta',
                              'tau0')


    } else if (n_z_tau >= 1 & n_z_alpha == 0) {
      colnames(thetaest) <- c('alpha0',
                              'beta',
                              'tau0',
                              paste('tau', 1:n_z_tau, colnames(z_tau), sep = "_"))


    } else if (n_z_tau == 0 & n_z_alpha == 0) {
      colnames(thetaest) <- c('alpha0','beta','tau0')
    }
  }

  thetaest2 <- thetaest

  id_beta <- which(colnames(thetaest) == 'beta')
  id_alpha0 <- which(colnames(thetaest) == "alpha0")
  id_tau0 <- which(colnames(thetaest) == "tau0")

  thetaest2[id_beta] <- exp(thetaest[id_beta]) + 1
  thetaest2[-id_beta] <- exp(thetaest[-id_beta])

  varcov <- try(varcov_star*c(thetaest2^2),TRUE)
  std_err <- try((sqrt(diag(varcov))), TRUE)

  if (inherits(std_err, "try-error")) {
  } else{
    std_err <- try(matrix(std_err, ncol = length(std_err)), TRUE)

  if (pophaz.alpha) {
    if (n_z_tau >= 1 & n_z_alpha >= 1) {
      colnames(std_err) <- c('alpha0',
                             paste('alpha', 1:n_z_alpha, colnames(z_alpha),
                                   sep = "_"),
                             'beta',
                             'tau0',
                             paste('tau', 1:n_z_tau, colnames(z_tau),
                                   sep = "_"),
                             'pophaz.alpha')
    } else if (n_z_tau == 0 & n_z_alpha >= 1) {
      colnames(std_err) <- c('alpha0',
                             paste('alpha', 1:n_z_alpha, colnames(z_alpha),
                                   sep = "_"),
                             'beta',
                             'tau0',
                             'pophaz.alpha')
    } else if (n_z_tau >= 1 & n_z_alpha == 0) {
      colnames(std_err) <- c('alpha0',
                             'beta',
                             'tau0',
                             paste('tau', 1:n_z_tau, colnames(z_tau),
                                   sep = "_"),
                             'pophaz.alpha')
    }
    else if (n_z_tau == 0 & n_z_alpha == 0) {
      colnames(std_err) <- c('alpha0','beta','tau0', 'pophaz.alpha')
    }
  }else{
    if (n_z_tau >= 1 & n_z_alpha >= 1) {
      colnames(std_err) <- c('alpha0',
                             paste('alpha', 1:n_z_alpha, colnames(z_alpha),
                                   sep = "_"),
                             'beta',
                             'tau0',
                             paste('tau', 1:n_z_tau, colnames(z_tau), sep = "_"))
    }else if (n_z_tau == 0 & n_z_alpha >= 1) {
      colnames(std_err) <- c('alpha0',
                             paste('alpha', 1:n_z_alpha, colnames(z_alpha),
                                   sep = "_"),
                             'beta',
                             'tau0')
    }else if (n_z_tau >= 1 & n_z_alpha == 0) {
      colnames(std_err) <- c('alpha0',
                             'beta',
                             'tau0',
                             paste('tau', 1:n_z_tau, colnames(z_tau),
                                   sep = "_"))
    }

    else if (n_z_tau >= 1 & n_z_alpha == 0) {
      colnames(std_err) <- c('alpha0','beta','tau0')
    }
  }



  colnames(thetaest) <- try(paste0(colnames(thetaest),"*"),TRUE)
  std_err_star <- matrix(std_err_star, ncol = length(std_err_star))
  colnames(std_err_star) <- try(paste0(colnames(thetaest),"*"), TRUE)
  }

  return(
    list(
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
