fit.opt.maxim_mixture_alpha <- function(x = x,
                                        d = d,
                                        z_pcured = z_pcured,
                                        z_ucured = z_ucured,
                                        theta_init = theta_init,
                                        theta_lower = theta_lower,
                                        theta_upper = theta_upper,
                                        pophaz = pophaz,
                                        pophaz.alpha = pophaz.alpha,
                                        z_pophaz.alpha = z_pophaz.alpha,
                                        z_pophaz.alpha_id = z_pophaz.alpha_id,
                                        cumpophaz = cumpophaz,
                                        method_opt = method_opt,
                                        maxit_opt = maxit_opt,
                                        optim_func = optim_func,
                                        gradient = gradient,
                                        hessian_varcov = hessian_varcov,
                                        model = model,
                                        ncoor_des = ncoor_des,
                                        trace = trace,
                                        optim_fixed = optim_fixed,
                                        optimizer = optimizer,
                                        sign_delta = sign_delta) {

  theta  <- theta_init
  n_z_pcured <- ncol(z_pcured)
  n_z_ucured <- ncol(z_ucured)

  fixed <- optim_fixed

  if (n_z_pcured > 0 & n_z_ucured > 0) {
    if (pophaz.alpha) {
      if (length(z_pophaz.alpha_id) > 0) {
        names(theta) <- c(
          'beta0',
          paste('beta', 1:n_z_pcured, colnames(z_pcured) , sep = "_"),
          'lambda', 'gamma',
          paste('delta', 1:n_z_ucured, colnames(z_ucured), sep = "_"),
          'pophaz.alpha0', paste('pophaz.alpha', 1:length(z_pophaz.alpha_id), colnames(z_pophaz.alpha) , sep = "_"))
      }else{
        names(theta) <- c(
          'beta0',
          paste('beta', 1:n_z_pcured, colnames(z_pcured) , sep = "_"),
          'lambda', 'gamma',
          paste('delta', 1:n_z_ucured, colnames(z_ucured), sep = "_"),
          'pophaz.alpha0')
      }


    } else{
      names(theta) <- c(
        'beta0',
        paste('beta', 1:n_z_pcured, colnames(z_pcured) , sep = "_"),
        'lambda', 'gamma',
        paste('delta', 1:n_z_ucured, colnames(z_ucured), sep = "_"))
    }



  } else if (n_z_pcured == 0 & n_z_ucured > 0) {
    if (pophaz.alpha) {
      if (length(z_pophaz.alpha_id) > 0) {
        names(theta) <- c(
          'beta0',
          'lambda', 'gamma',
          paste('delta', 1:n_z_ucured, colnames(z_ucured), sep = "_"),
          'pophaz.alpha0', paste('pophaz.alpha', 1:length(z_pophaz.alpha_id), colnames(z_pophaz.alpha) , sep = "_"))
      }else{
        names(theta) <- c(
          'beta0',
          'lambda', 'gamma',
          paste('delta', 1:n_z_ucured, colnames(z_ucured), sep = "_"),
          'pophaz.alpha0')
      }


    } else{
      names(theta) <- c(
        'beta0',
        'lambda', 'gamma',
        paste('delta', 1:n_z_ucured, colnames(z_ucured), sep = "_"))
    }



  } else if (n_z_pcured > 0 & n_z_ucured == 0) {
    if (pophaz.alpha) {
      if (length(z_pophaz.alpha_id) > 0) {
        names(theta) <- c(
          'beta0',
          paste('beta', 1:n_z_pcured, colnames(z_pcured) , sep = "_"),
          'lambda', 'gamma',
          'pophaz.alpha0', paste('pophaz.alpha', 1:length(z_pophaz.alpha_id), colnames(z_pophaz.alpha) , sep = "_"))
      }else{
        names(theta) <- c(
          'beta0',
          paste('beta', 1:n_z_pcured, colnames(z_pcured) , sep = "_"),
          'lambda', 'gamma',
          'pophaz.alpha0')
      }


    } else{
      names(theta) <- c(
        'beta0',
        paste('beta', 1:n_z_pcured, colnames(z_pcured) , sep = "_"),
        'lambda', 'gamma')
    }



  }  else if (n_z_pcured == 0 & n_z_ucured == 0) {

    if (pophaz.alpha) {

      if (length(z_pophaz.alpha_id) > 0) {
        names(theta) <- c('beta0',
                          'lambda',
                          'gamma',
                          'pophaz.alpha0',
                          paste('pophaz.alpha',
                                1:length(z_pophaz.alpha_id),
                                colnames(z_pophaz.alpha) ,
                                sep = "_"))
      } else{
        names(theta) <- c('beta0', 'lambda', 'gamma', 'pophaz.alpha0')
      }
    } else{
      names(theta) <- c('beta0', 'lambda', 'gamma')
    }
  }


  theta2 <- theta

  if (gradient == TRUE) {

    if (n_z_pcured > 0 & n_z_ucured > 0) {




      ll_fn_alphaweibull_deriv_i <-  function(beta0, betak, lambda, gamma, delta, x, d, pophaz,
                                              cumpophaz,
                                              alpha_pophaz,
                                              z_pcured, z_ucured,
                                              n_z_pcured, n_z_ucured)
      {
        .e1 <- exp(gamma)
        .e3 <- exp(-(beta0 + z_pcured %*% betak))
        .e4 <- 1 + .e3
        .e5 <- 1/.e4
        .e6 <- exp(lambda + z_ucured %*% delta)
        .e7 <- x^.e1
        .e9 <- exp(-(.e7 * .e6))
        .e10 <- 1 - .e5
        .e12 <- .e10 * .e9 + .e5
        .e14 <- x^(.e1 - 1)
        .e15 <- .e14 * .e10
        .e21 <- pophaz * exp(alpha_pophaz) + .e15 * .e9 * .e1 *
          .e6/.e12
        .e22 <- d * .e1
        .e24 <- .e6 * (x^(2 * .e1 - 1) * .e10 * .e9/.e12 - x^(2 *
                                                                .e1 - 1)) + .e14
        .e25 <- .e12 * .e4^2
        .e27 <- 1 - (1 + .e22 * .e6 * (.e14 + .e15 * (1 - .e9)/.e12)/.e21) *
          .e9
        .e30 <- .e22 * .e24/.e21 - .e7
        .e31 <- log(x)
        cbind(beta0 = -(.e27 * .e3/.e25),
              betak = -(z_pcured * c(.e27 * .e3/.e25)),
              lambda = -(.e10 * .e30 * .e9 * .e6/.e12),
              gamma = -(.e10 * (d * (.e1 * .e24 * .e31 + .e14)/.e21 -
                                  .e7 * .e31) * .e9 * .e1 * .e6/.e12),
              delta = -(z_ucured * c(.e10 * .e30 * .e9 * .e6/.e12)))
      }


      ll_fn_alphaweibull_deriv_all <- function(theta
      ) {
        beta0  <- theta[1] ; lambda  <- theta[(1 + n_z_pcured + 1)]
        gamma <- theta[(1 + n_z_pcured + 2)] ;
        betak <- theta[2:(1 + n_z_pcured)]
        delta <- theta[(1 + n_z_pcured + 3):(1 + n_z_pcured + 2 + n_z_ucured)]
        if (n_z_pcured >1) {
          names(betak) <- c(paste0("betak", 1:n_z_pcured))
        }
        if (n_z_ucured>1) {
          names(delta) <- c(paste0("delta", 1:n_z_ucured))

        }

        alpha_pophaz <- c(matrix(1, ncol = 1, nrow = length(pophaz)) %*% 0)
        pophaz <- c(pophaz)
        cumpophaz <- c(cumpophaz)

        gradval <- colSums(ll_fn_alphaweibull_deriv_i(beta0, betak, lambda,
                                                      gamma, delta,
                                                      x, d, pophaz,
                                                      cumpophaz,
                                                      alpha_pophaz,
                                                      z_pcured, z_ucured,
                                                      n_z_pcured, n_z_ucured))
        names(gradval) <- NULL
        return(gradval)
      }




    }
    else if (n_z_pcured > 0 & n_z_ucured == 0) {

      ll_fn_alphaweibull_deriv_i_z_p <- function(beta0, betak, lambda, gamma,
                                                 x, d, pophaz,
                                                 cumpophaz,
                                                 alpha_pophaz,
                                                 z_pcured)
      {
        .e1 <- exp(gamma)
        .e3 <- exp(-(beta0 + z_pcured %*% betak))
        .e4 <- 1 + .e3
        .e5 <- 1/.e4
        .e6 <- exp(lambda)
        .e7 <- x^.e1
        .e9 <- exp(-(.e7 * .e6))
        .e10 <- 1 - .e5
        .e12 <- .e10 * .e9 + .e5
        .e14 <- x^(.e1 - 1)
        .e15 <- .e14 * .e10
        .e20 <- pophaz * exp(alpha_pophaz) + .e15 * .e9 * .e1 *
          .e6/.e12
        .e22 <- d * .e1
        .e23 <- .e12 * .e4^2
        .e25 <- 1 - (1 + .e22 * .e6 * (.e14 + .e15 * (1 - .e9)/.e12)/.e20) *
          .e9
        .e27 <- .e6 * (x^(2 * .e1 - 1) * .e10 * .e9/.e12 - x^(2 *
                                                                .e1 - 1)) + .e14
        .e28 <- log(x)
        cbind(beta0 = -(.e25 * .e3/.e23),
              betak = -(z_pcured * c(.e25 * .e3/.e23)),
              lambda = -(.e10 * (.e22 * .e27/.e20 - .e7) * .e9 * .e6/.e12),
              gamma = -(.e10 * (d * (.e1 * .e27 * .e28 + .e14)/.e20 - .e7 * .e28) * .e9 * .e1 * .e6/.e12))
      }





      ll_fn_alphaweibull_deriv_all <- function(theta
      ) {
        beta0  <- theta[1] ; lambda  <- theta[(1 + n_z_pcured + 1)]
        gamma <- theta[(1 + n_z_pcured + 2)] ;
        betak <- theta[2:(1 + n_z_pcured)]
        alpha_pophaz <- matrix(1, ncol = 1, nrow = length(pophaz)) %*% 0

        gradval <- colSums(ll_fn_alphaweibull_deriv_i_z_p(beta0, betak, lambda, gamma,
                                                          x, d, pophaz,
                                                          cumpophaz,
                                                          alpha_pophaz,
                                                          z_pcured))
        names(gradval) <- NULL
        return(gradval)
      }




    }else if (n_z_pcured == 0 & n_z_ucured > 0) {

      ll_fn_alphaweibull_deriv_i_z_u <- function(beta0, lambda, gamma,
                                                 delta, x, d, pophaz,
                                                 cumpophaz,
                                                 alpha_pophaz,
                                                 z_ucured) {
        .e1 <- exp(gamma)
        .e3 <- exp(-beta0)
        .e4 <- 1 + .e3
        .e5 <- 1/.e4
        .e6 <- exp(lambda + z_ucured %*% delta)
        .e7 <- x^.e1
        .e9 <- exp(-(.e7 * .e6))
        .e10 <- 1 - .e5
        .e12 <- .e10 * .e9 + .e5
        .e14 <- x^(.e1 - 1)
        .e18 <- .e14 * .e10
        .e21 <- pophaz * exp(alpha_pophaz) + .e18 * .e9 * .e1 *
          .e6/.e12
        .e22 <- d * .e1
        .e24 <- .e6 * (x^(2 * .e1 - 1) * .e10 * .e9/.e12 - x^(2 *
                                                                .e1 - 1)) + .e14
        .e27 <- .e22 * .e24/.e21 - .e7
        .e28 <- log(x)
        c(beta0 = -((1 - (1 + .e22 * .e6 * (.e14 + .e18 * (1 - .e9)/.e12)/.e21) *
                       .e9) * .e3/(.e12 * .e4^2)),
          lambda = -(.e10 * .e27 * .e9 * .e6/.e12),
          gamma = -(.e10 * (d * (.e1 * .e24 * .e28 + .e14)/.e21 - .e7 * .e28) * .e9 * .e1 * .e6/.e12),
          delta = -(z_ucured * c(.e10 * .e27 * .e9 * .e6/.e12)))
      }


      ll_fn_alphaweibull_deriv_all <- function(theta
      ) {
        beta0  <- theta[1] ; lambda  <- theta[(1 + n_z_pcured + 1)]
        gamma <- theta[(1 + n_z_pcured + 2)]
        delta <- theta[(1 + n_z_pcured + 3):(1 + n_z_pcured + 2 + n_z_ucured)]
        alpha_pophaz <- matrix(1, ncol = 1, nrow = length(pophaz)) %*% 0

        gradval <- colSums(attr(ll_fn_alphaweibull_deriv_i_z_u(beta0, lambda, gamma,
                                                               delta, x, d, pophaz,
                                                               cumpophaz,
                                                               alpha_pophaz,
                                                               z_ucured),
                                "gradient"))
        names(gradval) <- NULL
        return(gradval)
      }



    }else if (n_z_pcured == 0 & n_z_ucured == 0) {
      ll_fn_alphaweibull_deriv_i <- deriv(
        ~ -(
          d * log((exp(alpha_pophaz)) * pophaz +
                    (((1 - (1/(1 + exp(-beta0))))*
                        (exp(gamma)*exp(lambda)*((x)^(exp(gamma) - 1))) *
                        (exp(-exp(lambda)*(x)^exp(gamma)))) /
                       ((1/(1 + exp(-beta0))) +
                          (1 - (1/(1 + exp(-beta0))))*(exp(-exp(lambda)*(x)^exp(gamma)))))) -
            -log((1/(1 + exp(-beta0))) + (1 - (1/(1 + exp(-beta0))))*
                   (exp(-exp(lambda)*(x)^exp(gamma)))) - cumpophaz * (exp(alpha_pophaz))),
        c("beta0", "lambda", "gamma"), function(beta0, lambda, gamma,
                                                x, d, pophaz, cumpophaz,
                                                alpha_pophaz
        )NULL)


      ll_fn_alphaweibull_deriv_all <- function(theta
      ) {
        beta0  <- theta[1] ; lambda  <- theta[2] ; gamma <- theta[3]
        alpha_pophaz <- matrix(1, ncol = 1, nrow = length(pophaz)) %*% 0

        gradval <- colSums(attr(ll_fn_alphaweibull_deriv_i(beta0, lambda, gamma,
                                                           x, d, pophaz, cumpophaz,
                                                           alpha_pophaz
        ),
        "gradient"))
        return(gradval)
      }

    }

  }

  ll_fn_alphaweibull <- function(theta) {

    if (pophaz.alpha) {
      if (length(z_pophaz.alpha_id) > 0) {
        alpha_pophaz <- matrix(1, ncol = 1,
                               nrow = length(pophaz)) %*%
          theta["pophaz.alpha0"] + z_pophaz.alpha %*% theta["pophaz.alpha_k"]
      } else {
        alpha_pophaz <- matrix(1,
                               ncol = 1,
                               nrow = length(pophaz)) %*% theta["pophaz.alpha0"]
      }

    }else{
      alpha_pophaz <- matrix(1, ncol = 1, nrow = length(pophaz)) %*% 0

    }

    if (anyNA(theta))
      return(Inf)

    loglikval <- c(sum(
      d * log(c(exp(alpha_pophaz)) * pophaz + lexc_alphaweibull(z_ucured, z_pcured, x, theta)) -
        cumLexc_alphaweibull(z_ucured, z_pcured, x, theta) - cumpophaz * c(exp(alpha_pophaz))))
    return(-loglikval)
  }




  ll_fn_alphaweibull3 <- function(beta0,
                                  beta_k,
                                  lambda,
                                  gamma,
                                  delta_k,
                                  pophaz.alpha0 = NULL,
                                  pophaz.alpha_k = NULL) {

    if (pophaz.alpha) {
      if (n_z_pcured > 0) {
        if (optim_func == "bbmle") {

          beta_k <- eval(
            parse(text = paste("c(",
                               paste(paste("`",paste('beta', 1:n_z_pcured,  colnames(z_pcured) , sep = "_"),"`", sep = ""),
                                     collapse = ","),
                               ")", sep = "")))

          delta_k <- eval(
            parse(text = paste("c(",
                               paste( paste("`",paste('delta', 1:n_z_ucured,  colnames(z_ucured), sep = "_"), "`", sep = ""),
                                      collapse = ","),
                               ")", sep = "")))
        }


        if (length(z_pophaz.alpha_id) > 0) {
          theta <- c(beta0, beta_k, lambda, gamma, delta_k, pophaz.alpha0, pophaz.alpha_k)
          alpha_pophaz <- matrix(1, ncol = 1, nrow = length(pophaz)) %*% pophaz.alpha0 + z_pophaz.alpha %*% pophaz.alpha_k
        } else {
          theta <- c(beta0, beta_k, lambda, gamma, delta_k, pophaz.alpha0)
          alpha_pophaz <- matrix(1, ncol = 1, nrow = length(pophaz)) %*% pophaz.alpha0
        }

      }else {
        if (length(z_pophaz.alpha_id) > 0) {
          theta <- c(beta0,  lambda, gamma,  pophaz.alpha0, pophaz.alpha_k)
          alpha_pophaz <- matrix(1, ncol = 1, nrow = length(pophaz)) %*% pophaz.alpha0 + z_pophaz.alpha %*% pophaz.alpha_k

        }else{
          theta <- c(beta0,  lambda, gamma,  pophaz.alpha0, pophaz.alpha_k)
          alpha_pophaz <- matrix(1, ncol = 1, nrow = length(pophaz)) %*% pophaz.alpha0

        }
      }

    }else{


      if (n_z_pcured > 0) {
        if (optim_func == "bbmle") {

          beta_k <- eval(
            parse(text = paste("c(",
                               paste(paste("`", paste('beta', 1:n_z_pcured,  colnames(z_pcured) , sep = "_"), "`",
                                           sep = ""),
                                     collapse = ","),
                               ")", sep = "")))

          delta_k <- eval(
            parse(text = paste("c(",
                               paste( paste("`", paste('delta', 1:n_z_ucured,  colnames(z_ucured), sep = "_"), "`",
                                            sep = ""),
                                      collapse = ","),
                               ")", sep = "")))
        }
        theta <- c(beta0, beta_k, lambda, gamma, delta_k)
        alpha_pophaz <- matrix(1, ncol = 1, nrow = length(pophaz)) %*% 0
      }else {
        theta <- c(beta0, lambda, gamma)
        alpha_pophaz <- matrix(1, ncol = 1, nrow = length(pophaz)) %*% 0
      }
    }

    if (anyNA(theta))
      return(1e10)


    loglikval <- c(-sum(
      d * log(c(exp(alpha_pophaz)) * pophaz + lexc_alphaweibull(z_ucured, z_pcured, x, theta)) -
        cumLexc_alphaweibull(z_ucured, z_pcured, x, theta) - cumpophaz * c(exp(alpha_pophaz))))



    if (anyNA(loglikval)) {
      return(1e10)
    } else {
      return(loglikval)
    }

  }

  if (optim_func == "optim") {
    if (method_opt == "L-BFGS-B") {

      if (gradient) {


        NCD <- ncoor_des
        if (is.null(NCD)) {
          optimized <- stats::optim(par = theta,
                                    fn = ll_fn_alphaweibull,
                                    gr = ll_fn_alphaweibull_deriv_all,
                                    method = method_opt,
                                    hessian = TRUE,
                                    lower = theta_lower,
                                    upper = theta_upper,
                                    control = list(maxit = maxit_opt)
          )
        } else {
          init.cd <- theta
          for (j in 1:NCD) {
            for (i in 1:length(init.cd)) {
              tempf <- Vectorize(function(par){
                val <- replace(init.cd, i, par)
                return(ll_fn_alphaweibull(val))
              })


              tempgr <- Vectorize(function(par){
                val <- replace(init.cd, i, par)
                return(ll_fn_alphaweibull_deriv_all(val)[i])
              })


              new.val <- stats::optim(par = init.cd[i],
                                      fn = tempf,
                                      gr = tempgr,
                                      method = method_opt,
                                      hessian = TRUE,
                                      lower = theta_lower[i],
                                      upper = theta_upper[i],
                                      control = list(maxit = maxit_opt,
                                                     trace = trace,
                                                     REPORT = 500))$par
              message("new initial value of parameter number ",i," is equal to: ",new.val,"\n")
              init.cd <- replace(init.cd, i,  new.val )
            }
          }

          if (any(is.na(init.cd))) {
            initG <- theta
          }  else {
            initG <- init.cd
          }

          optimized <- stats::optim(par = initG ,
                                    fn = ll_fn_alphaweibull,
                                    gr = ll_fn_alphaweibull_deriv_all,
                                    method = method_opt,
                                    hessian = TRUE,
                                    lower = theta_lower,
                                    upper = theta_upper,
                                    control = list(maxit = maxit_opt,
                                                   trace = trace,
                                                   REPORT = 500))
        }

      }else {

        NCD <- ncoor_des
        if (is.null(NCD)) {

          optimized <- stats::optim(par = theta,
                                    fn = ll_fn_alphaweibull,
                                    method = method_opt)
        } else {
          init.cd <- theta
          for (j in 1:NCD) {
            for (i in 1:length(init.cd)) {
              tempf <- Vectorize(function(par){
                val <- replace(init.cd, i, par)
                return(ll_fn_alphaweibull(val))
              })
              new.val <- stats::optim(par = init.cd[i],
                                      fn = tempf,
                                      method = method_opt,
                                      hessian = TRUE,
                                      lower = theta_lower[i],
                                      upper = theta_upper[i],
                                      control = list(maxit = maxit_opt,
                                                     trace = trace,
                                                     REPORT = 500))$par
              message("new initial value of parameter number ",i," is equal to: ",new.val,"\n")
              init.cd <- replace(init.cd, i,  new.val )
            }
          }

          if (any(is.na(init.cd))) {
            initG <- theta
          }  else {
            initG <- init.cd
          }

          optimized <- stats::optim(par = initG ,
                                    fn = ll_fn_alphaweibull,
                                    method = method_opt,
                                    hessian = TRUE,
                                    lower = theta_lower,
                                    upper = theta_upper,
                                    control = list(maxit = maxit_opt,
                                                   trace = trace,
                                                   REPORT = 500))
        }

      }
    } else {
      optimized <- stats::optim(par = theta,
                                fn = ll_fn_alphaweibull,
                                method = method_opt,
                                control = list(pgtol   = 1e-15,
                                               maxit   = maxit_opt)
                                )

    }

  }else if (optim_func == "bbmle") {

    if (method_opt == "L-BFGS-B" |
        method_opt == "Rvmmin" |
        method_opt == "Rcgmin" |
        method_opt == "nlminb" |
        method_opt == "bobyqa" ) {
      if (gradient) {
        stop("not yet implemented with gradient")
      }else {
        NCD <- ncoor_des
        if (pophaz.alpha) {
          alpha0 <- theta[length(theta)]

        }else {
          alpha0 <- 1
        }

        if (is.null(NCD)) {
          par_alphaweibull <- genpar_alphaweibull(theta, z_pcured,
                                                  pophaz.alpha,
                                                  z_pophaz.alpha,
                                                  z_pophaz.alpha_id)

          if (pophaz.alpha == TRUE & n_z_pcured > 0) {
            beta0 <- par_alphaweibull["beta0"]
            beta_k <- par_alphaweibull[
              c(which(stringr::str_detect(names(par_alphaweibull),"beta_k")))]
            lambda <- par_alphaweibull["lambda"]
            gamma <- par_alphaweibull["gamma"]
            delta_k <- par_alphaweibull[
              c(which(stringr::str_detect(names(par_alphaweibull),"delta_k")))]
            pophaz.alpha0 <- par_alphaweibull["pophaz.alpha0"]
            if (length(z_pophaz.alpha_id) > 0) {
              pophaz.alpha_k <- par_alphaweibull[
                c(which(stringr::str_detect(
                  names(par_alphaweibull),"pophaz.alpha_k")))]

              start <- c(beta0 = as.numeric(beta0),
                         beta_k = as.numeric(beta_k),
                         lambda = as.numeric(lambda),
                         gamma = as.numeric(gamma),
                         delta_k = as.numeric(delta_k),
                         pophaz.alpha0 = as.numeric(pophaz.alpha0),
                         pophaz.alpha_k = as.numeric(pophaz.alpha_k))
            }else{
              start <- c(beta0 = as.numeric(beta0),
                         beta_k = as.numeric(beta_k),
                         lambda = as.numeric(lambda),
                         gamma = as.numeric(gamma),
                         delta_k = as.numeric(delta_k),
                         pophaz.alpha0 = as.numeric(pophaz.alpha0))
            }



            names(start) <- names(theta2)



          }else if (pophaz.alpha == TRUE & n_z_pcured <= 0) {
            beta0 <- par_alphaweibull["beta0"]
            lambda <- par_alphaweibull["lambda"]
            gamma <- par_alphaweibull["gamma"]
            pophaz.alpha0 <- par_alphaweibull["pophaz.alpha0"]
            if (length(z_pophaz.alpha_id) > 0) {
              pophaz.alpha_k <- par_alphaweibull[
                c(which(stringr::str_detect(
                  names(par_alphaweibull),"pophaz.alpha_k")))]

              start <- c(beta0 = as.numeric(beta0),
                         lambda = as.numeric(lambda),
                         gamma = as.numeric(gamma),
                         pophaz.alpha0 = as.numeric(pophaz.alpha0),
                         pophaz.alpha_k = as.numeric(pophaz.alpha_k))
            }else{
              start <- c(beta0 = as.numeric(beta0),
                         lambda = as.numeric(lambda),
                         gamma = as.numeric(gamma),
                         pophaz.alpha0 = as.numeric(pophaz.alpha0))
            }



            names(start) <- names(theta2)



          }else if (pophaz.alpha == FALSE & n_z_pcured > 0) {
            beta0 <- par_alphaweibull["beta0"]
            beta_k <- par_alphaweibull[
              c(which(stringr::str_detect(names(par_alphaweibull),"beta_k")))]
            lambda <- par_alphaweibull["lambda"]
            gamma <- par_alphaweibull["gamma"]
            delta_k <- par_alphaweibull[
              c(which(stringr::str_detect(names(par_alphaweibull),"delta_k")))]

            start <- c(beta0 = as.numeric(beta0),
                       beta_k = as.numeric(beta_k),
                       lambda = as.numeric(lambda),
                       gamma = as.numeric(gamma),
                       delta_k = as.numeric(delta_k))

            names(start) <- names(theta2)



          }else if (pophaz.alpha == FALSE & n_z_pcured <= 0) {
            beta0 <- par_alphaweibull["beta0"]
            lambda <- par_alphaweibull["lambda"]
            gamma <- par_alphaweibull["gamma"]

            start <- c(beta0 = as.numeric(beta0),
                       lambda = as.numeric(lambda),
                       gamma = as.numeric(gamma))

            names(start) <- names(theta2)
          }


          add.arguments <- function(f, params) {
            t <- paste("arg <- alist(",
                       paste(sapply(1:length(params),
                                    function(i)
                                      paste(paste('`',names(params),"`", sep = "")[i], "=", sep = "")),
                             collapse = ","),
                       ")", sep = "")
            formals(f) <- eval(parse(text = t))
            f
          }

          ll_fn_alphaweibull3noalpha2 <- add.arguments(ll_fn_alphaweibull3, start)
          bbmle::parnames(ll_fn_alphaweibull3noalpha2) <- names(theta2)


          start2 <- paste("list(",
                          paste(sapply(1:length(start),
                                       function(i)
                                         paste(names(start)[i], "=",
                                               (start)[i], sep = "")),
                                collapse = ","),
                          ")", sep = "")
          start3 <- eval(parse(text = start2))

          theta_lower2 <- paste("list(",
                                paste(sapply(1:length(start),
                                             function(i)
                                               paste(names(start)[i], "=",
                                                     (theta_lower)[i],
                                                     sep = "")),
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
            minuslogl = ll_fn_alphaweibull3noalpha2,
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
                return(ll_fn_alphaweibull(val))
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
                                        fn = ll_fn_alphaweibull(),
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
      alpha0 <- theta[length(theta)]

    }else {
      alpha0 <- 1
    }

      if (method_opt == "all.methods") {

        NCD <- ncoor_des

        if (is.null(NCD)) {
          par_alphaweibull <- genpar_alphaweibull(theta, z_pcured,
                                                  pophaz.alpha,
                                                  z_pophaz.alpha,
                                                  z_pophaz.alpha_id)

          if (pophaz.alpha == TRUE & n_z_pcured > 0) {
            beta0 <- par_alphaweibull["beta0"]
            beta_k <- par_alphaweibull[
              c(which(stringr::str_detect(names(par_alphaweibull),"beta_k")))]
            lambda <- par_alphaweibull["lambda"]
            gamma <- par_alphaweibull["gamma"]
            delta_k <- par_alphaweibull[
              c(which(stringr::str_detect(names(par_alphaweibull),"delta_k")))]
            pophaz.alpha0 <- par_alphaweibull["pophaz.alpha0"]
            if (length(z_pophaz.alpha_id) > 0) {
              pophaz.alpha_k <- par_alphaweibull[
                c(which(stringr::str_detect(
                  names(par_alphaweibull),"pophaz.alpha_k")))]

              start <- c(beta0 = as.numeric(beta0),
                         beta_k = as.numeric(beta_k),
                         lambda = as.numeric(lambda),
                         gamma = as.numeric(gamma),
                         delta_k = as.numeric(delta_k),
                         pophaz.alpha0 = as.numeric(pophaz.alpha0),
                         pophaz.alpha_k = as.numeric(pophaz.alpha_k))
            }else{
              start <- c(beta0 = as.numeric(beta0),
                         beta_k = as.numeric(beta_k),
                         lambda = as.numeric(lambda),
                         gamma = as.numeric(gamma),
                         delta_k = as.numeric(delta_k),
                         pophaz.alpha0 = as.numeric(pophaz.alpha0))
            }



            names(start) <- names(theta2)



          }else if (pophaz.alpha == TRUE & n_z_pcured <= 0) {
            beta0 <- par_alphaweibull["beta0"]
            lambda <- par_alphaweibull["lambda"]
            gamma <- par_alphaweibull["gamma"]
            pophaz.alpha0 <- par_alphaweibull["pophaz.alpha0"]
            if (length(z_pophaz.alpha_id) > 0) {
              pophaz.alpha_k <- par_alphaweibull[
                c(which(stringr::str_detect(
                  names(par_alphaweibull),"pophaz.alpha_k")))]

              start <- c(beta0 = as.numeric(beta0),
                         lambda = as.numeric(lambda),
                         gamma = as.numeric(gamma),
                         pophaz.alpha0 = as.numeric(pophaz.alpha0),
                         pophaz.alpha_k = as.numeric(pophaz.alpha_k))
            }else{
              start <- c(beta0 = as.numeric(beta0),
                         lambda = as.numeric(lambda),
                         gamma = as.numeric(gamma),
                         pophaz.alpha0 = as.numeric(pophaz.alpha0))
            }



            names(start) <- names(theta2)



          }else if (pophaz.alpha == FALSE & n_z_pcured > 0) {
            beta0 <- par_alphaweibull["beta0"]
            beta_k <- par_alphaweibull[
              c(which(stringr::str_detect(names(par_alphaweibull),"beta_k")))]
            lambda <- par_alphaweibull["lambda"]
            gamma <- par_alphaweibull["gamma"]
            delta_k <- par_alphaweibull[
              c(which(stringr::str_detect(names(par_alphaweibull),"delta_k")))]

            start <- c(beta0 = as.numeric(beta0),
                       beta_k = as.numeric(beta_k),
                       lambda = as.numeric(lambda),
                       gamma = as.numeric(gamma),
                       delta_k = as.numeric(delta_k))

            names(start) <- names(theta2)



          }else if (pophaz.alpha == FALSE & n_z_pcured <= 0) {
            beta0 <- par_alphaweibull["beta0"]
            lambda <- par_alphaweibull["lambda"]
            gamma <- par_alphaweibull["gamma"]

            start <- c(beta0 = as.numeric(beta0),
                       lambda = as.numeric(lambda),
                       gamma = as.numeric(gamma))

            names(start) <- names(theta2)
          }


          add.arguments <- function(f, params) {
            t <- paste("arg <- alist(",
                       paste(sapply(1:length(params),
                                    function(i)
                                      paste(paste('`',names(params),"`", sep = "")[i], "=", sep = "")),
                             collapse = ","),
                       ")", sep = "")
            formals(f) <- eval(parse(text = t))
            f
          }

          ll_fn_alphaweibull3noalpha2 <- add.arguments(ll_fn_alphaweibull3, start)
          bbmle::parnames(ll_fn_alphaweibull3noalpha2) <- names(theta2)

          start2 <- paste("list(",
                          paste(sapply(1:length(start),
                                       function(i)
                                         paste(paste('`',names(start),"`", sep = "")[i], "=",(start)[i], sep = "")),
                                collapse = ","),
                          ")", sep = "")
          start3 <- eval(parse(text = start2))

          theta_lower2 <- paste("list(",
                                paste(sapply(1:length(start),
                                             function(i)
                                               paste(paste('`',names(start),"`", sep = "")[i], "=",
                                                     (theta_lower)[i], sep = "")),
                                      collapse = ","),
                                ")", sep = "")

          theta_upper2 <- paste("list(",
                                paste(sapply(1:length(start),
                                             function(i)
                                               paste(paste('`',names(start),"`", sep = "")[i], "=",
                                                     (theta_upper)[i], sep = "")),
                                      collapse = ","),
                                ")", sep = "")

          theta_lower3 <- eval(parse(text = theta_lower2))
          theta_upper3 <- eval(parse(text = theta_upper2))




          optimized2 <- optimx::optimx(
            par = unlist(start3),
            fn = ll_fn_alphaweibull,
            lower = unlist(theta_lower3),
            upper = unlist(theta_upper3),
            control = list(all.methods = TRUE,
                           maximize = FALSE,
                           trace = trace,
                           save.failures = TRUE))

          id_conver <- which(optimized2$convcode == 0 & optimized2$kkt1 == TRUE & optimized2$kkt2 == TRUE)
          id_minval <- which.min(optimized2$value[id_conver])
          start4 <- optimized2[id_conver[id_minval],1:attributes(optimized2)$npar]
          if (is.data.frame(start4) & nrow(start4)>0) {
            start3 <- c(start4)
          }

          optimized <- bbmle::mle2(
            start = start3,
            minuslogl = ll_fn_alphaweibull3noalpha2,
            lower = unlist(theta_lower3),
            upper = unlist(theta_upper3))


        } else {
          init.cd <- theta
          for (j in 1:NCD) {
            message("cycle = ", j, "\n")

            for (i in 1:length(init.cd)) {
              tempf <- Vectorize(function(par){
                val <- replace(init.cd, i, par)
                return(ll_fn_alphaweibull(val))
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
                                        fn = ll_fn_alphaweibull(),
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
          stop("not implemented for gradient")
        }else{
          par_alphaweibull <- genpar_alphaweibull(theta, z_pcured,
                                                  pophaz.alpha,
                                                  z_pophaz.alpha,
                                                  z_pophaz.alpha_id)

          if (pophaz.alpha == TRUE & n_z_pcured > 0) {
            beta0 <- par_alphaweibull["beta0"]
            beta_k <- par_alphaweibull[
              c(which(stringr::str_detect(names(par_alphaweibull),"beta_k")))]
            lambda <- par_alphaweibull["lambda"]
            gamma <- par_alphaweibull["gamma"]
            delta_k <- par_alphaweibull[
              c(which(stringr::str_detect(names(par_alphaweibull),"delta_k")))]
            pophaz.alpha0 <- par_alphaweibull["pophaz.alpha0"]
            if (length(z_pophaz.alpha_id) > 0) {
              pophaz.alpha_k <- par_alphaweibull[
                c(which(stringr::str_detect(
                  names(par_alphaweibull),"pophaz.alpha_k")))]

              start <- c(beta0 = as.numeric(beta0),
                         beta_k = as.numeric(beta_k),
                         lambda = as.numeric(lambda),
                         gamma = as.numeric(gamma),
                         delta_k = as.numeric(delta_k),
                         pophaz.alpha0 = as.numeric(pophaz.alpha0),
                         pophaz.alpha_k = as.numeric(pophaz.alpha_k))
            }else{
              start <- c(beta0 = as.numeric(beta0),
                         beta_k = as.numeric(beta_k),
                         lambda = as.numeric(lambda),
                         gamma = as.numeric(gamma),
                         delta_k = as.numeric(delta_k),
                         pophaz.alpha0 = as.numeric(pophaz.alpha0))
            }



            names(start) <- names(theta2)



          }else if (pophaz.alpha == TRUE & n_z_pcured <= 0) {
            beta0 <- par_alphaweibull["beta0"]
            lambda <- par_alphaweibull["lambda"]
            gamma <- par_alphaweibull["gamma"]
            pophaz.alpha0 <- par_alphaweibull["pophaz.alpha0"]
            if (length(z_pophaz.alpha_id) > 0) {
              pophaz.alpha_k <- par_alphaweibull[
                c(which(stringr::str_detect(
                  names(par_alphaweibull),"pophaz.alpha_k")))]

              start <- c(beta0 = as.numeric(beta0),
                         lambda = as.numeric(lambda),
                         gamma = as.numeric(gamma),
                         pophaz.alpha0 = as.numeric(pophaz.alpha0),
                         pophaz.alpha_k = as.numeric(pophaz.alpha_k))
            }else{
              start <- c(beta0 = as.numeric(beta0),
                         lambda = as.numeric(lambda),
                         gamma = as.numeric(gamma),
                         pophaz.alpha0 = as.numeric(pophaz.alpha0))
            }



            names(start) <- names(theta2)



          }else if (pophaz.alpha == FALSE & n_z_pcured > 0) {
            beta0 <- par_alphaweibull["beta0"]
            beta_k <- par_alphaweibull[
              c(which(stringr::str_detect(names(par_alphaweibull),"beta_k")))]
            lambda <- par_alphaweibull["lambda"]
            gamma <- par_alphaweibull["gamma"]
            delta_k <- par_alphaweibull[
              c(which(stringr::str_detect(names(par_alphaweibull),"delta_k")))]

            start <- c(beta0 = as.numeric(beta0),
                       beta_k = as.numeric(beta_k),
                       lambda = as.numeric(lambda),
                       gamma = as.numeric(gamma),
                       delta_k = as.numeric(delta_k))

            names(start) <- names(theta2)



          }else if (pophaz.alpha == FALSE & n_z_pcured <= 0) {
            beta0 <- par_alphaweibull["beta0"]
            lambda <- par_alphaweibull["lambda"]
            gamma <- par_alphaweibull["gamma"]

            start <- c(beta0 = as.numeric(beta0),
                       lambda = as.numeric(lambda),
                       gamma = as.numeric(gamma))

            names(start) <- names(theta2)
          }


          add.arguments <- function(f, params) {
            t <- paste("arg <- alist(",
                       paste(sapply(1:length(params),
                                    function(i)
                                      paste(paste('`',names(params),"`", sep = "")[i], "=", sep = "")),
                             collapse = ","),
                       ")", sep = "")
            formals(f) <- eval(parse(text = t))
            f
          }

          ll_fn_alphaweibull3noalpha2 <- add.arguments(ll_fn_alphaweibull3, start)
          bbmle::parnames(ll_fn_alphaweibull3noalpha2) <- names(theta2)

          if (pophaz.alpha) {
            pophaz.alpha0 <- theta[length(theta)]

          }else {
            pophaz.alpha0 <- 1
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
            minuslogl = ll_fn_alphaweibull3noalpha2,
            fixed = fixed,
            optimizer = optimizer,
            method = method_opt,
            control = list(maxit = maxit_opt,
                           trace = trace))

        }
      }
    }
  }


  if (inherits(optimized, "try-error")) {
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
    par_alphaweibull <- genpar_alphaweibull(par, z_pcured,
                                            pophaz.alpha,
                                            z_pophaz.alpha,
                                            z_pophaz.alpha_id)

    if (pophaz.alpha == TRUE & n_z_pcured > 0) {
      beta0 <- par_alphaweibull["beta0"]
      beta_k <- par_alphaweibull[
        c(which(stringr::str_detect(names(par_alphaweibull),"beta_k")))]
      lambda <- par_alphaweibull["lambda"]
      gamma <- par_alphaweibull["gamma"]
      delta_k <- par_alphaweibull[
        c(which(stringr::str_detect(names(par_alphaweibull),"delta_k")))]
      pophaz.alpha0 <- par_alphaweibull["pophaz.alpha0"]
      if (length(z_pophaz.alpha_id) > 0) {
        pophaz.alpha_k <- par_alphaweibull[
          c(which(stringr::str_detect(
            names(par_alphaweibull),"pophaz.alpha_k")))]

        start <- c(beta0 = as.numeric(beta0),
                   beta_k = as.numeric(beta_k),
                   lambda = as.numeric(lambda),
                   gamma = as.numeric(gamma),
                   delta_k = as.numeric(delta_k),
                   pophaz.alpha0 = as.numeric(pophaz.alpha0),
                   pophaz.alpha_k = as.numeric(pophaz.alpha_k))
      }else{
        start <- c(beta0 = as.numeric(beta0),
                   beta_k = as.numeric(beta_k),
                   lambda = as.numeric(lambda),
                   gamma = as.numeric(gamma),
                   delta_k = as.numeric(delta_k),
                   pophaz.alpha0 = as.numeric(pophaz.alpha0))
      }



      names(start) <- names(theta2)



    }else if (pophaz.alpha == TRUE & n_z_pcured <= 0) {
      beta0 <- par_alphaweibull["beta0"]
      lambda <- par_alphaweibull["lambda"]
      gamma <- par_alphaweibull["gamma"]
      pophaz.alpha0 <- par_alphaweibull["pophaz.alpha0"]
      if (length(z_pophaz.alpha_id) > 0) {
        pophaz.alpha_k <- par_alphaweibull[
          c(which(stringr::str_detect(
            names(par_alphaweibull),"pophaz.alpha_k")))]

        start <- c(beta0 = as.numeric(beta0),
                   lambda = as.numeric(lambda),
                   gamma = as.numeric(gamma),
                   pophaz.alpha0 = as.numeric(pophaz.alpha0),
                   pophaz.alpha_k = as.numeric(pophaz.alpha_k))
      }else{
        start <- c(beta0 = as.numeric(beta0),
                   lambda = as.numeric(lambda),
                   gamma = as.numeric(gamma),
                   pophaz.alpha0 = as.numeric(pophaz.alpha0))
      }



      names(start) <- names(theta2)



    }else if (pophaz.alpha == FALSE & n_z_pcured > 0) {
      beta0 <- par_alphaweibull["beta0"]
      beta_k <- par_alphaweibull[
        c(which(stringr::str_detect(names(par_alphaweibull),"beta_k")))]
      lambda <- par_alphaweibull["lambda"]
      gamma <- par_alphaweibull["gamma"]
      delta_k <- par_alphaweibull[
        c(which(stringr::str_detect(names(par_alphaweibull),"delta_k")))]

      start <- c(beta0 = as.numeric(beta0),
                 beta_k = as.numeric(beta_k),
                 lambda = as.numeric(lambda),
                 gamma = as.numeric(gamma),
                 delta_k = as.numeric(delta_k))

      names(start) <- names(theta2)



    }else if (pophaz.alpha == FALSE & n_z_pcured <= 0) {
      beta0 <- par_alphaweibull["beta0"]
      lambda <- par_alphaweibull["lambda"]
      gamma <- par_alphaweibull["gamma"]

      start <- c(beta0 = as.numeric(beta0),
                 lambda = as.numeric(lambda),
                 gamma = as.numeric(gamma))

      names(start) <- names(theta2)
    }

    varcov_star <- optimized@vcov
    if (any(inherits(optimized@vcov, "try-error"))) {
      varcov_star <- try(solve(optimized@details$hessian), TRUE)
    }
    std_err_star <- try(sqrt(diag(varcov_star)), TRUE)

  } else {
    par <- optimized$par
    AIC <- 2*optimized$value + 2*length(par)
    loglik <- -optimized$value
    evaluations <- optimized$counts
    iterations <- evaluations
    convergence <- optimized$convergence
    message <- optimized$message
    SD <- numDeriv::hessian(ll_fn_alphaweibull, par)

    if (hessian_varcov) {
      varcov_star <- try(solve(SD), TRUE)
      std_err_star <- try(sqrt(diag(varcov_star)), TRUE)

    } else{
      stop("optimization process requires the use of  'hessian_varcov =TRUE'")
    }

    if (any(inherits(varcov_star, "try-error"))) {
      varcov_star   <- try(solve(optimized$hessian), TRUE)
      std_err_star <- try(sqrt(diag(varcov_star)), TRUE)
    }

  }

  thetaest <- matrix(par, nrow = 1)
  n_z_pcured <- ncol(z_pcured)
  n_z_ucured <- ncol(z_ucured)

  if (n_z_pcured > 0 & n_z_ucured > 0) {
    if (pophaz.alpha) {
      if (length(z_pophaz.alpha_id) > 0) {
        colnames(thetaest) <- c(
          'beta0',
          paste('beta', 1:n_z_pcured, colnames(z_pcured) , sep = "_"),
          'lambda', 'gamma',
          paste('delta', 1:n_z_ucured, colnames(z_ucured), sep = "_"),
          'pophaz.alpha0', paste('pophaz.alpha', 1:length(z_pophaz.alpha_id), colnames(z_pophaz.alpha) , sep = "_"))
      }else{
        colnames(thetaest) <- c(
          'beta0',
          paste('beta', 1:n_z_pcured, colnames(z_pcured) , sep = "_"),
          'lambda', 'gamma',
          paste('delta', 1:n_z_ucured, colnames(z_ucured), sep = "_"),
          'pophaz.alpha0')
      }


    } else{
      colnames(thetaest) <- c(
        'beta0',
        paste('beta', 1:n_z_pcured, colnames(z_pcured) , sep = "_"),
        'lambda', 'gamma',
        paste('delta', 1:n_z_ucured, colnames(z_ucured), sep = "_"))
    }



  } else if (n_z_pcured == 0 & n_z_ucured > 0) {
    if (pophaz.alpha) {
      if (length(z_pophaz.alpha_id) > 0) {
        colnames(thetaest) <- c(
          'beta0',
          'lambda', 'gamma',
          paste('delta', 1:n_z_ucured, colnames(z_ucured), sep = "_"),
          'pophaz.alpha0', paste('pophaz.alpha', 1:length(z_pophaz.alpha_id), colnames(z_pophaz.alpha) , sep = "_"))
      }else{
        colnames(thetaest) <- c(
          'beta0',
          'lambda', 'gamma',
          paste('delta', 1:n_z_ucured, colnames(z_ucured), sep = "_"),
          'pophaz.alpha0')
      }


    } else{
      colnames(thetaest) <- c(
        'beta0',
        'lambda', 'gamma',
        paste('delta', 1:n_z_ucured, colnames(z_ucured), sep = "_"))
    }



  } else if (n_z_pcured > 0 & n_z_ucured == 0) {
    if (pophaz.alpha) {
      if (length(z_pophaz.alpha_id) > 0) {
        colnames(thetaest) <- c(
          'beta0',
          paste('beta', 1:n_z_pcured, colnames(z_pcured) , sep = "_"),
          'lambda', 'gamma',
          'pophaz.alpha0', paste('pophaz.alpha', 1:length(z_pophaz.alpha_id), colnames(z_pophaz.alpha) , sep = "_"))
      }else {
        colnames(thetaest) <- c(
          'beta0',
          paste('beta', 1:n_z_pcured, colnames(z_pcured) , sep = "_"),
          'lambda', 'gamma',
          'pophaz.alpha0')
      }


    } else{
      colnames(thetaest) <- c(
        'beta0',
        paste('beta', 1:n_z_pcured, colnames(z_pcured) , sep = "_"),
        'lambda', 'gamma')
    }



  } else if (n_z_pcured == 0 & n_z_ucured == 0) {

    if (pophaz.alpha) {

      if (length(z_pophaz.alpha_id) > 0) {
        colnames(thetaest) <- c('beta0',
                                'lambda',
                                'gamma',
                                'pophaz.alpha0',
                                paste('pophaz.alpha',
                                      1:length(z_pophaz.alpha_id),
                                      colnames(z_pophaz.alpha) ,
                                      sep = "_"))
      } else{
        colnames(thetaest) <- c('beta0', 'lambda', 'gamma', 'pophaz.alpha0')
      }
    } else{
      colnames(thetaest) <- c('beta0', 'lambda', 'gamma')
    }
  }


  if (pophaz.alpha) {
    expvar <- c("lambda", "gamma", "pophaz.alpha0")
  } else {
    expvar <- c("lambda", "gamma")
  }

  idexpvar <- sapply(1:length(expvar), function(i)
    which(stringr::str_detect(string = colnames(thetaest),
                              pattern  = expvar[i])))

  estimates <- (thetaest)

  estimates[1] <- 1/c(1 + exp(-(thetaest[1])))
  if (n_z_pcured > 0) {
    estimates[2:(n_z_pcured + 1)] <-
      1 / c(1 + exp(-(thetaest[1] + thetaest[2:(n_z_pcured + 1)])))

  }

  estimates[idexpvar] <- exp(thetaest[idexpvar])

  dfestimates <- estimates
  dfestimates[-c(idexpvar)] <- 1
  dfestimates[1] <-  exp(-(thetaest[1]))/c(1 + exp(-(thetaest[1])))^2
  if (n_z_pcured > 0) {
    dfestimates[2:(n_z_pcured + 1)] <-
      exp(-(thetaest[1] + thetaest[2:(n_z_pcured + 1)])) / c(1 + exp(-(thetaest[1] +
                                                                         thetaest[2:(n_z_pcured + 1)]))) ^ 2
  }
  varcov <- varcov_star * c(dfestimates)^2
  std_err <- sqrt(diag(varcov))

  if (method_opt == "all.methods") {
    return(
      list(
        coefficients = thetaest,
        estimates = estimates,
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
      list(
        coefficients = thetaest,
        estimates = estimates,
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
  invisible()
}
