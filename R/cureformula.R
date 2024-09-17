#' @import stats stringr
cureformula <- function(f = f,
                        temp = temp,
                        Terms = Terms,
                        model = model,
                        link_tau = link_tau,
                        dist = dist,
                        pophaz = pophaz,
                        cumpophaz = cumpophaz,
                        method_opt = method_opt,
                        pophaz.alpha = pophaz.alpha,
                        maxit_opt = maxit_opt,
                        optim_func = optim_func,
                        gradient = gradient,
                        hessian_varcov = hessian_varcov,
                        ncoor_des = ncoor_des,
                        trace = trace,
                        optim_fixed = optim_fixed,
                        optimizer = optimizer,
                        init = init,
                        nvalues = nvalues,
                        iter_eps = iter_eps,
                        clustertype = clustertype,
                        nproc = nproc,
                        subset = subset,
                        na.action = na.action,
                        sign_delta = sign_delta) {
  m_eval <- eval(temp, parent.frame())
  data <- model.extract(m_eval, "response")
  event <- model.extract(m_eval, "response")[, "status"]
  time <- model.extract(m_eval, "response")[, "time"]
  if(any(time==0)){
    time[time==0]<-1/730
    message("Some individuals have their time equal to 0. \n For them we set time=1/730 (half a day)")
  }
  xmax <- max(time)
  myvarnames <- colnames(model.matrix(Terms, m_eval)[,-1, drop = FALSE])
  z_alpha_id <- which(str_detect(c(myvarnames),
                                          pattern = "z_alpha"))
  z_tau_id <- which(str_detect(c(myvarnames),
                                        pattern = "z_tau"))

  z_c_id <- which(str_detect(c(myvarnames),
                                      pattern = "z_c"))



  if (length(z_alpha_id) > 0) {
    z_alpha <- as.data.frame(
      model.matrix(Terms, m_eval)[,-1, drop = FALSE][, c(z_alpha_id)])
    colnames(z_alpha) <- c(str_remove(myvarnames[c(z_alpha_id)], "z_alpha"))
    z_alpha <- as.matrix(z_alpha)
  }else{
    z_alpha <- matrix(
      nrow = nrow(model.matrix(Terms, m_eval)[,-1, drop = FALSE]),
      ncol = 0)
  }

  if (length(z_tau_id) > 0) {
    z_tau <- as.data.frame(
      model.matrix(Terms, m_eval)[,-1, drop = FALSE][, c(z_tau_id)])
    colnames(z_tau) <- c(str_remove(myvarnames[c(z_tau_id)],"z_tau"))
    z_tau <- as.matrix(z_tau)
  }
  else{
    z_tau <- matrix(
      nrow = nrow(model.matrix(Terms, m_eval)[,-1, drop = FALSE]),
      ncol = 0)
  }
  if (length(z_c_id) > 0) {
    z_c <- as.data.frame(
      model.matrix(Terms, m_eval)[,-1, drop = FALSE][, c(z_c_id)])
    colnames(z_c) <- c(str_remove(myvarnames[c(z_c_id)],"z_c"))
    z_c <- as.matrix(z_c)

  }else{
    z_c <- matrix(nrow = nrow(model.matrix(Terms, m_eval)[,-1, drop = FALSE]),
                  ncol = 0)
  }


  if (dist == "tneh"  &
      link_tau == "loglinear" &
      model == "nmixture" ) {
    n_z_tau <- ncol(z_tau)
    n_z_alpha <- ncol(z_alpha)
    n_z_tau_ad <- n_z_tau - 1
    n_z_alpha_ad <- n_z_alpha - 1

    if (pophaz.alpha) {
      n_par <- 1 + n_z_alpha + 2 + n_z_tau + 1
    }else{
      n_par <- 1 + n_z_alpha + 2 + n_z_tau
    }

    valinit <- is.null(init)

    if (valinit) {

      if (pophaz.alpha) {

        theta_lower <- c(-100,
                         rep(-100, n_z_alpha),
                         -100, -100,
                         rep(-100, n_z_tau),
                         -100)
        theta_upper <- c(
          rep(100, (1 + n_z_alpha + 1)),
          100,
          rep(100, (n_z_tau + 1)))


        theta_lower2 <- c(-10,
                          rep(-10, n_z_alpha),
                          -10, -10,
                          rep(-10, n_z_tau),
                          -10)
        theta_upper2 <- c(
          rep(10, (1 + n_z_alpha + 1)),
          10,
          rep(10, (n_z_tau + 1)))


        if (!is.null(theta_lower) && !is.null(theta_upper)) {
          theta_init <- gen_start_values(nvalues = 1,
                                         min = rep(-10, length(theta_lower)),
                                         max = rep(10, length(theta_upper)),
                                         npar = length(theta_lower))
        } else if (ncoor_des == "iter") {
          stop("add bounds for initial values")
        }

      } else{
        theta_lower <- c(-100,
                         rep(-100, n_z_alpha),
                         -100, -100,
                         rep(-100, n_z_tau))
        theta_upper <- c(
          rep(100, (1 + n_z_alpha + 1)),
          100,
          rep(100, (n_z_tau )))



        theta_lower2 <- c(-10,
                          rep(-10, n_z_alpha),
                          -10, -10,
                          rep(-10, n_z_tau))
        theta_upper2 <- c(
          rep(10, (1 + n_z_alpha + 1)),
          10,
          rep(10, (n_z_tau )))


        if (!is.null(theta_lower) && !is.null(theta_upper)) {
          theta_init <- gen_start_values(nvalues = 1,
                                         min = rep(-100, length(theta_lower)),
                                         max = rep(100, length(theta_upper)),
                                         npar = length(theta_lower))
        } else if (ncoor_des == "iter") {
          stop("add bounds for initial values")
        }

      }


    } else {

      if (!is.null(init$theta_init)) {
        if (n_par != length(init$theta_init)) {
          stop("check the number of initial values provided!")
        }else {
          theta_init = init$theta_init
        }
      }

      if (!is.null(init$theta_lower)) {
        if (n_par != length(init$theta_lower)) {
          stop("check the number of initial values provided (lower bound)!")
        }else {
          theta_lower = init$theta_lower
        }
      }

      if (!is.null(init$theta_upper)) {
        if (n_par != length(init$theta_upper)) {
          stop("check the number of initial values provided (lower bound)!")
        }else {
          theta_upper = init$theta_upper
        }
      }



    }

    if (valinit) {
      mod <- list()
      par <- list()
      cptM <- 0
      nConverg <- 0
      num_NConv <- vector("numeric")
      i <- 1


      start_val <- gen_start_values(nvalues = nvalues,
                                    min = theta_lower2,
                                    max = theta_upper2,
                                    npar = length(theta_lower2))

      while (i <= nvalues && (cptM < 1)) {
        #cat("init", start_val[i, ],"lower", theta_lower, "upper",theta_upper,"\n")

        mod[[i]] <- try(fit.opt.maxim.multneh(x = time,
                                              d = event,
                                              z_tau = z_tau,
                                              z_alpha = z_alpha,
                                              theta_init =  start_val[i, ],
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
                                              trace = trace), TRUE)

        if (
          !anyNA(mod[[i]]$std_err_star) &&
          !anyNA(mod[[i]]$coefficients) &&
          !any(mod[[i]]$std_err_star > 1000)) {
          par <- mod[[i]]
          cptM <- cptM + 1
          num_NConv[i] <- NA
        }else{
          par <- mod[[i]]
          #cat(paste("non convergence with inititial values", i,"\n", sep = " "))
          nConverg <- nConverg + 1
          num_NConv[i] <- i
        }
        #cat(paste("next evaluation with initial values =  ", c(i + 1), "\n", sep = " "))
        i <- i + 1
      }

    }else {
      par <- fit.opt.maxim.multneh(x = time,
                                   d = event,
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
                                   trace = trace)
    }




  }
  else if (dist == "tneh" & link_tau == "loglinear2" & model == "nmixture" ) {

    n_z_tau <- ncol(z_tau)
    n_z_alpha <- ncol(z_alpha)
    n_z_tau_ad <- n_z_tau - 1
    n_z_alpha_ad <- n_z_alpha - 1
    if (pophaz.alpha) {
      n_par <- 1 + n_z_alpha + 2 + n_z_tau + 1
    }else{
      n_par <- 1 + n_z_alpha + 2 + n_z_tau
    }

    valinit <- is.null(init)

    if (valinit) {
      if (pophaz.alpha) {

        theta_lower <- c(0,
                         rep(-5, n_z_alpha),
                         1, 0,
                         rep(-5, n_z_tau),
                         -100)
        theta_upper <- c(
          rep(100, (1 + n_z_alpha)),
          100,
          c(100),
          rep(100, (n_z_tau + 1)))



        theta_lower2 <- c(0,
                          rep(-5, n_z_alpha),
                          1, 0,
                          rep(-5, n_z_tau),
                          -10)

        theta_upper2 <- c(
          rep(10, (1 + n_z_alpha)),
          10,
          c(10),
          rep(10, (n_z_tau + 1)))


        if (!is.null(theta_lower) && !is.null(theta_upper)) {
          theta_init <- gen_start_values(nvalues = 1,
                                         min = theta_lower,
                                         max = theta_upper,
                                         npar = length(theta_lower))
        } else if (ncoor_des == "iter") {
          stop("add bounds for initial values")
        }

      } else{
        theta_lower <- c(0,
                         rep(-5, n_z_alpha),
                         1, 0,
                         rep(-5, n_z_tau))
        theta_upper <- c(
          rep(100, (1 + n_z_alpha)),
          100,
          (100),
          rep(100, (n_z_tau)))

        theta_lower2 <- c(0,
                          rep(-5, n_z_alpha),
                          1, 0,
                          rep(-18, n_z_tau))


        theta_upper2 <- c(
          rep(10, (1 + n_z_alpha)),
          10,
          (10),
          rep(10, (n_z_tau)))

        if (!is.null(theta_lower) && !is.null(theta_upper)) {
          theta_init <- gen_start_values(nvalues = 1,
                                         min = theta_lower,
                                         max = theta_upper,
                                         npar = length(theta_lower))
        } else if (ncoor_des == "iter") {
          stop("add bounds for initial values")
        }

      }

    }else{
      if (!is.null(init$theta_init)) {
        if (n_par != length(init$theta_init)) {
          stop("check the number of initial values provided!")
        }else {
          theta_init = init$theta_init
        }
      }

      if (!is.null(init$theta_lower)) {
        if (n_par != length(init$theta_lower)) {
          stop("check the number of initial values provided (lower bound)!")
        }else {
          theta_lower = init$theta_lower
        }
      }

      if (!is.null(init$theta_upper)) {
        if (n_par != length(init$theta_upper)) {
          stop("check the number of initial values provided (lower bound)!")
        }else {
          theta_upper = init$theta_upper
        }
      }

    }

    if (valinit) {
      mod <- list()
      par <- list()
      cptM <- 0
      nConverg <- 0
      num_NConv <- vector("numeric")
      i <- 1


      start_val <- gen_start_values(nvalues = nvalues,
                                    min = theta_lower2,
                                    max = theta_upper2,
                                    npar = length(theta_lower))

      while (i <= nvalues && (cptM < 1)) {
        #cat("init", start_val[i, ],"lower", theta_lower, "upper",theta_upper,"\n")

        mod[[i]] <- try(fit.opt.maxim.adtneh2bis(x = time,
                                                 d = event,
                                                 z_tau = z_tau,
                                                 z_alpha = z_alpha,
                                                 theta_init = start_val[i, ],
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
                                                 optimizer = optimizer), TRUE)

        if (!anyNA(mod[[i]]$std_err_star) &&
          !anyNA(mod[[i]]$coefficients) &&
          !any(mod[[i]]$std_err_star > 1000)) {
          par <- mod[[i]]
          cptM <- cptM + 1
          num_NConv[i] <- NA
        }else{
          par <- mod[[i]]
          #cat(paste("non convergence with inititial values", i,"\n", sep = " "))
          nConverg <- nConverg + 1
          num_NConv[i] <- i
        }
        #cat(paste("next evaluation with initial values =  ",c(i + 1), "\n", sep = " "))
        i <- i + 1
      }

    }else{
      par <- fit.opt.maxim.adtneh2bis(x = time,
                                      d = event,
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
                                      optimizer = optimizer)

    }


  }else if (dist == "tneh" & link_tau == "linear"  & model == "nmixture") {

    n_z_tau <- ncol(z_tau)
    n_z_alpha <- ncol(z_alpha)
    n_z_tau_ad <- n_z_tau - 1
    n_z_alpha_ad <- n_z_alpha - 1
    if (pophaz.alpha) {
      n_par <- 1 + n_z_alpha + 2 + n_z_tau + 1
    }else{
      n_par <- 1 + n_z_alpha + 2 + n_z_tau
    }

    valinit <- is.null(init)

    if (valinit) {
      if (pophaz.alpha) {

        theta_lower <- c(0,
                         rep(-5, n_z_alpha),
                         1, 0,
                         rep(-5, n_z_tau),
                         -100)
        theta_upper <- c(
          rep(100, (1 + n_z_alpha)),
          100,
          c(100),
          rep(100, (n_z_tau + 1)))



        theta_lower2 <- c(0,
                          rep(-5, n_z_alpha),
                          1, 0,
                          rep(-5, n_z_tau),
                          -10)

        theta_upper2 <- c(
          rep(10, (1 + n_z_alpha)),
          10,
          c(10),
          rep(10, (n_z_tau + 1)))


        if (!is.null(theta_lower) && !is.null(theta_upper)) {
          theta_init <- gen_start_values(nvalues = 1,
                                         min = theta_lower,
                                         max = theta_upper,
                                         npar = length(theta_lower))
        } else if (ncoor_des == "iter") {
          stop("add bounds for initial values")
        }

      } else{
        theta_lower <- c(0,
                         rep(-5, n_z_alpha),
                         1, 0,
                         rep(-5, n_z_tau))
        theta_upper <- c(
          rep(100, (1 + n_z_alpha)),
          100,
          (100),
          rep(100, (n_z_tau)))

        theta_lower2 <- c(0,
                          rep(-5, n_z_alpha),
                          1, 0,
                          rep(-18, n_z_tau))


        theta_upper2 <- c(
          rep(10, (1 + n_z_alpha)),
          10,
          (10),
          rep(10, (n_z_tau)))

        if (!is.null(theta_lower) && !is.null(theta_upper)) {
          theta_init <- gen_start_values(nvalues = 1,
                                         min = theta_lower,
                                         max = theta_upper,
                                         npar = length(theta_lower))
        } else if (ncoor_des == "iter") {
          stop("add bounds for initial values")
        }

      }

    }else{
      if (!is.null(init$theta_init)) {
        if (n_par != length(init$theta_init)) {
          stop("check the number of initial values provided!")
        }else {
          theta_init = init$theta_init
        }
      }

      if (!is.null(init$theta_lower)) {
        if (n_par != length(init$theta_lower)) {
          stop("check the number of initial values provided (lower bound)!")
        }else {
          theta_lower = init$theta_lower
        }
      }

      if (!is.null(init$theta_upper)) {
        if (n_par != length(init$theta_upper)) {
          stop("check the number of initial values provided (lower bound)!")
        }else {
          theta_upper = init$theta_upper
        }
      }

    }

    if (valinit) {
      mod <- list()
      par <- list()
      cptM <- 0
      nConverg <- 0
      num_NConv <- vector("numeric")
      i <- 1


      start_val <- gen_start_values(nvalues = nvalues,
                                    min = theta_lower2,
                                    max = theta_upper2,
                                    npar = length(theta_lower))

      while (i <= nvalues && (cptM < 1)) {
        #cat("init", start_val[i, ],"lower", theta_lower, "upper",theta_upper,"\n")

        mod[[i]] <- try(fit.opt.maxim.adtneh2(x = time,
                                              d = event,
                                              z_tau = z_tau,
                                              z_alpha = z_alpha,
                                              theta_init = start_val[i, ],
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
                                              optimizer = optimizer), TRUE)

        if (
          !anyNA(mod[[i]]$std_err_star) &&
          !anyNA(mod[[i]]$coefficients) &&
          !any(mod[[i]]$std_err_star > 1000)) {
          par <- mod[[i]]
          cptM <- cptM + 1
          num_NConv[i] <- NA
        }else{
          par <- mod[[i]]
          #cat(paste("non convergence with inititial values", i,"\n", sep = " "))
          nConverg <- nConverg + 1
          num_NConv[i] <- i
        }
        #cat(paste("next evaluation with initial values =  ",c(i + 1), "\n", sep = " "))
        i <- i + 1
      }

    }
    else{
      par <- fit.opt.maxim.adtneh2(x = time,
                                   d = event,
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
                                   optimizer = optimizer)

    }


  }   else if (model == "mixture" & dist == "weib")
    {

    z_pophaz.alpha_id <- which(str_detect(c(myvarnames),
                                                   pattern = "z_pophaz.alpha"))
    if (length(z_pophaz.alpha_id) > 0) {
      z_pophaz.alpha <- as.data.frame(
        model.matrix(Terms, m_eval)[, -1, drop = FALSE][, c(z_pophaz.alpha_id)])
      colnames(z_pophaz.alpha) <- c(str_remove(
        myvarnames[c(z_pophaz.alpha_id)], "z_pophaz.alpha"))
      z_pophaz.alpha <- as.matrix(z_pophaz.alpha)
      z_ucured <-  model.matrix(f, m_eval, rhs = 1)[,-1,drop = FALSE]

      z_pcured <- model.matrix(f, m_eval, rhs = 2)[,-1,drop = FALSE]

    }else{
      z_pophaz.alpha <- matrix(nrow = nrow(model.matrix(Terms, m_eval)[,-1, drop = FALSE]), ncol = 0)
      z_ucured <-  model.matrix(f, m_eval, rhs = 1)[,-1,drop = FALSE]

      z_pcured <- model.matrix(f, m_eval, rhs = 2)[,-1,drop = FALSE]

    }

    if (pophaz.alpha) {
      n_par <- (ncol(z_pcured) + ncol(z_ucured) + 3) + 1 + length(z_pophaz.alpha_id)
    }else {
      n_par <- (ncol(z_pcured) + ncol(z_ucured) + 3)
    }

    if (is.null(init)) {
      theta_init = rep(0, n_par)
      theta_lower = rep(-Inf, n_par)
      theta_upper = rep(Inf, n_par)
    }else{
      if (!is.null(init$theta_init)) {
        if (n_par != length(init$theta_init)) {
          stop("check the number of initial values provided!")
        }else {
          theta_init = init$theta_init
        }
      }

      if (!is.null(init$theta_lower)) {
        if (n_par != length(init$theta_lower)) {
          stop("check the number of initial values provided (lower bound)!")
        }else {
          theta_lower = init$theta_lower
        }
      }

      if (!is.null(init$theta_upper)) {
        if (n_par != length(init$theta_upper)) {
          stop("check the number of initial values provided (lower bound)!")
        }else {
          theta_upper = init$theta_upper
        }
      }

    }
    par <- fit.opt.maxim_mixture_alpha(x = time,
                                       d = event,
                                       z_pcured = z_pcured,
                                       z_ucured = z_ucured,
                                       z_pophaz.alpha = z_pophaz.alpha,
                                       z_pophaz.alpha_id = z_pophaz.alpha_id,
                                       theta_init = theta_init,
                                       theta_lower = theta_lower,
                                       theta_upper = theta_upper,
                                       pophaz = pophaz,
                                       cumpophaz = cumpophaz,
                                       method_opt = method_opt,
                                       pophaz.alpha = pophaz.alpha,
                                       maxit_opt = maxit_opt,
                                       optim_func = optim_func,
                                       gradient = gradient,
                                       hessian_varcov = hessian_varcov,
                                       ncoor_des = ncoor_des,
                                       trace = trace,
                                       optim_fixed = optim_fixed,
                                       optimizer = optimizer)

  } else if (model == "mixture" & dist == "eweib")
    {
    z_pophaz.alpha_id <- which(str_detect(c(myvarnames),
                                                   pattern = "z_pophaz.alpha"))
    if (length(z_pophaz.alpha_id) > 0) {
      z_pophaz.alpha <- as.data.frame(
        model.matrix(Terms, m_eval)[, -1, drop = FALSE][, c(z_pophaz.alpha_id)])
      colnames(z_pophaz.alpha) <- c(str_remove(
        myvarnames[c(z_pophaz.alpha_id)], "z_pophaz.alpha"))
      z_pophaz.alpha <- as.matrix(z_pophaz.alpha)
      z_ucured <-  model.matrix(f, m_eval, rhs = 1)[,-1,drop = FALSE]

      z_pcured <- model.matrix(f, m_eval, rhs = 2)[,-1,drop = FALSE]

    }else{
      z_pophaz.alpha <- matrix(nrow = nrow(model.matrix(Terms, m_eval)[,-1, drop = FALSE]), ncol = 0)

      z_ucured <-  model.matrix(f, m_eval, rhs = 1)[,-1,drop = FALSE]

      z_pcured <- model.matrix(f, m_eval, rhs = 2)[,-1,drop = FALSE]
    }

    if (pophaz.alpha) {
      n_par <- (ncol(z_pcured) + ncol(z_ucured) + 3) + 2 + length(z_pophaz.alpha_id)
    }else {
      n_par <- (ncol(z_pcured) + ncol(z_ucured) + 4)
    }

    if (is.null(init)) {
      theta_init = rep(0, n_par)
      theta_lower = rep(-Inf, n_par)
      theta_upper = rep(Inf, n_par)
    }else{
      if (!is.null(init$theta_init)) {
        if (n_par != length(init$theta_init)) {
          stop("check the number of initial values provided!")
        }else {
          theta_init = init$theta_init
        }
      }

      if (!is.null(init$theta_lower)) {
        if (n_par != length(init$theta_lower)) {
          stop("check the number of initial values provided (lower bound)!")
        }else {
          theta_lower = init$theta_lower
        }
      }

      if (!is.null(init$theta_upper)) {
        if (n_par != length(init$theta_upper)) {
          stop("check the number of initial values provided (lower bound)!")
        }else {
          theta_upper = init$theta_upper
        }
      }

    }

    par1 <- fit.opt.maxim_mixture_alpha(x = time,
                                        d = event,
                                        z_pcured = z_pcured,
                                        z_ucured = z_ucured,
                                        z_pophaz.alpha = z_pophaz.alpha,
                                        z_pophaz.alpha_id = z_pophaz.alpha_id,
                                        theta_init = theta_init[-4],
                                        theta_lower = theta_lower[-4],
                                        theta_upper = theta_upper[-4],
                                        pophaz = pophaz,
                                        cumpophaz = cumpophaz,
                                        method_opt = method_opt,
                                        pophaz.alpha = pophaz.alpha,
                                        maxit_opt = maxit_opt,
                                        optim_func = optim_func,
                                        gradient = gradient,
                                        hessian_varcov = hessian_varcov,
                                        ncoor_des = ncoor_des,
                                        trace = trace,
                                        optim_fixed = optim_fixed,
                                        optimizer = optimizer,
                                        sign_delta = sign_delta)
    if (length(par1$coef) > 3) {
      parint <- c(par1$coef[1:3], 0.1, par1$coef[4:length(par1$coef)])

    }else {
      parint <- c(par1$coef[1:3], 0.1)

    }

    par <- fit_optim.maxim_mixture_eweib_alpha(x = time,
                                               d = event,
                                               z_pcured = z_pcured,
                                               z_ucured = z_ucured,
                                               z_pophaz.alpha = z_pophaz.alpha,
                                               z_pophaz.alpha_id = z_pophaz.alpha_id,
                                               theta_init = parint,
                                               theta_lower = theta_lower,
                                               theta_upper = theta_upper,
                                               pophaz = pophaz,
                                               cumpophaz = cumpophaz,
                                               method_opt = method_opt,
                                               pophaz.alpha = pophaz.alpha,
                                               maxit_opt = maxit_opt,
                                               optim_func = optim_func,
                                               gradient = gradient,
                                               hessian_varcov = hessian_varcov,
                                               ncoor_des = ncoor_des,
                                               trace = trace,
                                               optim_fixed = optim_fixed,
                                               optimizer = optimizer,
                                               sign_delta = sign_delta)


  }


  par$n.event <- sum(event)
  par$n.obs <- nrow(data)
  par$model <- model
  par$Terms <- Terms
  #par$f <- f
  par$pophaz.alpha <- pophaz.alpha
  par$pophaz <- pophaz
  par$cumpophaz <- cumpophaz
  par$dist <- dist
  par$xmax <- xmax

  if (model == "mixture") {
    par$z_pcured <- z_pcured
    par$z_ucured <- z_ucured
  }else{
    if (dist == "tneh") {
      par$z_tau <- z_tau
      par$link_tau <- link_tau
      par$z_alpha <- z_alpha
      par$z_c <- z_c
    }

  }
par
}
