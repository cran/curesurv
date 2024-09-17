predictform <- function(object,
                        newdata = NULL,
                        xmax = NULL,
                        level = 0.975,
                        epsilon = 0.05, sign_delta = 1,
                        ...) {
  if (!inherits(object, "curesurv"))
    stop("Primary argument much be a curesurv object")
  Call <- match.call()
  newcall <- Call[c(1, match(c("newdata"),
                             names(Call), nomatch = 0))]
  names(newcall)[2] <- "data"
  newcall[[1]] <- as.name("model.frame")

  Terms <- newcall$formula <- object$Terms

  newcall$formula <- Terms
  m_eval <- eval(newcall, parent.frame())
  event <- stats::model.extract(m_eval, "response")[, "status"]
  time <- stats::model.extract(m_eval, "response")[, "time"]
  myvarnames <-
    colnames(stats::model.matrix(newcall$formula, m_eval)[, -1, drop = FALSE])


  f<-Formula::Formula(object$formula)





  if (object$model == "nmixture") {

    if ((object$dist == "tneh" &
         object$link_tau == "linear") | (object$dist == "tneh" &
                                         object$link_tau == "loglinear")) {
      z_alpha_id <- which(stringr::str_detect(c(myvarnames),
                                              pattern = "z_alpha"))

      z_tau_id <- which(stringr::str_detect(c(myvarnames),
                                            pattern = "z_tau"))




      if (length(z_alpha_id) > 0) {
        z_alpha <- as.data.frame(stats::model.matrix(Terms, m_eval)[, -1, drop = FALSE][, c(z_alpha_id)])
        if(ncol(z_alpha)!=length(z_alpha_id)){
          rownames(z_alpha) <- c(stringr::str_remove(myvarnames[c(z_alpha_id)], "z_alpha"))
        }else{
          colnames(z_alpha) <- c(stringr::str_remove(myvarnames[c(z_alpha_id)], "z_alpha"))
        }
        z_alpha <- as.matrix(z_alpha)
      } else{
        z_alpha <- matrix(nrow = nrow(stats::model.matrix(Terms,
                                                          m_eval)[, -1, drop = FALSE]),
                          ncol = 0)
      }

      if (length(z_tau_id) > 0) {
        z_tau <- as.data.frame(stats::model.matrix(Terms,
                                                   m_eval)[, -1, drop = FALSE][, c(z_tau_id)])
        if(ncol(z_tau)!=length(z_tau_id)){
          rownames(z_tau) <-c(stringr::str_remove(myvarnames[c(z_tau_id)], "z_tau"))
        }else{
          colnames(z_tau) <-c(stringr::str_remove(myvarnames[c(z_tau_id)], "z_tau"))
        }
        z_tau <- as.matrix(z_tau)

      } else{
        z_tau <- matrix(nrow = nrow(stats::model.matrix(Terms,
                                                        m_eval)[, -1, drop = FALSE]),
                        ncol = 0)
      }
    }

  } else if (object$model == "mixture") {
    if ( object$dist == "weib" |  object$dist == "eweib") {
      z_pcured <- stats::model.matrix(Terms,
                                      m_eval)[, -1, drop = FALSE]



      z_pophaz.alpha_id <- which(stringr::str_detect(c(myvarnames),
                                                     pattern = "z_pophaz.alpha"))
      if (length(z_pophaz.alpha_id) > 0) {
        z_pophaz.alpha <- as.data.frame(
          stats::model.matrix(Terms, m_eval)[, -1, drop = FALSE][, c(z_pophaz.alpha_id)])
        colnames(z_pophaz.alpha) <- c(stringr::str_remove(
          myvarnames[c(z_pophaz.alpha_id)], "z_pophaz.alpha"))
        z_pophaz.alpha <- as.matrix(z_pophaz.alpha)

        z_ucured <-  model.matrix(f, m_eval, rhs = 1)[,-1,drop = FALSE]
        z_pcured <- model.matrix(f, m_eval, rhs = 2)[,-1,drop = FALSE]

      }else{
        z_pophaz.alpha <- matrix(nrow = nrow(stats::model.matrix(Terms, m_eval)[,-1, drop = FALSE]), ncol = 0)

        z_ucured <-  model.matrix(f, m_eval, rhs = 1)[,-1,drop = FALSE]
        z_pcured <- model.matrix(f, m_eval, rhs = 2)[,-1,drop = FALSE]
      }

      if (object$pophaz.alpha) {
        n_par <- (ncol(z_pcured) + ncol(z_ucured) + 3) + 1 + length(z_pophaz.alpha_id)
      }else {
        n_par <- (ncol(z_pcured) + ncol(z_ucured) + 3)
      }



    }
  }

  theta <- object$coefficients



  if (object$model == "nmixture") {
    if ((object$dist == "tneh" &
         object$link_tau == "loglinear")) {

      cumLexctopred<-cumLexc_mul_topred(z_tau,
                         z_alpha,
                         time,
                         theta)
      netsurv <- exp(-cumLexctopred$cumhaz)
      netsurv_tau <-  exp(-cumLexctopred$cumhaz2)
      ex_haz <- lexc_mul(z_tau, z_alpha, time, theta)
      pt_cure <- netsurv_tau / netsurv
      tau <- cumLexctopred$tau

      if(length(z_alpha_id)>0&length(z_tau_id)==0){
        pred <- data.frame(time = time,
                           Z=z_alpha,
                           ex_haz = ex_haz,
                           netsurv = netsurv,
                           netsurv_tau = netsurv_tau,
                           pt_cure=pt_cure,
                           tau = tau)
      }else if(length(z_alpha_id)==0&length(z_tau_id)>0){
        pred <- data.frame(time = time,
                           Z=z_tau,
                           ex_haz = ex_haz,
                           netsurv = netsurv,
                           netsurv_tau = netsurv_tau,
                           pt_cure=pt_cure,
                           tau = tau)
      }
      else if(length(z_alpha_id)>0&length(z_tau_id)>0){
        rep_col<-intersect(colnames(z_tau),colnames(z_alpha))
        unique_col<-setdiff(colnames(z_tau),rep_col)
        Z<-cbind(z_alpha,z_tau[,unique_col,drop=FALSE])
        pred <- data.frame(time = time,
                           Z=Z,
                           ex_haz = ex_haz,
                           netsurv = netsurv,
                           netsurv_tau = netsurv_tau,
                           pt_cure=pt_cure,
                           tau = tau)
      }else{
        pred <- data.frame(time = time,
                           ex_haz = ex_haz,
                           netsurv = netsurv,
                           netsurv_tau = netsurv_tau,
                           pt_cure=pt_cure,
                           tau = tau)
      }
      pred <- predcall_tneh(object,
                            pred,
                            z_tau = z_tau,
                            z_alpha = z_alpha,
                            x = time,
                            xmax = xmax,
                            level = level,
                            epsilon = epsilon,cumLexc_topred=cumLexctopred)

    }else if ((object$dist == "tneh" &
               object$link_tau == "linear")) {

      cumLexc_topred<-cumLexc_ad2_topred(z_tau,
                                         z_alpha,
                                         time,
                                         theta)

      netsurv <- exp(-cumLexc_topred$cumhaz)
      netsurv_tau <-  exp(-cumLexc_topred$cumhaz2)
      ex_haz <- lexc_ad2(z_tau, z_alpha, time, theta)
      pt_cure <- netsurv_tau / netsurv
      tau <- cumLexc_topred$tau

      if(length(z_alpha_id)>0&length(z_tau_id)==0){
        pre_pred <- data.frame(time = time,
                           Z=z_alpha,
                           ex_haz = ex_haz,
                           netsurv = netsurv,
                           netsurv_tau = netsurv_tau,
                           pt_cure=pt_cure,
                           tau = tau)
      }else if(length(z_alpha_id)==0&length(z_tau_id)>0){
        pre_pred <- data.frame(time = time,
                           Z=z_tau,
                           ex_haz = ex_haz,
                           netsurv = netsurv,
                           netsurv_tau = netsurv_tau,
                           pt_cure=pt_cure,
                           tau = tau)
      }
      else if(length(z_alpha_id)>0&length(z_tau_id)>0){
        rep_col<-intersect(colnames(z_tau),colnames(z_alpha))
        unique_col<-setdiff(colnames(z_tau),rep_col)
        Z<-cbind(z_alpha,z_tau[,unique_col,drop=FALSE])
        pre_pred <- data.frame(time = time,
                           Z=Z,
                           ex_haz = ex_haz,
                           netsurv = netsurv,
                           netsurv_tau = netsurv_tau,
                           pt_cure=pt_cure,
                           tau = tau)
      }else{
        pre_pred <- data.frame(time = time,
                           ex_haz = ex_haz,
                           netsurv = netsurv,
                           netsurv_tau = netsurv_tau,
                           pt_cure=pt_cure,
                           tau = tau)
      }


      if (is.null(xmax)) {
        xmax <- object$xmax
      }
      pred <- predcall_tneh(object,
                            pre_pred,
                            z_tau = z_tau,
                            z_alpha = z_alpha,
                            x = time,
                            xmax = xmax,
                            level = level,
                            epsilon = epsilon,
                            cumLexc_topred=cumLexc_topred)

    } else if ((object$dist == "tneh" &
                object$link_tau == "loglinear")) {
      stop("Not implemented yet for this model")

    }

  } else if (object$model == "mixture") {
    if (object$dist == "weib") {
      cumLexctopred<-cumLexc_alphaweibull_topred(z_ucured =  z_ucured,
                                  z_pcured = z_pcured,
                                  x = time,
                                  theta, sign_delta = sign_delta)

      netsurv <- exp(-cumLexctopred$cumhaz)
      ex_haz <- lexc_alphaweibull(z_ucured =  z_ucured,
                                  z_pcured = z_pcured,
                                  x = time,
                                  theta,  sign_delta = sign_delta)

      pt_cure <- (cumLexctopred$cured) / netsurv

      cumHaze <- cumLexctopred$cumhaz
      log_cumHaze <- log(cumHaze)


      if(ncol(z_pcured)==0&ncol(z_ucured)==0){
        pre_pred <- data.frame(
          time = time,
          ex_haz = ex_haz,
          netsurv = netsurv,
          cumHaze = cumHaze,
          cured=cumLexctopred$cured,
          log_cumHaze = log_cumHaze,
          pt_cure = pt_cure
        )
      }else if(ncol(z_pcured)>0&ncol(z_ucured)==0){
        pre_pred <- data.frame(
          time = time,
          Z=z_pcured,
          ex_haz = ex_haz,
          netsurv = netsurv,
          cumHaze = cumHaze,
          cured=cumLexctopred$cured,
          log_cumHaze = log_cumHaze,
          pt_cure = pt_cure
        )
      }else if(ncol(z_pcured)==0&ncol(z_ucured)>0){
        pre_pred <- data.frame(
          time = time,
          Z=z_ucured,
          ex_haz = ex_haz,
          netsurv = netsurv,
          cumHaze = cumHaze,
          cured=cumLexctopred$cured,
          log_cumHaze = log_cumHaze,
          pt_cure = pt_cure
        )
      }else {
        rep_col<-intersect(colnames(z_pcured),colnames(z_ucured))
        unique_cols<-setdiff(colnames(z_ucured),rep_col)
        Z<-cbind(z_pcured,z_ucured[,unique_cols,drop=FALSE])
        pre_pred <- data.frame(
          time = time,
          Z=Z,
          ex_haz = ex_haz,
          netsurv = netsurv,
          cumHaze = cumHaze,
          cured=cumLexctopred$cured,
          log_cumHaze = log_cumHaze,
          pt_cure = pt_cure
        )
      }


      pred <- predcall_wei(object,
                           pre_pred,
                           z_pcured = z_pcured,
                           z_ucured = z_ucured,
                           x = time,
                           level = level,
                           epsilon = epsilon, sign_delta = sign_delta,
                           cumLexctopred=cumLexctopred)




    } else if (object$dist == "eweib") {
      stop("Not implemented yet for this model")

    }
  }

  class(pred) <- c("predCuresurv", "data.frame")
  return(pred)

}
