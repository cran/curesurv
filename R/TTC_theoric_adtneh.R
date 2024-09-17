#' @title TTC_theoric_adtneh function
#'
#' @description Product the value of TTC with TNEH model
#'
#' @param object An object of class \code{curesurv}.
#'
#' @param newdata newdata
#'
#' @keywords internal


TTC_theoric_adtneh <- function(newdata, object){

  if (!inherits(object, "curesurv"))
    stop("Primary argument much be a curesurv object")
  Call <- match.call()
  newcall <- Call[c(1, match(c("newdata"),
                             names(Call), nomatch = 0))]
  names(newcall)[2] <- "data"
  Terms <- newcall$formula <- object$Terms
  newcall[[1L]] <- quote(stats::model.frame)
  m_eval <- eval(newcall, parent.frame())
  event <- stats::model.extract(m_eval, "response")[, "status"]
  time <- stats::model.extract(m_eval, "response")[, "time"]
  myvarnames <- colnames(stats::model.matrix(Terms, m_eval)[,-1, drop = FALSE])
  z_alpha_id <- which(stringr::str_detect(c(myvarnames),
                                          pattern = "z_alpha"))
  z_tau_id <- which(stringr::str_detect(c(myvarnames),
                                        pattern = "z_tau"))

  if (length(z_alpha_id) > 0) {
    z_alpha <- as.data.frame(
      stats::model.matrix(Terms, m_eval)[,-1, drop = FALSE][, c(z_alpha_id)])
    colnames(z_alpha) <- c(stringr::str_remove(myvarnames[c(z_alpha_id)], "z_alpha"))
    z_alpha <- as.matrix(z_alpha)
  }
  else{
    z_alpha <- matrix(nrow = nrow(stats::model.matrix(Terms, m_eval)[,-1, drop = FALSE]),
                      ncol = 0)
  }

  if (length(z_tau_id) > 0) {
    z_tau <- as.data.frame(
      stats::model.matrix(Terms, m_eval)[,-1, drop = FALSE][, c(z_tau_id)])
    colnames(z_tau) <- c(stringr::str_remove(myvarnames[c(z_tau_id)],"z_tau"))
    z_tau <- as.matrix(z_tau)
  }
  else{
    z_tau <- matrix(nrow = nrow(stats::model.matrix(Terms, m_eval)[,-1, drop = FALSE]),
                    ncol = 0)
  }


  n_z_tau <- ncol(z_tau)
  n_z_alpha <- ncol(z_alpha)
  n_z_tau_ad <- n_z_tau - 1
  n_z_alpha_ad <- n_z_alpha - 1

  alpha0 <- 2.3
  alpha_k <- -0.1
  beta <- 4.8
  tau0 <- 5.5
  tau_k <- 0.9


  if (n_z_tau == 0 & n_z_alpha == 0) {

    alpha <- alpha0
    tau <- tau0

  }


  else if (n_z_tau > 0 & n_z_alpha > 0) {

    alpha <- alpha0 + z_alpha %*% alpha_k
    tau <- tau0 + z_tau %*% tau_k

  }


  else if (n_z_tau > 0 & n_z_alpha == 0) {

    alpha <- alpha0
    tau <- tau0 + z_tau %*% tau_k

  }


  else if (n_z_tau == 0 & n_z_alpha > 0) {

    alpha <- alpha0 + z_alpha %*% alpha_k
    tau <- tau0

  }


  funca <- function(x) {
    x ^ (alpha[1] - 1) * (1 - x) ^ (beta - 1)
  }
  funcb <-
    function(t) {
      exp(tau[1] * stats::integrate(funca, 1, t / tau[1])$value)
    }

  inverse <- function(f, lower, upper) {
    function(y)
      uniroot(function(x)
        f(x) - 0.95, lower = lower, upper = upper)[["root"]]
  }

  TTC <- inverse(funcb, 0, tau[1])
  time_to_cure <- TTC(funcb(1))


  return(time_to_cure)
}
