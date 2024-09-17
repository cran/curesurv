
dpcuredtheta <- function(newdata, object) {
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
  }else{
    z_alpha <- stats::model.matrix(Terms, m_eval)[,-1, drop = FALSE]
    if (length(z_tau_id) > 0)
      colnames(z_alpha) <- c(stringr::str_remove(myvarnames,"z_tau"))
  }

  if (length(z_tau_id) > 0) {
    z_tau <- as.data.frame(
      stats::model.matrix(Terms, m_eval)[,-1, drop = FALSE][, c(z_tau_id)])
    colnames(z_tau) <- c(stringr::str_remove(myvarnames[c(z_tau_id)],"z_tau"))
    z_tau <- as.matrix(z_tau)

  }else{
    z_tau <- stats::model.matrix(Terms, m_eval)[,-1, drop = FALSE]
    if (length(z_alpha_id) > 0)
      colnames(z_tau) <- c(stringr::str_remove(myvarnames, "z_alpha"))

  }

  theta <- object$coefficients
  n_z_tau <- ncol(z_tau)
  n_z_alpha <- ncol(z_alpha)
  n_z_tau_ad <- n_z_tau - 1
  n_z_alpha_ad <- n_z_alpha - 1
  alpha0 <- theta[1]
  if (n_z_tau > 0) {
    alpha_k <- theta[2:(n_z_alpha + 1)]
    beta <- (theta[n_z_alpha + 2])

    tau0 <- theta[n_z_tau + 3]
    tau_z <- theta[(n_z_tau + 4):(n_z_tau + 4 + n_z_tau_ad)]
    alpha <- exp(alpha0 + z_alpha %*% alpha_k)

    tau <- exp(tau0 + z_tau %*% tau_z)

  } else {

    beta <- (theta[2])

    tau0 <- theta[3]
    alpha <- exp(alpha0)
    tau <- exp(tau0)
  }
  beta2 <- (exp(beta - 1) + 1)
  D <- matrix(0, 1, (n_z_tau + 4 + n_z_tau_ad))
  aux <- -beta(alpha, beta2) * pcured(newdata, theta)
  D[, 1] <- aux * tau * (digamma(alpha) - digamma(alpha + beta2))
  D[, 2:(n_z_alpha + 1)] <- aux * tau * (digamma(alpha) - digamma(alpha + beta2)) * z_alpha#a1
  D[, (n_z_alpha + 2)] <- aux * tau * (digamma(beta2) - digamma(alpha + beta2))
  D[, (n_z_tau + 3)] <- aux
  D[, (n_z_tau + 4):(n_z_tau + 4 + n_z_tau_ad)] <- aux * z_tau# a2

  return(D)
  }
