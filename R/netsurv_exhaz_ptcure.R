netsurv_exhaz_ptcure <- function(z_pcured, time = time, theta) {
  netsurv <- exp(-cumLexc_alphaweibull_topred(z_pcured,
                                            x = time,
                                            theta)$cumhaz)
  ex_haz <- lexc_alphaweibull(z_pcured,
                            x = time,
                            theta)

  pt_cure <- cumLexc_alphaweibull_topred(z_pcured,
                                       x = time,
                                       theta)$cured / netsurv

  res = c(netsurv, ex_haz, pt_cure)

return(res)
}
