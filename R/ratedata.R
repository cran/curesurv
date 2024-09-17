ratedata <- function(pophaz = pophaz, cumpophaz = cumpophaz, data = data) {
  if (is.null(pophaz) & is.null(cumpophaz)) {
  pophaz <- rep(0, nrow(data))
  cumpophaz <- rep(0, nrow(data))

}else if (!is.null(pophaz) & is.null(cumpophaz)) {
  pophaz <- data[,eval(pophaz)]
  cumpophaz <- rep(0, nrow(data))

}else if (!is.null(pophaz) & !is.null(cumpophaz)) {
  pophaz <- data[,eval(pophaz)]
  cumpophaz <- data[,eval(cumpophaz)]

}else if (is.null(pophaz) & !is.null(cumpophaz)) {
  stop("please provide pophaz variable in the data")
}

  rate_data <- data.frame(pophaz, cumpophaz)
  return(rate_data)
}
