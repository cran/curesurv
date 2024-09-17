gen_start_values <- function(nvalues,
                             npar,
                             min, max,
                             seed = 1980) {
  set.seed(seed)
  start <- randtoolbox::sobol(n = nvalues, dim = npar)
  rescale_start <- sapply(1:ncol(start),
                          function(i) {
                            min[i] +
                              (max[i] - min[i]) * start[,i]
                          })
  return(rescale_start)
}


