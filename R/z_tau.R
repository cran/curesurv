#' @title z_tau function identifying variables acting on tau parameter
#'
#' @description variables adjusted on tau parameter in
#'  non-mixture cure model with "tneh" specified for the distribution.
#'
#' @param x the name of the column in the dataset representing the variable that will act on tau parameter of the "tneh" model
#'
#' @keywords z_tau
#'
#' @return the variable x
#'
#' @author Juste Goungounga, Judith Breaud, Eugenie Blandin, Olayide Boussari, Valerie Jooste
#'
#' @references Boussari O, Bordes L, Romain G, Colonna M, Bossard N,
#'             Remontet L, Jooste V. Modeling excess hazard with time-to-cure
#'             as a parameter. Biometrics. 2020 Aug 31.
#'             doi: 10.1111/biom.13361. Epub ahead of print.
#'             PMID: 32869288.
#' (\href{https://pubmed.ncbi.nlm.nih.gov/32869288/}{pubmed})
#' @export
z_tau <- function(x){
  m.eval <- (quote(x))
  return(eval(m.eval))
}

