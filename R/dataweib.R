#' Simulated data with vital status information from Weibull mixture cure model
#'
#' Simulated data
#'
#'
#' @docType data
#'
#' @usage data(dataweib)
#'
#' @format This dataset contains the following variables:
#' \describe{
#'  \item{age}{Age at diagnosis}
#'  \item{age_cr}{centered and scaled age at diagnosis}
#'  \item{age_classe}{"<45", "45_59" and  ">=60" age groups }
#'  \item{sexe}{"male", "female" gender groups }
#'  \item{stage}{"<0", "1" , "2" and  "3" for stage I-IV groups }
#'  \item{time_obs}{Follow-up time (years)}
#'  \item{event}{Vital status}
#'  \item{cumehazard}{individual cumulative expected hazard}
#'  \item{ehazard}{individual instantaneous expected hazard}
#' }
#'
#'
#' @keywords datasets
#'
#'
#' @examples
#' data(dataweib)
#' summary(dataweib)
"dataweib"
