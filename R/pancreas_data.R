#' Simulated pancreas data with vital status information
#'
#' Simulated data
#'
#'
#' @docType data
#'
#' @usage data(pancreas_data)
#'
#' @format This dataset contains the following variables:
#' \describe{
#'  \item{age}{Age at diagnosis}
#'  \item{age_cr}{centered and scaled age at diagnosis}
#'  \item{age_classe}{"<45", "45_59" and  ">=60" age groups }
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
#' data(pancreas_data)
#' summary(pancreas_data)
"pancreas_data"
