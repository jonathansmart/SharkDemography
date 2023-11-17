#' SharkDemography.
#'
#' @name SharkDemography.
#' @author Jonathan Smart
#' @description A package to conduct Shark demographic analyses using Leslie matrix models.
#' @docType package
NULL


utils::globalVariables(c(
  'rnorm','runif','txtProgressBar','setTxtProgressBar','quantile','Method','Value','parse_number','info','F.','mean.lambda',
  'AVG','.','Lambda','MaxAge', 'MinAge', 'Range', 'na.omit','.data'
))

#' Example life history data for Silky Sharks
#'
#' An example dataset of a completed set of life history data for Silky Sharks
#' in a `create_data_input()` template.
#' @docType data
#' @keywords datasets
#' @name Silky_data
#' @usage data(Silky_data)
#' @format A list of 10 with class "Demography.inputs"
NULL
