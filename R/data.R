
#' Example NHANES Data for Survey-Weighted Cox Spline Analyses
#'
#' @description
#' Example data included with the svycoxspline package for
#' demonstrating survey-weighted Cox proportional hazards
#' models using natural cubic splines.
#'
#' @details
#' The dataset includes variables from the National Health
#' and Nutrition Examination Survey (NHANES) 2011-2012 and 2013-2014, including
#' survey design variables, body mass index, follow-up
#' time, mortality status, age, and sex.
#'
#' Survey design variables include:
#' \itemize{
#'   \item \code{sdmvpsu}: Masked variance pseudo-PSU.
#'   \item \code{sdmvstra}: Masked variance pseudo-stratum.
#'   \item \code{mec4yr}: Four-year MEC examination weight.
#' }
#'
#' The data are provided for demonstrating the functions
#' and methods implemented in the svycoxspline package.
#'
#' The package vignette loads the example CSV using:
#'
#' \preformatted{
#' dat <- data.table::fread(
#'   system.file(
#'     "data",
#'     "nhanes_spline_data.csv",
#'     package = "svycoxspline"
#'   ),
#'   data.table = FALSE
#' )
#' }
#'
#' @format An example dataset containing NHANES survey
#'   design variables, BMI, follow-up time, mortality
#'   status, age, and sex.
#'
#' @source
#' National Health and Nutrition Examination Survey
#' (NHANES), National Center for Health Statistics,
#' Centers for Disease Control and Prevention.
#'
#' Original NHANES datasets and documentation:
#' \url{https://wwwn.cdc.gov/nchs/nhanes/}
#'
#' NHANES Linked Mortality Files:
#' \url{https://www.cdc.gov/nchs/data-linkage/mortality-public.htm}
#'
#' @name nhanes_spline_data
#' @docType data
#' @keywords datasets
NULL
