#' Convert g/cm² to Mg/ha
#'
#' @description
#' Converts areal mass units from grams per square centimeter (`g cm^-2`)
#' to megagrams per hectare (`Mg ha^-1`), as commonly used for SOC or nutrient stocks.
#'
#' @param x Numeric vector in g/cm².
#'
#' @details
#' Conversion factor:
#' 1 g/cm² = 100 Mg/ha
#'
#' This function is unit-agnostic; it simply multiplies by 100.
#'
#' @return Numeric vector in Mg/ha.
#'
#' @examples
#' as_Mg_ha(2.5)   # returns 250
#'
#' @export
as_Mg_ha <- function(x) {
  if (!is.numeric(x))
    stop("`x` must be numeric.", call. = FALSE)
  x * 100
}
