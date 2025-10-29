#' Aggregate ESM or FD results by treatment and depth
#'
#' @description
#' Summarises the outputs from [esm()] (either ESM or fixed-depth) by treatment,
#' depth, and optionally replicate or site. This helper function computes means,
#' standard deviations, and standard errors for stocks or properties (e.g., SOC, N, P).
#'
#' @param data A data frame produced by [esm()].
#' @param var Character string. The response variable name (e.g., `"SOC"`).
#' @param group_by Optional character vector with columns to aggregate by.
#'   Defaults to `c("trt", "lower")`.
#' @param se Logical. If `TRUE`, computes standard error of the mean. Default = `TRUE`.
#' @param na.rm Logical. Remove `NA` before summarising. Default = `TRUE`.
#'
#' @details
#' The function automatically identifies columns ending with `_g_cm2` or
#' `Cum_..._g_cm2` to aggregate in the same units.
#'
#' @return
#' A tibble with mean ± sd (and se if requested) by the specified groups.
#'
#' @seealso [esm()], [as_Mg_ha()]
#'
#' @examples
#' \dontrun{
#' esm <- esm(df, som = FALSE, soc_col = "SOC", rvar = "SOC", output = "esm")
#' esm_aggregate(esm, var = "SOC")
#' }
#'
#' @export
#' @importFrom dplyr group_by summarise across all_of ungroup n
esm_aggregate <- function(data, var, group_by = c("trt", "lower"), se = TRUE, na.rm = TRUE) {

  if (!requireNamespace("dplyr", quietly = TRUE))
    stop("Package 'dplyr' is required.", call. = FALSE)

  # select relevant columns dynamically
  pattern <- paste0("^(", var, "|", var, "_g_cm2|Cum_", var, "_g_cm2)$")
  vars_to_summarise <- grep(pattern, names(data), value = TRUE)
  if (length(vars_to_summarise) == 0) stop("No columns matching '", var, "' found.")

  out <- data |>
    dplyr::group_by(dplyr::across(dplyr::all_of(group_by))) |>
    dplyr::summarise(
      dplyr::across(
        dplyr::all_of(vars_to_summarise),
        list(
          mean = ~mean(.x, na.rm = na.rm),
          sd   = ~sd(.x, na.rm = na.rm),
          se   = ~if (se) sd(.x, na.rm = na.rm)/sqrt(dplyr::n()) else NA_real_
        ),
        .names = "{.col}_{.fn}"
      ),
      .groups = "drop"
    )

  out
}
