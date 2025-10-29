#' Check integrity of soil layer profiles (contiguity, zero start, min length)
#'
#' @description
#' Validates that each (profile) defined by `group` columns has:
#' 1) layers sorted and strictly increasing,
#' 2) starts at 0 cm (if `require_zero_start = TRUE`),
#' 3) contiguity (`upper == lag(lower)`),
#' 4) maximum depth >= `min_core_length_cm`.
#'
#' @param data A data.frame/tibble with at least: `upper`, `lower` and grouping columns.
#' @param group Character vector of grouping columns defining a profile. Default: `c("trt","rep")`.
#' @param require_zero_start Logical; require the first layer to start at 0 cm. Default: `TRUE`.
#' @param min_core_length_cm Numeric; minimum required depth (cm). Default: `0`.
#'
#' @return A tibble with one row per profile and diagnostics:
#' * `n_layers`, `zero_start_ok`, `strict_increase_ok`, `contiguous_ok`,
#' * `max_lower`, `meets_min_core`, and `ok` (all checks passed).
#' Additionally, `first_break_at_row` indicates where contiguity breaks (if any).
#'
#' @examples
#' \dontrun{
#' diag <- check_profile_integrity(df, group = c("trt","rep"))
#' dplyr::filter(diag, !ok)
#' }
#'
#' @export
#' @importFrom dplyr arrange group_by ungroup summarise mutate across all_of n first last
#' @importFrom tibble tibble
check_profile_integrity <- function(data,
                                    group = c("trt","rep"),
                                    require_zero_start = TRUE,
                                    min_core_length_cm = 0) {
  if (!all(c("upper","lower") %in% names(data))) {
    stop("`data` must contain columns `upper` and `lower`.", call. = FALSE)
  }
  miss <- setdiff(group, names(data))
  if (length(miss)) {
    stop("Grouping columns not found in `data`: ", paste(miss, collapse = ", "), call. = FALSE)
  }

  # sort by group and depth
  d <- dplyr::arrange(data, dplyr::across(dplyr::all_of(group)), upper, lower)

  # per-profile diagnostics
  out <- d %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(group))) %>%
    dplyr::summarise(
      n_layers = dplyr::n(),
      first_upper = if (dplyr::n() > 0) upper[1] else NA_real_,
      last_lower  = if (dplyr::n() > 0) lower[dplyr::n()] else NA_real_,
      zero_start_ok = if (dplyr::n() > 0) (!require_zero_start || isTRUE(all.equal(first_upper, 0))) else FALSE,
      # strictly increasing:
      strict_increase_ok = if (dplyr::n() > 1) all(diff(lower) > 0 & diff(upper) >= 0) else TRUE,
      # contiguity: upper[i] == lower[i-1]
      contiguous_ok = if (dplyr::n() > 1) all(upper[-1] == lower[-length(lower)]) else TRUE,
      # where contiguity first breaks
      first_break_at_row = {
        if (dplyr::n() > 1) {
          brk <- which(upper[-1] != lower[-length(lower)])
          if (length(brk)) brk[1] + 1 else NA_integer_
        } else NA_integer_
      },
      max_lower = if (dplyr::n() > 0) max(lower, na.rm = TRUE) else NA_real_,
      meets_min_core = is.finite(max_lower) && !is.na(max_lower) &&
        max_lower >= min_core_length_cm,
      ok = zero_start_ok && strict_increase_ok && contiguous_ok && meets_min_core,
      .groups = "drop"
    )

  out
}
