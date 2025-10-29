#' Assert integrity of soil profiles and fail with a readable message
#'
#' @description
#' Runs [check_profile_integrity()] and stops with an informative error
#' if any profile fails. Useful at the top of pipelines or tests.
#'
#' @inheritParams check_profile_integrity
#'
#' @return Invisibly returns the diagnostics tibble when all profiles pass.
#' @export
assert_profile_integrity <- function(data,
                                     group = c("trt","rep"),
                                     require_zero_start = TRUE,
                                     min_core_length_cm = 0) {
  diag <- check_profile_integrity(data, group, require_zero_start, min_core_length_cm)
  bad <- diag[!diag$ok, , drop = FALSE]
  if (nrow(bad) > 0) {
    # build a concise message per bad profile
    msgs <- apply(bad, 1, function(row) {
      # row is a named vector; extract flags
      flags <- c(
        if (!isTRUE(as.logical(row[["zero_start_ok"]]))) "zero_start" else NULL,
        if (!isTRUE(as.logical(row[["strict_increase_ok"]]))) "strict_increase" else NULL,
        if (!isTRUE(as.logical(row[["contiguous_ok"]]))) paste0("contiguity@row ", row[["first_break_at_row"]]) else NULL,
        if (!isTRUE(as.logical(row[["meets_min_core"]]))) "min_core" else NULL
      )
      grp_vals <- paste(sprintf("%s=%s", names(row)[seq_along(group)], row[seq_along(group)]), collapse = ", ")
      paste0("{ ", grp_vals, " } failed: ", paste(flags, collapse = ", "))
    })
    stop("Profile integrity check failed for:\n - ", paste(msgs, collapse = "\n - "),
         call. = FALSE)
  }
  invisible(diag)
}
