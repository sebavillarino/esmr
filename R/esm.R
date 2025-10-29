#' Equivalent Soil Mass (ESM) for soil stocks and properties
#'
#' @description
#' Computes per-interval and cumulative quantities for a response variable
#' (e.g., %C, %N, pH) using the Equivalent Soil Mass (ESM) approach, or
#' at Fixed Depth (FD).
#'
#' @details
#' Expected units:
#' * `upper`, `lower`: cm
#' * `bd`: g cm^-3 (lowercase). If your data has `BD`, it is mapped to `bd`.
#' * `som`: percent by mass (%) when `som = TRUE`.
#' * When `som = FALSE`, you must provide `soc_col` (percent by mass, %), and
#'   SOM is estimated as `2 * SOC`.
#' * `rvar`: percent by mass (%) if elemental concentration; non-additive props
#'   (e.g., pH) are re-interpolated onto the SMM grid (no "pH stock").
#'
#' @param data Data frame with: `trt`, `rep`, `upper`, `lower`, `bd`, and `rvar`.
#'   Optionally `som` (if `som = TRUE`) or a SOC column named in `soc_col` (if `som = FALSE`).
#'   `ref_trt` is only required when `output = "esm"` and
#'   `reference_mode = "treatment_mean"`.
#' @param rvar String. Name of the response variable column (percent).
#' @param som Logical. If `TRUE`, use provided `som` column; if `FALSE`, estimate
#'   SOM as `2 * SOC` using the column named in `soc_col`. Default: `TRUE`.
#' @param soc_col String or `NULL`. Required when `som = FALSE`. Name of SOC (%) column.
#' @param min_core_length_cm Numeric. Minimum core length (cm). Default 0.
#' @param output "esm" (default) or "fd".
#' @param reference_mode "treatment_mean" (default), "min_smm_by_depth", or "custom".
#' @param custom_ld,custom_smm Custom reference grid when `reference_mode = "custom"`.
#' @param extrapolation "none" (default) or "to_depth".
#' @param reference_by Optional character vector of column names to build the
#'   reference grid by group (e.g., "Site", or c("Site","yr")).
#'
#' @return Tibble with per-interval and cumulative outputs (always with `bd` in lowercase).
#'   For `output="esm"`, the attribute `"SMM_reference"` contains the SMM target grid used
#'   (including `reference_by` columns if provided).
#' @export
esm <- function(data,
                rvar,
                som = TRUE,
                soc_col = NULL,
                min_core_length_cm = 0,
                output = c("esm", "fd"),
                reference_mode = c("treatment_mean", "min_smm_by_depth", "custom"),
                custom_ld = NULL,
                custom_smm = NULL,
                extrapolation = c("none", "to_depth"),
                reference_by = NULL) {

  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("Package 'dplyr' is required.", call. = FALSE)
  }

  output <- match.arg(output)
  reference_mode <- match.arg(reference_mode)
  extrapolation <- match.arg(extrapolation)

  if (!is.data.frame(data)) stop("`data` must be a data.frame/tibble.", call. = FALSE)

  # ---- standardize bd to lowercase ----
  if (!"bd" %in% names(data)) {
    if ("BD" %in% names(data)) {
      data <- dplyr::rename(data, bd = BD)
    } else {
      stop("Missing column: `bd` (g cm^-3).", call. = FALSE)
    }
  }

  # ---- validate reference_by if used ----
  if (!is.null(reference_by)) {
    if (!is.character(reference_by)) stop("`reference_by` must be character vector.", call. = FALSE)
    miss <- setdiff(reference_by, names(data))
    if (length(miss)) stop("`reference_by` columns not in data: ", paste(miss, collapse = ", "), call. = FALSE)
  }

  # ---- conditional column requirements ----
  req_base <- c("trt", "rep", "upper", "lower", "bd", rvar)
  missing_cols <- setdiff(req_base, names(data))
  if (length(missing_cols) > 0) {
    stop(paste("Missing columns:", paste(missing_cols, collapse = ", ")), call. = FALSE)
  }

  need_ref_trt <- (output == "esm" && reference_mode == "treatment_mean")
  if (need_ref_trt) {
    if (!"ref_trt" %in% names(data)) stop("Missing column: `ref_trt` for treatment_mean ESM.", call. = FALSE)
    if (any(is.na(data$ref_trt)))     stop("NAs in `ref_trt`.", call. = FALSE)
  }

  # ---- numeric checks ----
  num_ok <- sapply(data[, c("upper", "lower", "bd", rvar)], is.numeric)
  if (!all(num_ok)) stop("`upper`, `lower`, `bd`, and `rvar` must be numeric.", call. = FALSE)
  if (min_core_length_cm < 0) stop("`min_core_length_cm` must be >= 0.", call. = FALSE)

  # ---- SOM handling ----
  if (isTRUE(som)) {
    if (!"som" %in% names(data)) stop("`som = TRUE` but `som` column is missing.", call. = FALSE)
    if (!is.numeric(data$som))   stop("`som` must be numeric (%).", call. = FALSE)
    SOM_vec <- data$som
  } else {
    if (is.null(soc_col) || !(soc_col %in% names(data))) {
      stop("`som = FALSE` requires a valid `soc_col` present in `data`.", call. = FALSE)
    }
    if (!is.numeric(data[[soc_col]])) stop("`soc_col` must be numeric (%).", call. = FALSE)
    if (any(is.na(data[[soc_col]])))  stop("`soc_col` contains NAs; cannot estimate SOM.", call. = FALSE)
    SOM_vec <- 2 * as.numeric(data[[soc_col]])
  }

  # ---- build working df ----
  keep_cols <- c("trt", "rep", "upper", "lower", "bd", rvar,
                 intersect("ref_trt", names(data)),
                 reference_by)
  d <- data[, keep_cols, drop = FALSE]
  if (!"ref_trt" %in% names(d)) d$ref_trt <- NA_character_
  names(d)[names(d) == rvar] <- "RVAR"
  d$SOM  <- SOM_vec
  d$Type <- "FD"

  # must start at 0 cm
  d <- dplyr::group_by(d, trt, rep) |>
    dplyr::filter(min(upper, na.rm = TRUE) == 0) |>
    dplyr::ungroup()

  # drop NAs in key fields
  d <- dplyr::filter(d, !is.na(upper) & !is.na(lower) &
                       !is.na(SOM) & !is.na(bd) & !is.na(RVAR))

  # order & enforce layer contiguity
  d <- dplyr::arrange(d, trt, rep, upper, lower)
  contiguous <- function(u, l) all(u == dplyr::lag(l, default = 0))
  d <- dplyr::group_by(d, trt, rep) |>
    dplyr::filter(contiguous(upper, lower)) |>
    dplyr::ungroup()

  # ---- enforce minimum core length (robust: summarise & inner_join) ----
  max_depth <- d %>%
    dplyr::group_by(trt, rep) %>%
    dplyr::summarise(
      max_lower = if (all(is.na(lower))) NA_real_ else max(lower, na.rm = TRUE),
      .groups = "drop"
    )

  d <- d %>%
    dplyr::inner_join(
      max_depth %>%
        dplyr::filter(!is.na(max_lower) & is.finite(max_lower) &
                        max_lower >= min_core_length_cm),
      by = c("trt", "rep")
    )

  if (nrow(d) == 0) {
    stop("No observations remain after filters (start at 0 cm, contiguity, min core length).", call. = FALSE)
  }

  # add 0-mass row at (0,0) for initial interpolation
  zeros <- dplyr::distinct(d, dplyr::across(c(trt, rep, ref_trt, dplyr::all_of(reference_by))))
  zeros$Type <- "FD"; zeros$upper <- 0; zeros$lower <- 0; zeros$SOM <- 0; zeros$bd <- 0; zeros$RVAR <- 0
  d <- dplyr::bind_rows(d, zeros) |>
    dplyr::arrange(trt, rep, upper, lower)

  # ---- fixed-depth per-interval calculations ----
  d$Soil_g_cm2     <- (d$lower - d$upper) * d$bd
  d$RVAR_g_cm2     <- (d$RVAR / 100) * d$Soil_g_cm2
  d$SOM_g_cm2      <- (d$SOM  / 100) * d$Soil_g_cm2
  d$Min_Soil_g_cm2 <- d$Soil_g_cm2 - d$SOM_g_cm2

  # cumulative (FD)
  d_fd <- dplyr::group_by(d, trt, rep) |>
    dplyr::mutate(
      Cum_Soil_g_cm2     = cumsum(Soil_g_cm2),
      Cum_RVAR_g_cm2     = cumsum(RVAR_g_cm2),
      Cum_SOM_g_cm2      = cumsum(SOM_g_cm2),
      Cum_Min_Soil_g_cm2 = cumsum(Min_Soil_g_cm2)
    ) |>
    dplyr::ungroup()

  # remove (0,0)
  d_fd <- dplyr::filter(d_fd, !(upper == 0 & lower == 0))
  d_fd$Type <- "FD"

  # ---- fast return for FD ----
  if (output == "fd") {
    out <- d_fd
    names(out)[names(out) == "RVAR"]       <- rvar
    names(out)[names(out) == "RVAR_g_cm2"] <- paste0(rvar, "_g_cm2")
    return(dplyr::relocate(out, Type, trt, rep, ref_trt))
  }

  # ---- build ESM reference grid ----
  if (reference_mode == "custom") {
    if (is.null(custom_ld) || is.null(custom_smm)) {
      stop("`reference_mode = 'custom'` requires `custom_ld` and `custom_smm`.", call. = FALSE)
    }
    if (length(custom_ld) != length(custom_smm)) {
      stop("`custom_ld` and `custom_smm` must have equal length.", call. = FALSE)
    }
    if (is.unsorted(custom_ld, strictly = TRUE)) {
      stop("`custom_ld` must be strictly increasing.", call. = FALSE)
    }
    ref_tbl <- dplyr::tibble(
      upper = c(0, head(custom_ld, -1)),
      lower = custom_ld,
      Cum_Min_Soil_g_cm2 = custom_smm
    )
    if (!is.null(reference_by)) {
      for (nm in reference_by) ref_tbl[[nm]] <- NA
    }
    ref_tbl$ref_trt <- "custom"

  } else if (reference_mode == "min_smm_by_depth") {
    # EXACT: FD → min by `lower` (optionally by reference_by) → custom grid
    group_vars <- c(reference_by, "lower")

    min_smm <- d_fd |>
      dplyr::group_by(dplyr::across(dplyr::all_of(group_vars))) |>
      dplyr::summarise(
        Cum_Min_Soil_g_cm2 = min(Cum_Min_Soil_g_cm2, na.rm = TRUE),
        .groups = "drop"
      ) |>
      dplyr::arrange(dplyr::across(dplyr::all_of(reference_by)), lower)

    ref_tbl <- min_smm |>
      dplyr::group_by(dplyr::across(dplyr::all_of(reference_by))) |>
      dplyr::mutate(upper = dplyr::lag(lower, default = 0)) |>
      dplyr::ungroup() |>
      dplyr::select(dplyr::all_of(reference_by), upper, lower, Cum_Min_Soil_g_cm2)

    # (Optional) enforce monotonic SMM per group:
    # ref_tbl <- ref_tbl |>
    #   dplyr::group_by(dplyr::across(dplyr::all_of(reference_by))) |>
    #   dplyr::mutate(Cum_Min_Soil_g_cm2 = cummax(Cum_Min_Soil_g_cm2)) |>
    #   dplyr::ungroup()

    ref_tbl$ref_trt <- "min_smm_by_depth"

  } else { # "treatment_mean"
    ref_rows <- dplyr::filter(d_fd, trt == ref_trt)
    if (nrow(ref_rows) == 0) stop("No rows for the reference treatment found in `data`.", call. = FALSE)

    ref_depths <- dplyr::distinct(dplyr::select(ref_rows, upper, lower))
    ref_tbl <- ref_depths |>
      dplyr::left_join(
        ref_rows |>
          dplyr::group_by(upper, lower) |>
          dplyr::summarise(Cum_Min_Soil_g_cm2 = mean(Cum_Min_Soil_g_cm2, na.rm = TRUE), .groups = "drop"),
        by = c("upper", "lower")
      )
    if (!is.null(reference_by)) {
      for (nm in reference_by) ref_tbl[[nm]] <- NA
    }
    ref_tbl$ref_trt <- unique(d_fd$ref_trt)
  }

  # ---- ESM interpolation by (trt, rep) ----
  by_keys <- dplyr::distinct(dplyr::select(d_fd, trt, rep, ref_trt, dplyr::all_of(reference_by)))
  make_mono_spline <- function(x, y) stats::splinefun(x, y, method = "monoH.FC")

  esm_list <- lapply(seq_len(nrow(by_keys)), function(i) {
    key <- by_keys[i, ]
    cur <- dplyr::filter(d_fd, trt == key$trt, rep == key$rep)

    # pick reference rows from the corresponding group (if any)
    ref_subset <- if (is.null(reference_by)) ref_tbl else {
      out <- ref_tbl
      for (nm in reference_by) out <- dplyr::filter(out, get(nm) == key[[nm]])
      out
    }

    # apply extrapolation rule
    if (identical(extrapolation, "none")) {
      ref_use <- dplyr::filter(ref_subset, Cum_Min_Soil_g_cm2 <= max(cur$Cum_Min_Soil_g_cm2, na.rm = TRUE))
    } else { # "to_depth"
      ref_use <- dplyr::filter(ref_subset, lower <= max(cur$lower, na.rm = TRUE))
    }
    if (nrow(ref_use) == 0) return(NULL)

    s_r   <- make_mono_spline(cur$Cum_Min_Soil_g_cm2, cur$Cum_RVAR_g_cm2)
    s_som <- make_mono_spline(cur$Cum_Min_Soil_g_cm2, cur$Cum_SOM_g_cm2)

    ref_use$Cum_RVAR_g_cm2 <- s_r(ref_use$Cum_Min_Soil_g_cm2)
    ref_use$Cum_SOM_g_cm2  <- s_som(ref_use$Cum_Min_Soil_g_cm2)

    ref_use <- ref_use |>
      dplyr::mutate(
        Cum_RVAR_g_cm2 = cummax(Cum_RVAR_g_cm2),
        Cum_SOM_g_cm2  = cummax(Cum_SOM_g_cm2)
      )

    out <- ref_use |>
      dplyr::mutate(
        trt     = key$trt,
        ref_trt = key$ref_trt,
        rep     = key$rep,
        Min_Soil_g_cm2 = Cum_Min_Soil_g_cm2 - dplyr::lag(Cum_Min_Soil_g_cm2, default = 0),
        RVAR_g_cm2     = Cum_RVAR_g_cm2     - dplyr::lag(Cum_RVAR_g_cm2, default = 0),
        SOM_g_cm2      = Cum_SOM_g_cm2      - dplyr::lag(Cum_SOM_g_cm2, default = 0),
        Soil_g_cm2     = Min_Soil_g_cm2 + SOM_g_cm2,
        bd             = Soil_g_cm2 / (lower - upper),
        RVAR           = ifelse(Soil_g_cm2 > 0, RVAR_g_cm2 / Soil_g_cm2 * 100, NA_real_),
        SOM            = ifelse(Soil_g_cm2 > 0, SOM_g_cm2  / Soil_g_cm2 * 100, NA_real_),
        Type           = "ESM"
      ) |>
      dplyr::mutate(
        Cum_Soil_g_cm2 = cumsum(Soil_g_cm2),
        Cum_RVAR_g_cm2 = cumsum(RVAR_g_cm2),
        Cum_SOM_g_cm2  = cumsum(SOM_g_cm2)
      )

    out
  })

  esm_out <- dplyr::bind_rows(esm_list)
  if (nrow(esm_out) == 0) {
    stop("Failed to build ESM series for any (trt, rep) given the reference and extrapolation mode.", call. = FALSE)
  }

  # ---- column order + names (bd lowercase) ----
  base_cols <- c("Type", "trt", "rep", "ref_trt", reference_by, "upper", "lower",
                 "RVAR", "SOM", "bd",
                 "Soil_g_cm2", "SOM_g_cm2", "Min_Soil_g_cm2", "RVAR_g_cm2",
                 "Cum_Soil_g_cm2", "Cum_RVAR_g_cm2", "Cum_SOM_g_cm2", "Cum_Min_Soil_g_cm2")
  base_cols <- intersect(base_cols, names(esm_out))
  esm_out <- esm_out[, base_cols]

  names(esm_out)[names(esm_out) == "RVAR"]            <- rvar
  names(esm_out)[names(esm_out) == "RVAR_g_cm2"]      <- paste0(rvar, "_g_cm2")
  names(esm_out)[names(esm_out) == "Cum_RVAR_g_cm2"]  <- paste0("Cum_", rvar, "_g_cm2")

  # attach reference grid (with reference_by columns if present)
  ref_cols <- c(reference_by, "upper", "lower", "Cum_Min_Soil_g_cm2", "ref_trt")
  ref_cols <- intersect(ref_cols, names(ref_tbl))
  attr(esm_out, "SMM_reference") <- ref_tbl[, ref_cols, drop = FALSE]

  dplyr::arrange(esm_out, trt, rep, upper, lower)
}
