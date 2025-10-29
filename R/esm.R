#' Equivalent Soil Mass (ESM) for soil stocks and properties
#'
#' @description
#' Computes per-interval and cumulative quantities for a response variable
#' (e.g., \%C, \%N, pH expressed as %) using the Equivalent Soil Mass (ESM)
#' approach, or returns Fixed-Depth (FD) quantities.
#'
#' @param data Data frame with at least: `trt`, `rep`, `upper`, `lower`, `bd`, and the column named in `rvar`.
#'             If `som = TRUE`, also requires a `som` column.
#' @param rvar String. Name of the response variable column (percent).
#' @param som Logical. If `TRUE`, use provided `som` column; if `FALSE`, estimate SOM as `2 * SOC` using `soc_col`. Default `TRUE`.
#' @param soc_col String or `NULL`. Required when `som = FALSE`. Name of SOC (%) column used to estimate SOM.
#' @param min_core_length_cm Numeric. Minimum core length (cm). Default `0`.
#' @param output `"esm"` (default) or `"fd"`.
#' @param reference_mode `"ref_trt"`, `"min"`, or `"custom"`.
#'        Backward compatible alias: `"min_smm_by_depth"` → `"min"`.
#' @param custom_ld,custom_smm Custom reference grid when `reference_mode = "custom"`.
#'        `custom_ld` must be strictly increasing; both vectors must have equal length.
#' @param extrapolation `"none"` (default) or `"to_depth"`.
#'   - `"none"`: never extrapolate beyond the maximum *sampled mineral mass* of each profile.
#'   - `"to_depth"`: allow extrapolation up to the maximum *sampled depth* (lower) of each profile.
#' @param reference_by Optional character vector of columns that define groups **over which the reference is computed**.
#'        If `NULL`, `"min"` reference is global by depth (entire dataset).
#' @param core_by Optional character vector of columns that identify **unique cores/profiles** in addition to `trt` and `rep`
#'        (e.g., `"Site"`, `"Plot"`, `"Block"`). If `NULL`, only `trt + rep` identify a core.
#' @param ... Ignored. Accepts legacy `group_by` for backward compatibility (mapped to `reference_by`).
#'
#' @return A tibble with per-interval and cumulative outputs. Column `bd` is always lowercase.
#'         For `output="esm"`, the attribute `"SMM_reference"` contains the reference grid used.
#' @export
#' @importFrom dplyr %>% group_by ungroup filter arrange mutate lag distinct select left_join summarise
#' @importFrom dplyr across bind_rows relocate inner_join n
#' @importFrom stats splinefun
esm <- function(data,
                rvar,
                som = TRUE,
                soc_col = NULL,
                min_core_length_cm = 0,
                output = c("esm", "fd"),
                reference_mode = c("ref_trt", "min", "custom", "min_smm_by_depth"),
                custom_ld = NULL,
                custom_smm = NULL,
                extrapolation = c("none", "to_depth"),
                reference_by = NULL,
                core_by = NULL,
                ...) {

  output <- match.arg(output)
  reference_mode <- match.arg(reference_mode, c("ref_trt", "min", "custom", "min_smm_by_depth"))
  extrapolation  <- match.arg(extrapolation)

  # Back-compat aliases
  if (identical(reference_mode, "min_smm_by_depth")) {
    message("`reference_mode = 'min_smm_by_depth'` is deprecated; using `reference_mode = 'min'`.")
    reference_mode <- "min"
  }
  dots <- list(...)
  if (is.null(reference_by) && !is.null(dots$group_by)) {
    message("`group_by` is deprecated; using its value as `reference_by`.")
    reference_by <- dots$group_by
  }

  if (!is.data.frame(data)) stop("`data` must be a data.frame/tibble.", call. = FALSE)

  # Enforce bd lowercase
  if (!"bd" %in% names(data)) {
    if ("BD" %in% names(data)) {
      data <- dplyr::rename(data, bd = BD)
    } else {
      stop("Missing column: `bd` (g cm^-3).", call. = FALSE)
    }
  }

  # Required columns
  req_base <- c("trt", "rep", "upper", "lower", "bd", rvar)
  miss <- setdiff(req_base, names(data))
  if (length(miss)) stop("Missing columns: ", paste(miss, collapse = ", "), call. = FALSE)

  need_ref_trt <- (output == "esm" && reference_mode == "ref_trt")
  if (need_ref_trt) {
    if (!"ref_trt" %in% names(data)) stop("Missing `ref_trt` for treatment_mean ESM.", call. = FALSE)
    if (any(is.na(data$ref_trt)))     stop("NAs in `ref_trt`.", call. = FALSE)
  } else {
    if (!"ref_trt" %in% names(data)) data$ref_trt <- NA_character_
  }

  if (!is.null(reference_by)) {
    if (!is.character(reference_by)) stop("`reference_by` must be a character vector.", call. = FALSE)
    miss_g <- setdiff(reference_by, names(data))
    if (length(miss_g)) stop("`reference_by` columns not in data: ", paste(miss_g, collapse = ", "), call. = FALSE)
  }

  if (!is.null(core_by)) {
    if (!is.character(core_by)) stop("`core_by` must be a character vector.", call. = FALSE)
    miss_c <- setdiff(core_by, names(data))
    if (length(miss_c)) stop("`core_by` columns not in data: ", paste(miss_c, collapse = ", "), call. = FALSE)
  }

  # Compose core identifier columns (no autodetect)
  core_vars <- unique(c(core_by, "trt", "rep"))

  # Numeric checks
  num_ok <- sapply(data[, c("upper", "lower", "bd", rvar)], is.numeric)
  if (!all(num_ok)) stop("`upper`, `lower`, `bd`, and `rvar` must be numeric.", call. = FALSE)
  if (min_core_length_cm < 0) stop("`min_core_length_cm` must be >= 0.", call. = FALSE)

  # SOM handling
  if (isTRUE(som)) {
    if (!"som" %in% names(data)) stop("`som = TRUE` but `som` column is missing.", call. = FALSE)
    if (!is.numeric(data$som))   stop("`som` must be numeric (%).", call. = FALSE)
    SOM_vec <- data$som
  } else {
    if (is.null(soc_col) || !(soc_col %in% names(data))) {
      stop("`som = FALSE` requires `soc_col` present in `data`.", call. = FALSE)
    }
    if (!is.numeric(data[[soc_col]])) stop("`soc_col` must be numeric (%).", call. = FALSE)
    if (any(is.na(data[[soc_col]])))  stop("`soc_col` contains NAs; cannot estimate SOM.", call. = FALSE)
    SOM_vec <- 2 * as.numeric(data[[soc_col]])
  }

  # Working df
  keep_cols <- unique(c(core_vars, "ref_trt", "upper", "lower", "bd", rvar, reference_by))
  d <- data[, keep_cols, drop = FALSE]
  names(d)[names(d) == rvar] <- "RVAR"
  d$SOM  <- SOM_vec
  d$Type <- "FD"

  # Must start at 0 cm
  d <- d %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(core_vars))) %>%
    dplyr::filter(min(upper, na.rm = TRUE) == 0) %>%
    dplyr::ungroup()

  # Drop NAs in key fields
  d <- dplyr::filter(d, !is.na(upper) & !is.na(lower) &
                       !is.na(SOM) & !is.na(bd) & !is.na(RVAR))

  # Order & contiguity
  d <- d %>% dplyr::arrange(dplyr::across(dplyr::all_of(core_vars)), upper, lower)
  contiguous <- function(u, l) all(u == dplyr::lag(l, default = 0))
  d <- d %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(core_vars))) %>%
    dplyr::filter(contiguous(upper, lower)) %>%
    dplyr::ungroup()

  # Enforce minimum core length
  max_depth <- d %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(core_vars))) %>%
    dplyr::summarise(
      max_lower = if (all(is.na(lower))) NA_real_ else max(lower, na.rm = TRUE),
      .groups = "drop"
    )

  d <- d %>%
    dplyr::inner_join(
      max_depth %>%
        dplyr::filter(!is.na(max_lower) & is.finite(max_lower) &
                        max_lower >= min_core_length_cm),
      by = core_vars
    )

  if (nrow(d) == 0) {
    stop("No observations remain after filters (start at 0 cm, contiguity, min core length).", call. = FALSE)
  }

  # Add 0-mass row at (0,0) per core
  zero_keys <- unique(c(core_vars, "ref_trt"))
  zeros <- d %>%
    dplyr::distinct(dplyr::across(dplyr::all_of(zero_keys)))
  zeros$Type <- "FD"; zeros$upper <- 0; zeros$lower <- 0; zeros$SOM <- 0; zeros$bd <- 0; zeros$RVAR <- 0

  d <- dplyr::bind_rows(d, zeros) %>%
    dplyr::arrange(dplyr::across(dplyr::all_of(core_vars)), upper, lower)

  # FD calculations
  d$Soil_g_cm2     <- (d$lower - d$upper) * d$bd
  d$RVAR_g_cm2     <- (d$RVAR / 100) * d$Soil_g_cm2
  d$SOM_g_cm2      <- (d$SOM  / 100) * d$Soil_g_cm2
  d$Min_Soil_g_cm2 <- d$Soil_g_cm2 - d$SOM_g_cm2

  d_fd <- d %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(core_vars))) %>%
    dplyr::mutate(
      Cum_Soil_g_cm2     = cumsum(Soil_g_cm2),
      Cum_RVAR_g_cm2     = cumsum(RVAR_g_cm2),
      Cum_SOM_g_cm2      = cumsum(SOM_g_cm2),
      Cum_Min_Soil_g_cm2 = cumsum(Min_Soil_g_cm2)
    ) %>%
    dplyr::ungroup()

  d_fd <- dplyr::filter(d_fd, !(upper == 0 & lower == 0))
  d_fd$Type <- "FD"

  if (output == "fd") {
    out <- d_fd
    names(out)[names(out) == "RVAR"]       <- rvar
    names(out)[names(out) == "RVAR_g_cm2"] <- paste0(rvar, "_g_cm2")
    return(dplyr::relocate(out, Type, dplyr::all_of(c(core_vars, reference_by, "ref_trt"))))
  }

  # Info: global min with multiple sites (solo informativo)
  if (reference_mode == "min" && is.null(reference_by) && "Site" %in% names(data)) {
    message("Using a single (global) reference across all sites. Set `reference_by = 'Site'` for site-specific minima.")
  }

  # Build reference grid
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
    if (!is.null(reference_by)) for (nm in reference_by) ref_tbl[[nm]] <- NA
    ref_tbl$ref_trt <- "custom"

  } else if (reference_mode == "min") {
    grp_vars <- c(reference_by, "lower")

    min_smm <- d_fd %>%
      dplyr::group_by(dplyr::across(dplyr::all_of(grp_vars))) %>%
      dplyr::summarise(
        Cum_Min_Soil_g_cm2 = min(Cum_Min_Soil_g_cm2, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      dplyr::arrange(dplyr::across(dplyr::all_of(reference_by)), lower)

    ref_tbl <- min_smm %>%
      dplyr::group_by(dplyr::across(dplyr::all_of(reference_by))) %>%
      dplyr::mutate(
        Cum_Min_Soil_g_cm2 = pmax(Cum_Min_Soil_g_cm2, 1e-9),
        Cum_Min_Soil_g_cm2 = cummax(Cum_Min_Soil_g_cm2),
        upper = dplyr::lag(lower, default = 0)
      ) %>%
      dplyr::ungroup() %>%
      dplyr::select(dplyr::all_of(reference_by), upper, lower, Cum_Min_Soil_g_cm2)

    ref_check <- ref_tbl %>%
      dplyr::group_by(dplyr::across(dplyr::all_of(reference_by))) %>%
      dplyr::summarise(any_mass = any(diff(Cum_Min_Soil_g_cm2) > 0, na.rm = TRUE), .groups = "drop")
    if (nrow(ref_check) == 0 || !all(ref_check$any_mass)) {
      stop("Reference grid 'min' has zero mass increments; check FD inputs.", call. = FALSE)
    }

    ref_tbl$ref_trt <- "min"

  } else {  # treatment_mean
    ref_rows <- dplyr::filter(d_fd, trt == ref_trt)
    if (nrow(ref_rows) == 0) stop("No rows for the reference treatment found in `data`.", call. = FALSE)
    ref_depths <- dplyr::distinct(dplyr::select(ref_rows, dplyr::all_of(reference_by), upper, lower))
    ref_tbl <- ref_depths %>%
      dplyr::left_join(
        ref_rows %>%
          dplyr::group_by(dplyr::across(dplyr::all_of(c(reference_by, "upper", "lower")))) %>%
          dplyr::summarise(Cum_Min_Soil_g_cm2 = mean(Cum_Min_Soil_g_cm2, na.rm = TRUE), .groups = "drop"),
        by = c(reference_by, "upper", "lower")
      )
    ref_tbl$ref_trt <- unique(d_fd$ref_trt)
  }

  # ESM interpolation across all cores
  key_cols <- unique(c(core_vars, reference_by, "ref_trt"))
  by_keys  <- dplyr::distinct(dplyr::select(d_fd, dplyr::all_of(key_cols)))

  make_mono_spline <- function(x, y) stats::splinefun(x, y, method = "monoH.FC")

  esm_list <- lapply(seq_len(nrow(by_keys)), function(i) {
    key <- by_keys[i, , drop = FALSE]

    # FD series for this core
    cur <- d_fd
    for (nm in core_vars) {
      if (length(nm) && nm %in% names(cur)) cur <- cur[cur[[nm]] == key[[nm]], , drop = FALSE]
    }
    if (nrow(cur) == 0) return(NULL)

    # Appropriate reference subset (global or grouped)
    ref_subset <- ref_tbl
    if (!is.null(reference_by)) {
      for (nm in reference_by) {
        ref_subset <- ref_subset[is.na(ref_subset[[nm]]) | ref_subset[[nm]] == key[[nm]], , drop = FALSE]
      }
    }
    if (nrow(ref_subset) == 0) return(NULL)

    # Extrapolation rule
    if (identical(extrapolation, "none")) {
      ref_use <- dplyr::filter(ref_subset, Cum_Min_Soil_g_cm2 <= max(cur$Cum_Min_Soil_g_cm2, na.rm = TRUE))
    } else {
      ref_use <- dplyr::filter(ref_subset, lower <= max(cur$lower, na.rm = TRUE))
    }
    if (nrow(ref_use) == 0) return(NULL)

    # Monotone splines in cumulative space
    s_r   <- make_mono_spline(cur$Cum_Min_Soil_g_cm2, cur$Cum_RVAR_g_cm2)
    s_som <- make_mono_spline(cur$Cum_Min_Soil_g_cm2, cur$Cum_SOM_g_cm2)

    # Evaluate on reference mass grid
    ref_use <- ref_use %>%
      dplyr::mutate(
        Cum_Min_Soil_g_cm2 = as.numeric(Cum_Min_Soil_g_cm2),
        Cum_RVAR_g_cm2     = s_r(Cum_Min_Soil_g_cm2),
        Cum_SOM_g_cm2      = s_som(Cum_Min_Soil_g_cm2)
      ) %>%
      dplyr::mutate(
        Cum_RVAR_g_cm2 = cummax(Cum_RVAR_g_cm2),
        Cum_SOM_g_cm2  = cummax(Cum_SOM_g_cm2)
      ) %>%
      dplyr::mutate(
        Min_Soil_g_cm2 = pmax(Cum_Min_Soil_g_cm2 - dplyr::lag(Cum_Min_Soil_g_cm2, default = 0), 0),
        RVAR_g_cm2     = pmax(Cum_RVAR_g_cm2     - dplyr::lag(Cum_RVAR_g_cm2,     default = 0), 0),
        SOM_g_cm2      = pmax(Cum_SOM_g_cm2      - dplyr::lag(Cum_SOM_g_cm2,      default = 0), 0),
        Soil_g_cm2     = pmax(Min_Soil_g_cm2 + SOM_g_cm2, 0),
        bd             = ifelse(lower > upper, Soil_g_cm2 / (lower - upper), NA_real_),
        RVAR           = ifelse(Soil_g_cm2 > 0, RVAR_g_cm2 / Soil_g_cm2 * 100, NA_real_),
        SOM            = ifelse(Soil_g_cm2 > 0, SOM_g_cm2  / Soil_g_cm2 * 100, NA_real_),
        Type           = "ESM",
        Cum_Soil_g_cm2 = cumsum(Soil_g_cm2),
        Cum_RVAR_g_cm2 = cumsum(RVAR_g_cm2),
        Cum_SOM_g_cm2  = cumsum(SOM_g_cm2)
      )

    # Attach identifiers
    for (nm in core_vars)     ref_use[[nm]] <- key[[nm]]
    for (nm in reference_by)  ref_use[[nm]] <- key[[nm]]
    ref_use[["ref_trt"]]      <- key[["ref_trt"]]

    ref_use
  })

  esm_out <- dplyr::bind_rows(esm_list)
  if (nrow(esm_out) == 0) {
    stop("Failed to build ESM series for any core given the reference and extrapolation mode.", call. = FALSE)
  }

  # Column order & names
  base_cols <- c("Type", core_vars, reference_by, "ref_trt",
                 "upper", "lower", "RVAR", "SOM", "bd",
                 "Soil_g_cm2", "SOM_g_cm2", "Min_Soil_g_cm2", "RVAR_g_cm2",
                 "Cum_Soil_g_cm2", "Cum_RVAR_g_cm2", "Cum_SOM_g_cm2", "Cum_Min_Soil_g_cm2")
  base_cols <- intersect(base_cols, names(esm_out))
  esm_out <- esm_out[, base_cols]

  names(esm_out)[names(esm_out) == "RVAR"]           <- rvar
  names(esm_out)[names(esm_out) == "RVAR_g_cm2"]     <- paste0(rvar, "_g_cm2")
  names(esm_out)[names(esm_out) == "Cum_RVAR_g_cm2"] <- paste0("Cum_", rvar, "_g_cm2")

  # Attach reference grid used
  ref_cols <- c(reference_by, "upper", "lower", "Cum_Min_Soil_g_cm2", "ref_trt")
  ref_cols <- intersect(ref_cols, names(ref_tbl))
  attr(esm_out, "SMM_reference") <- ref_tbl[, ref_cols, drop = FALSE]

  dplyr::arrange(esm_out, dplyr::across(dplyr::all_of(c(core_vars, reference_by, "lower"))))
}
