test_that("ref_trt requiere ref_trt en ESM, pero no en FD", {
  skip_on_cran()

  df <- tibble::tribble(
    ~trt, ~rep, ~upper, ~lower, ~bd, ~SOC, ~som, ~ref_trt,
    "cs", 1   , 0     , 5     , 1.1, 2.0 , 4.0, "cs",
    "cs", 1   , 5     , 15    , 1.2, 1.8 , 3.6, "cs",
    "ext",1   , 0     , 5     , 1.0, 2.2 , 4.4, "cs",
    "ext",1   , 5     , 15    , 1.1, 1.9 , 3.8, "cs"
  )

  # FD no necesita ref_trt
  fd <- esm(dplyr::select(df, -ref_trt), rvar = "SOC", som = TRUE, output = "fd")
  expect_s3_class(fd, "tbl_df")

  # ESM ref_trt sí necesita ref_trt
  expect_error(
    esm(dplyr::select(df, -ref_trt), rvar = "SOC", som = TRUE, reference_mode = "ref_trt"),
    "Missing `ref_trt`"
  )

  # Y con ref_trt presente debería funcionar
  trm <- esm(df, rvar = "SOC", som = TRUE, reference_mode = "ref_trt")
  expect_s3_class(trm, "tbl_df")
})
