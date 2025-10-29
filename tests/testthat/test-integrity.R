test_that("check/assert profile integrity works", {
  d_ok <- tibble::tibble(
    trt = "A", rep = 1,
    upper = c(0, 5, 15),
    lower = c(5, 15, 30),
    bd = 1.3, SOC = c(2, 1.5, 1)
  )
  # should pass
  diag <- check_profile_integrity(d_ok)
  expect_true(all(diag$ok))
  expect_silent(assert_profile_integrity(d_ok))

  # make a non-contiguous profile (gap)
  d_bad <- tibble::tibble(
    trt = "A", rep = 2,
    upper = c(0, 6, 20),
    lower = c(5, 15, 30),
    bd = 1.3, SOC = c(2, 1.5, 1)
  )
  diag2 <- check_profile_integrity(d_bad)
  expect_false(all(diag2$ok))
  expect_error(assert_profile_integrity(d_bad), "contiguity")
})
