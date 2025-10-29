test_that("as_Mg_ha() converts correctly", {
  expect_equal(as_Mg_ha(1), 100)
  expect_equal(as_Mg_ha(c(0.5, 2)), c(50, 200))
  expect_error(as_Mg_ha("text"))
})

test_that("esm_aggregate() summarises mean and se correctly", {
  d <- tibble::tibble(
    trt = rep("A", 4),
    lower = rep(15, 4),
    SOC_g_cm2 = c(1, 2, 3, 4)
  )
  agg <- esm_aggregate(d, var = "SOC")
  expect_true("SOC_g_cm2_mean" %in% names(agg))
  expect_equal(agg$SOC_g_cm2_mean, mean(d$SOC_g_cm2))
})
