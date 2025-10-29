test_that("min_smm_by_depth reproduce FD→min(lower)→custom", {
  df_test <- tibble::tibble(
    trt   = rep(c("A","B"), each = 6),
    rep   = rep(rep(1:2, each = 3), times = 2),  # 0–5–15 para rep=1 y rep=2
    upper = rep(c(0,5,15), times = 4),
    lower = rep(c(5,15,30), times = 4),
    bd    = rep(c(1.2,1.3,1.35, 1.1,1.25,1.3), times = 2),
    SOC   = c(2.0,1.5,1.0, 2.1,1.4,0.9, 1.8,1.2,0.8, 1.9,1.1,0.7)
  )

  fd <- esm(df_test, som = FALSE, soc_col = "SOC", rvar = "SOC", output = "fd")

  d.ref_min <- fd %>%
    dplyr::group_by(lower) %>%
    dplyr::reframe(Cum_Min_Soil_g_cm2 = min(Cum_Min_Soil_g_cm2))

  em.cust <- esm(df_test, reference_mode = "custom",
                 rvar = "SOC", som = FALSE, soc_col = "SOC",
                 custom_ld = d.ref_min$lower,
                 custom_smm = d.ref_min$Cum_Min_Soil_g_cm2)

  esmn2 <- esm(df_test, som = FALSE, soc_col = "SOC",
               rvar = "SOC", reference_mode = "min_smm_by_depth")

  ref_custom <- attr(em.cust, "SMM_reference") %>%
    dplyr::select(upper, lower, Cum_Min_Soil_g_cm2)
  ref_min <- attr(esmn2, "SMM_reference") %>%
    dplyr::select(upper, lower, Cum_Min_Soil_g_cm2)

  testthat::expect_equal(ref_custom$Cum_Min_Soil_g_cm2,
                         ref_min$Cum_Min_Soil_g_cm2, tolerance = 1e-6)

  ok_by_profile <- esmn2 %>%
    dplyr::arrange(trt, rep, lower) %>%
    dplyr::group_by(trt, rep) %>%
    dplyr::summarise(ok = all(diff(Cum_SOC_g_cm2) >= -1e-6), .groups = "drop") %>%
    dplyr::pull(ok)

  testthat::expect_true(all(ok_by_profile))
})
