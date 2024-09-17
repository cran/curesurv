context("Cure regression")

test_that("curesurv correctly checks inits values", {
  # Test the presence of init values
  expect_error(
    curesurv::curesurv(Surv(time = time_obs, event = event) ~ age_cr,
                         data = pancreas_data,
                         riskpop =  "ehazard",
                         riskpop.alpha = TRUE,
                       method_opt = "L-BFGS-B",
                       init = c(0,0,0,0,0))
    )
})
