test_that("mda.beta works", {
  expect_no_error({
    suppressMessages(mda.beta(D, beta.parallel = 1))
    })
})
