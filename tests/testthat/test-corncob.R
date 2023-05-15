test_that("mda.corncob works", {
  expect_no_error({
    suppressMessages(mda.corncob(D))
  })
})