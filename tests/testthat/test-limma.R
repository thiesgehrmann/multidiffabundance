test_that("mda.limma works", {
  expect_no_error({
    suppressMessages(mda.limma(D))
  })
})