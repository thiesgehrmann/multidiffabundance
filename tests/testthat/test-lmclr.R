test_that("mda.lmclr works", {
  expect_no_error({
  suppressWarnings(suppressMessages(mda.lmclr(D)))
  })
})