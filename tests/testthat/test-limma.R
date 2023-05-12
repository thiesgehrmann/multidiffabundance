test_that("mda.limma works", {
  res = suppressMessages(mda.limma(D))
  expect_snapshot(res)
})