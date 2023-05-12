test_that("mda.corncob works", {
  res = suppressMessages(mda.corncob(D))
  expect_snapshot(res)
})