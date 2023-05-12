test_that("mda.aldex2 works", {
  res = suppressMessages(mda.aldex2(D))
  expect_snapshot(res)
})