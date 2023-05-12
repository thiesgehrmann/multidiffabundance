test_that("mda.lmclr works", {

  res = suppresWarnings(suppressMessages(mda.lmclr(D)))
  expect_snapshot(res)
})