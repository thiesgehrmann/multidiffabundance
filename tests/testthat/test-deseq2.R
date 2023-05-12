test_that("mda.deseq2 works", {
  res = suppressMessages(mda.deseq2(D))
  expect_snapshot(res)
})