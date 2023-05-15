test_that("mda.deseq2 works", {
  expect_no_error({
    suppressMessages(mda.deseq2(D))
  })
})