test_that("mda.deseq2 works", {
  expect_no_error({
    suppressMessages(mda.deseq2(Ds) -> res_deseq2)
  })
  # no comments expected
  expect_equal(unlist(res_deseq2$res.full$comment)[1], "")
})