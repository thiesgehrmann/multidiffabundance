test_that("mda.limma works", {
  expect_no_error({
    suppressMessages(mda.zicoseq(Ds) -> res_zicoseq)
    expect_equal(unlist(res_zicoseq$res.full$comment)[1], "")
  })
})