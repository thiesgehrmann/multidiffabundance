test_that("mda.limma works", {
  expect_no_error({
    expect_warning(
      suppressMessages(mda.limma(D) -> res_limma),
      "One or more quantiles are zero"
    )
  })
  expect_equal(unlist(res_limma$res.full$comment)[1], "")
})