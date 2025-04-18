test_that("mda.aldex2 works", {
  expect_no_error({
    suppressMessages(mda.aldex2(Ds) -> res_aldex)
    })
  expect_equal(unlist(res_aldex$res.full$comment)[1], "")
})
