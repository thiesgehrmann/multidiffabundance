test_that("mda.aldex2 works", {
  expect_no_error({
    suppressMessages(mda.aldex2(D))
    })
})
