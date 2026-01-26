test_that("mda.all works", {
  expect_no_error({
    suppressMessages(mda.all(D))
    })
})