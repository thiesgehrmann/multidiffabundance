test_that("mda.ancombc works", {
  expect_no_error({
    suppressMessages(mda.ancombc(D))
  })
})

test_that("mda.ancombc2 works", {
  expect_no_error({
    suppressMessages(mda.ancombc2(D))
  })
})

