test_that("mda.alpha works", {
  expect_no_error({
    suppressMessages(mda.alpha(D))
    })
})

test_that("mda.alpha throws error when non-existing index given", {
  expect_error({
    suppressMessages(mda.alpha(D, alpha.index = "nope"),
    "[MDA] mda.alpha: No valid index value(s) defined. Options must be [shannon, simpson, invsimpson]."
    )
    })
})