test_that("mda.ancombc works", {
  res = suppressMessages(mda.ancombc(D))
  expect_snapshot(res)
})

test_that("mda.ancombc2 works", {
  res = suppressMessages(mda.ancombc2(D))
  expect_snapshot(res)
})

