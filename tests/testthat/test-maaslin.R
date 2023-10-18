test_that("mda.maaslin2 works", {
  sink(nullfile()) # I can't seem to shut maaslin up any other way...
  res = mda.maaslin2(D)
  sink()
  expect_named(res)
})
