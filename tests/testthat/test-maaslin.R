test_that("mda.maaslin2 works", {
  sink(nullfile()) # I can't seem to shut maaslin up any other way...
  res = mda.maaslin2(D)
  sink()
  expect_named(res)
})

# edgecase where maaslin2 fails when only testing one variable
test_that("mda.maaslin2 works for single covariate", {
  expect_no_error({
    res = mda.maaslin2(Ds)
  })
})