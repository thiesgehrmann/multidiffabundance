test_that("Can convert tidytacos object and use it downstream", {
  suppressMessages(library(tidytacos))
  expect_no_error({
  fml.urt <- mda.permute_formula(~location+plate)
  D.urt <- mda.from_tidytacos(urt, fml.urt, usecache=FALSE)
  suppressMessages(mda.deseq2(D.urt))
  })
})
