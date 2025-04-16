test_that("mda.ancombc2 works", {
  expect_no_error({
    suppressMessages(mda.ancombc2(D) -> res_ancomb)
  })


  # expect the function to output atleast one result
  #TODO: uncomment when you feel like fixing this expect_true(any(!is.na(res_ancomb$res.full$stat)))
})
