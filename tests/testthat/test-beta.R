# devtools::check fails on requesting >2 cores
# so we provide a switch here
chk <- Sys.getenv("_R_CHECK_LIMIT_CORES_", "")

if (nzchar(chk) && chk == "TRUE") {
    # use 1 core, we can't have nice things here
    num_workers <- 1L
} else {
    # use all cores in devtools::test()
    num_workers <- parallel::detectCores()
}

test_that("mda.beta works", {
  expect_no_error({
    suppressMessages(mda.beta(D, beta.parallel = num_workers))
    })
})
