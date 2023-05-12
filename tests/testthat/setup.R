load(test_path("mda.example.RData"))

D <- mda.create(
  mda.example$count_data,
  mda.example$meta_data,
  mda.example$formulas,
  usecache = FALSE,
  recache = FALSE
)