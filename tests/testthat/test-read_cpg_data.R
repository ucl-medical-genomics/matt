test_that("read_cpg_data without cache", {
  
  # Load external data
  bed_2f <- system.file("extdata", "bismark_coverage_CpG.bedgraph", package="matt")
  
  # Read CpG data
  dt <- read_cpg_data(
    bed_2f, pipeline="bedgraph",
    align_to_reference = F, use_cache=F,
    collapse_strands = F
  )
  
  # Assert that data.table with expected fields is returned
  expect_true("data.table" %in% class(dt))
  expect_true(nrow(dt) > 0)
  expect_true(
    all(
      c("chr","start","cov","meth") %in% names(dt)
    )
  )
})
