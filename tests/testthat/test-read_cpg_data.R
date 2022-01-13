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

test_that("read_cpg_data using cache", {
  # Test: cache file is written when use_cache=TRUE AND dataset_id is not NULL
  
  # Set cache directory to temporary location
  # TODO: Place this in general testthat setup/teardown
  original_cache <- .matt_env[["cache_dir"]]
  .matt_env[["cache_dir"]] <- tempdir()
  test_id = "test_file"
  expected_cache_file = file.path(.matt_env[["cache_dir"]], paste0(test_id, ".fst"))
  
  
  # Load external data
  bed_2f <- system.file("extdata", "bismark_coverage_CpG.bedgraph", package="matt")
  expect_true(!file.exists(expected_cache_file))
  
  # Read CpG data from file
  dt <- read_cpg_data(
    bed_2f, dataset_id=test_id, pipeline="bedgraph",
    align_to_reference = F, use_cache=T,
    collapse_strands = F
  )
  expect_true(file.exists(expected_cache_file))
  
  # Read CpG data from cache
  dt2 <- read_cpg_data(NULL, dataset_id=test_id, use_cache=T)
  expect_true("data.table" %in% class(dt2))
  
  # Cleanup and reset cache
  file.remove(expected_cache_file)
  .matt_env[["cache_dir"]] <- original_cache
})

# Test: correct data is retrieved when use_cache=TRUE AND dataset_id is not NULL
