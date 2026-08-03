test_that("setup_parallel works as expected", {
  # Test setup_parallel with default cores (7)
  cl <- setup_parallel(num_cores = 2)  # Use 2 cores to avoid heavy load in tests
  expect_true(inherits(cl, "cluster"))  # Check if the object is a parallel cluster
  expect_equal(length(cl), 2)  # Check if the cluster has 2 cores
  parallel::stopCluster(cl)  # Stop the cluster after the test
})
