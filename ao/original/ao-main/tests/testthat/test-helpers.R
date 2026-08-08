test_that("estimate can be split by target", {
  estimate <- 1:5
  target <- c("a", "b", "c")
  npar <- c(2, 2, 1)
  split <- split_by_target(estimate, target, npar)
  checkmate::expect_list(split, types = "integer", len = length(npar))
  expect_equal(split_by_target(estimate, NULL, npar), NULL)
})

test_that("random partition can be generated", {
  min <- 2
  partition <- generate_random_partition(1:5, 0.5, min)
  checkmate::expect_list(partition, types = "integer", min.len = min)
  expect_equal(sort(unlist(partition)), 1:5)
  expect_equal(generate_random_partition(1:5, 0.5, 5), as.list(1:5))
})

test_that("random partition generation is deterministic with a seed", {
  set.seed(123)
  partition_1 <- generate_random_partition(1:10, 0.5, 2)
  set.seed(123)
  partition_2 <- generate_random_partition(1:10, 0.5, 2)
  expect_identical(partition_1, partition_2)
})
