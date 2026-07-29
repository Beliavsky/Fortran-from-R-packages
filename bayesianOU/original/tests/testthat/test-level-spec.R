# level_spec resolution / named configurations (no Stan).

test_that("ou_level_spec returns the documented named configurations", {
  full <- list(cubic = TRUE,  sv = TRUE,  student_t = TRUE,  hierarchy = TRUE)
  lean <- list(cubic = FALSE, sv = FALSE, student_t = FALSE, hierarchy = TRUE)
  expect_equal(ou_level_spec("canonical"), list(level1 = full, level2 = lean))
  expect_equal(ou_level_spec("both_full"), list(level1 = full, level2 = full))
  expect_equal(ou_level_spec("both_lean"), list(level1 = lean, level2 = lean))
  expect_equal(ou_level_spec("n1_lean"),   list(level1 = lean, level2 = lean))
  expect_error(ou_level_spec("nonsense"))
})

test_that(".resolve_level_spec validates structure and maps to 8 flags", {
  spec <- bayesianOU:::.resolve_level_spec(NULL, 2L)
  fl <- bayesianOU:::.level_spec_flags(spec)
  expect_equal(unlist(fl[c("l1_cubic","l1_sv","l1_studentt","l1_hier")]),
               c(l1_cubic=1L, l1_sv=1L, l1_studentt=1L, l1_hier=1L))
  expect_equal(unlist(fl[c("l2_cubic","l2_sv","l2_studentt","l2_hier")]),
               c(l2_cubic=0L, l2_sv=0L, l2_studentt=0L, l2_hier=1L))

  fl2 <- bayesianOU:::.level_spec_flags(
    bayesianOU:::.resolve_level_spec(ou_level_spec("both_full"), 2L))
  expect_equal(fl2$l2_cubic, 1L); expect_equal(fl2$l2_sv, 1L)
  expect_equal(fl2$l2_studentt, 1L)

  fl3 <- bayesianOU:::.level_spec_flags(
    bayesianOU:::.resolve_level_spec(ou_level_spec("both_lean"), 2L))
  expect_equal(fl3$l1_cubic, 0L); expect_equal(fl3$l1_sv, 0L)
  expect_equal(fl3$l1_studentt, 0L); expect_equal(fl3$l1_hier, 1L)
})

test_that(".resolve_level_spec rejects malformed specs and lean L1 in single-level", {
  expect_error(bayesianOU:::.resolve_level_spec(list(level1 = list()), 2L),
               "level1")
  bad <- ou_level_spec("both_full"); bad$level1$cubic <- "yes"
  expect_error(bayesianOU:::.resolve_level_spec(bad, 2L), "TRUE/FALSE")
  # Single-level requires a full Level 1.
  expect_error(bayesianOU:::.resolve_level_spec(ou_level_spec("both_lean"), 1L),
               "single-level")
})
