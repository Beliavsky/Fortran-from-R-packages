
# 1. 로젠브록 함수 정의
rosenbrock <- function(x) {
  100 * (x[2] - x[1]^2)^2 + (1 - x[1])^2
}

# 2. 통합 테스트
test_that("General optimization suite works on Rosenbrock", {
  algorithms <- list(
    bfgs            = bfgs,
    dogleg          = dogleg,
    double_dogleg   = double_dogleg,
    modified_newton = modified_newton,
    newton_raphson  = newton_raphson,
    l_bfgs_b        = l_bfgs_b
  )
  
  start_val <- c(-1.2, 1.0)
  target <- c(1, 1)
  
  for (name in names(algorithms)) {
    res <- algorithms[[name]](
      start     = start_val,
      objective = rosenbrock,
      control   = list(max_iter = 10000, initial_delta = 2.0, tol_rel_x = 1e-6)
    )
    
    expect_true(res$converged, info = paste("Convergence failed:", name))
    expect_equal(res$par, target, tolerance = 1e-2, 
                 info = paste("Parameter mismatch:", name))
  }
})