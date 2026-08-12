library(scip)

## Test with dgCMatrix input
if (requireNamespace("Matrix", quietly = TRUE)) {
    A_dense <- matrix(c(1, 1, 1, 0, 0, 1), nrow = 3, byrow = TRUE)
    A_sparse <- Matrix::Matrix(A_dense, sparse = TRUE)
    b <- c(10, 6, 8)

    res_dense <- scip_solve(obj = c(-1, -1), A = A_dense, b = b,
                            sense = rep("<=", 3),
                            control = list(verbose = FALSE))
    res_sparse <- scip_solve(obj = c(-1, -1), A = A_sparse, b = b,
                             sense = rep("<=", 3),
                             control = list(verbose = FALSE))

    expect_equal(res_dense$status, res_sparse$status)
    expect_equal(res_dense$objval, res_sparse$objval)
    expect_equal(res_dense$x, res_sparse$x)
}

## Test with simple_triplet_matrix input
if (requireNamespace("slam", quietly = TRUE)) {
    A_stm <- slam::simple_triplet_matrix(
        i = c(1, 1, 2, 3), j = c(1, 2, 1, 2),
        v = c(1, 1, 1, 1), nrow = 3, ncol = 2
    )
    b <- c(10, 6, 8)
    res_stm <- scip_solve(obj = c(-1, -1), A = A_stm, b = b,
                          sense = rep("<=", 3),
                          control = list(verbose = FALSE))

    expect_equal(res_stm$status, "optimal")
    expect_equal(res_stm$objval, -10)
}
