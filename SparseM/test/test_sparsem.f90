program test_sparsem
    use sparsem, only: dp, csr_matrix, csc_matrix, coo_matrix, ssr_matrix, ssc_matrix
    use sparsem, only: csr_cholesky, csr_cholesky_stats, slm_result, slm_summary_result
    use sparsem, only: hb_matrix, mslm_result, mslm_summary_result
    use sparsem, only: csr_from_dense, csr_to_dense, csc_from_csr, csr_from_csc
    use sparsem, only: coo_from_csr, csr_from_coo, ssr_from_csr, csr_from_ssr
    use sparsem, only: ssc_from_csr, csr_from_ssc, csr_transpose, csr_matmul
    use sparsem, only: csr_diag, csr_set_diag, csr_norm, csr_rbind, csr_cbind
    use sparsem, only: csr_submatrix, csr_assign_submatrix, csr_kronecker
    use sparsem, only: csr_chol, chol_statistics, chol_as_csr, chol_determinant
    use sparsem, only: chol_forwardsolve, chol_backsolve, chol_solve
    use sparsem, only: chol_forwardsolve_matrix, chol_backsolve_matrix, chol_solve_matrix
    use sparsem, only: chol_forwardsolve_csr, chol_backsolve_csr, chol_solve_csr
    use sparsem, only: csr_solve, csr_solve_matrix, csr_solve_csr_rhs, csr_inverse, csr_determinant, slm_fit_csr
    use sparsem, only: slm_deviance, slm_extract_aic, slm_summary
    use sparsem, only: read_matrix_hb, hb_model_matrix, hb_model_response
    use sparsem, only: hb_model_guess, hb_model_xexact, mslm_fit_csr, mslm_summary
    use sparsem, only: minimum_degree_order
    implicit none

    real(dp), parameter :: tol = 1.0e-11_dp
    real(dp) :: a_dense(3, 3)
    real(dp) :: b_dense(3, 3)
    real(dp) :: x_dense(4, 2)
    real(dp) :: y(4)
    real(dp) :: ym(4, 2)
    real(dp) :: tri_dense(8, 8)
    real(dp) :: tri_rhs(8)
    real(dp) :: fill_dense(4, 4)
    real(dp) :: fill_rhs(4)
    real(dp) :: star_dense(10, 10)
    real(dp) :: star_rhs(10)
    real(dp) :: star_rhs_matrix(10, 2)
    real(dp), allocatable :: dense(:,:)
    real(dp), allocatable :: sol(:)
    real(dp), allocatable :: forward_sol(:)
    real(dp), allocatable :: composed_sol(:)
    real(dp), allocatable :: solm(:,:)
    real(dp), allocatable :: diagv(:)
    real(dp), allocatable :: hb_rhs(:,:)
    real(dp), allocatable :: hb_guess(:,:)
    real(dp), allocatable :: hb_xexact(:,:)
    real(dp) :: det
    real(dp) :: aic(2)
    type(csr_matrix) :: a
    type(csr_matrix) :: b
    type(csr_matrix) :: c
    type(csr_matrix) :: sub
    type(csr_matrix) :: replacement
    type(csc_matrix) :: csc
    type(coo_matrix) :: coo
    type(ssr_matrix) :: ssr
    type(ssc_matrix) :: ssc
    type(csr_cholesky) :: factor
    type(csr_cholesky) :: level_factor
    type(csr_cholesky) :: natural_factor
    type(csr_cholesky_stats) :: factor_stats
    type(slm_result) :: fit
    type(slm_summary_result) :: smry
    type(hb_matrix) :: hb
    type(mslm_result) :: mfit
    type(mslm_summary_result) :: msmry
    character(len=256) :: hb_path
    logical :: exists
    integer :: info
    integer :: det_sign
    integer :: i
    integer :: lev
    integer :: nofsub
    integer, allocatable :: mmd_perm(:)
    integer, allocatable :: mmd_inv_perm(:)

    a_dense = 0.0_dp
    a_dense(1, :) = [4.0_dp, 1.0_dp, 0.0_dp]
    a_dense(2, :) = [1.0_dp, 3.0_dp, 1.0_dp]
    a_dense(3, :) = [0.0_dp, 1.0_dp, 2.0_dp]
    a = csr_from_dense(a_dense)
    call assert_true(a%nrow == 3 .and. a%ncol == 3, 'csr dimensions')
    call assert_true(size(a%values) == 7, 'csr nnz')
    dense = csr_to_dense(a)
    call assert_matrix_close(dense, a_dense, tol, 'dense round trip')

    csc = csc_from_csr(a)
    c = csr_from_csc(csc)
    call assert_matrix_close(csr_to_dense(c), a_dense, tol, 'csc round trip')
    coo = coo_from_csr(a)
    c = csr_from_coo(coo)
    call assert_matrix_close(csr_to_dense(c), a_dense, tol, 'coo round trip')
    ssr = ssr_from_csr(a, info=info)
    call assert_true(info == 0, 'ssr conversion status')
    c = csr_from_ssr(ssr)
    call assert_matrix_close(csr_to_dense(c), a_dense, tol, 'ssr round trip')
    ssc = ssc_from_csr(a, info=info)
    call assert_true(info == 0, 'ssc conversion status')
    c = csr_from_ssc(ssc)
    call assert_matrix_close(csr_to_dense(c), a_dense, tol, 'ssc round trip')

    c = csr_transpose(a)
    call assert_matrix_close(csr_to_dense(c), transpose(a_dense), tol, 'transpose')
    c = csr_matmul(a, a, info)
    call assert_true(info == 0, 'matmul status')
    call assert_matrix_close(csr_to_dense(c), matmul(a_dense, a_dense), tol, 'matmul')

    diagv = csr_diag(a)
    call assert_vector_close(diagv, [4.0_dp, 3.0_dp, 2.0_dp], tol, 'diag')
    c = csr_set_diag(a, [5.0_dp], info)
    call assert_true(info == 0, 'set diag status')
    call assert_vector_close(csr_diag(c), [5.0_dp, 5.0_dp, 5.0_dp], tol, 'set diag')
    call assert_close(csr_norm(a, 'F', info), sqrt(33.0_dp), tol, 'F norm')
    call assert_true(info == 0, 'norm status')

    b_dense = 0.0_dp
    b_dense(1, 1) = 2.0_dp
    b_dense(2, 2) = 3.0_dp
    b_dense(3, 3) = 4.0_dp
    b = csr_from_dense(b_dense)
    c = csr_rbind(a, b, info)
    call assert_true(info == 0 .and. c%nrow == 6 .and. c%ncol == 3, 'rbind')
    c = csr_cbind(a, b, info)
    call assert_true(info == 0 .and. c%nrow == 3 .and. c%ncol == 6, 'cbind')

    sub = csr_submatrix(a, [1, 3], [1, 2], info)
    call assert_true(info == 0, 'submatrix status')
    call assert_matrix_close(csr_to_dense(sub), reshape([4.0_dp, 0.0_dp, 1.0_dp, 1.0_dp], [2, 2]), &
        tol, 'submatrix')
    replacement = csr_from_dense(reshape([9.0_dp], [1, 1]))
    c = csr_assign_submatrix(a, [1], [3], replacement, info)
    call assert_true(info == 0, 'submatrix assignment status')
    dense = csr_to_dense(c)
    call assert_close(dense(1, 3), 9.0_dp, tol, 'submatrix assignment')

    c = csr_kronecker(csr_from_dense(reshape([1.0_dp, 0.0_dp, 0.0_dp, 2.0_dp], [2, 2])), &
        csr_from_dense(reshape([3.0_dp, 4.0_dp, 5.0_dp, 6.0_dp], [2, 2])))
    call assert_true(c%nrow == 4 .and. c%ncol == 4, 'kronecker dimensions')

    factor = csr_chol(a, info)
    call assert_true(info == 0, 'cholesky status')
    sol = chol_solve(factor, [1.0_dp, 2.0_dp, 3.0_dp], info)
    call assert_true(info == 0, 'solve status')
    call assert_vector_close(sol, [2.0_dp / 9.0_dp, 1.0_dp / 9.0_dp, 13.0_dp / 9.0_dp], tol, 'solve')
    solm = csr_solve_matrix(a, reshape([1.0_dp, 2.0_dp, 3.0_dp, 3.0_dp, 2.0_dp, 1.0_dp], [3, 2]), info)
    call assert_true(info == 0 .and. all(shape(solm) == [3, 2]), 'matrix solve status')
    call assert_vector_close(matmul(a_dense, solm(:, 1)), [1.0_dp, 2.0_dp, 3.0_dp], 1.0e-10_dp, 'matrix solve rhs 1')
    call assert_vector_close(matmul(a_dense, solm(:, 2)), [3.0_dp, 2.0_dp, 1.0_dp], 1.0e-10_dp, 'matrix solve rhs 2')
    replacement = csr_from_dense(reshape([1.0_dp, 2.0_dp, 3.0_dp, 3.0_dp, 2.0_dp, 1.0_dp], [3, 2]))
    solm = csr_solve_csr_rhs(a, replacement, info)
    call assert_true(info == 0 .and. all(shape(solm) == [3, 2]), 'sparse-RHS matrix solve status')
    call assert_vector_close(matmul(a_dense, solm(:, 1)), [1.0_dp, 2.0_dp, 3.0_dp], 1.0e-10_dp, &
        'sparse-RHS matrix solve first column')
    solm = chol_solve_csr(factor, replacement, info)
    call assert_true(info == 0, 'factor sparse-RHS full solve status')
    call assert_vector_close(matmul(a_dense, solm(:, 2)), [3.0_dp, 2.0_dp, 1.0_dp], 1.0e-10_dp, &
        'factor sparse-RHS full solve second column')
    solm = chol_forwardsolve_csr(factor, replacement, info)
    call assert_true(info == 0, 'factor sparse-RHS forwardsolve status')
    replacement = csr_from_dense(solm)
    solm = chol_backsolve_csr(factor, replacement, info)
    call assert_true(info == 0, 'factor sparse-RHS backsolve status')
    call assert_vector_close(matmul(a_dense, solm(:, 1)), [1.0_dp, 2.0_dp, 3.0_dp], 1.0e-10_dp, &
        'factor sparse-RHS triangular composition')
    det = csr_determinant(a, logarithm=.false., info=info)
    call assert_true(info == 0, 'determinant status')
    call assert_close(det, 18.0_dp, 1.0e-10_dp, 'determinant')
    det = chol_determinant(factor, logarithm=.false., sign=det_sign, info=info)
    call assert_true(info == 0 .and. det_sign == 1, 'Cholesky-factor determinant status')
    call assert_close(det, sqrt(18.0_dp), 1.0e-10_dp, 'Cholesky-factor determinant')
    call assert_close(factor%factor_determinant, sqrt(18.0_dp), 1.0e-10_dp, &
        'cached Cholesky-factor determinant')
    det = chol_determinant(factor, sign=det_sign, info=info)
    call assert_close(det, 0.5_dp * log(18.0_dp), 1.0e-10_dp, 'Cholesky-factor log determinant')
    call assert_close(factor%factor_log_determinant, 0.5_dp * log(18.0_dp), 1.0e-10_dp, &
        'cached Cholesky-factor log determinant')
    c = chol_as_csr(factor, info=info)
    call assert_true(info == 0, 'Cholesky-factor upper CSR conversion status')
    call assert_matrix_close(csr_to_dense(c), csr_to_dense(factor%upper), tol, 'Cholesky-factor upper CSR conversion')
    c = chol_as_csr(factor, upper_tri=.false., info=info)
    call assert_true(info == 0, 'Cholesky-factor lower CSR conversion status')
    call assert_matrix_close(csr_to_dense(c), csr_to_dense(factor%lower), tol, 'Cholesky-factor lower CSR conversion')
    c = csr_inverse(a, info)
    call assert_true(info == 0, 'CSR inverse status')
    dense = matmul(a_dense, csr_to_dense(c))
    call assert_matrix_close(dense, reshape([1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, &
        0.0_dp, 0.0_dp, 1.0_dp], [3, 3]), 1.0e-10_dp, 'CSR inverse')

    tri_dense = 0.0_dp
    do i = 1, 8
        tri_dense(i, i) = 4.0_dp
    end do
    do i = 1, 7
        tri_dense(i + 1, i) = -1.0_dp
        tri_dense(i, i + 1) = -1.0_dp
    end do
    factor = csr_chol(csr_from_dense(tri_dense), info)
    call assert_true(info == 0, 'sparse tridiagonal cholesky status')
    call assert_true(size(factor%lower%values) == 15, 'sparse tridiagonal cholesky nnz')
    call assert_true(size(factor%upper%values) == 15, 'sparse transpose cholesky nnz')
    tri_rhs = [(real(i, dp), i = 1, 8)]
    sol = chol_solve(factor, tri_rhs, info)
    call assert_true(info == 0, 'sparse tridiagonal solve status')
    call assert_vector_close(matmul(tri_dense, sol), tri_rhs, 1.0e-10_dp, 'sparse tridiagonal solve')

    fill_dense = reshape([ &
        4.0_dp, 1.0_dp, 0.0_dp, 1.0_dp, &
        1.0_dp, 4.0_dp, 1.0_dp, 0.0_dp, &
        0.0_dp, 1.0_dp, 4.0_dp, 1.0_dp, &
        1.0_dp, 0.0_dp, 1.0_dp, 4.0_dp], [4, 4])
    factor = csr_chol(csr_from_dense(fill_dense), info)
    call assert_true(info == 0, 'sparse fill-in cholesky status')
    call assert_true(size(factor%lower%values) == 9, 'sparse cholesky structural fill')
    fill_rhs = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp]
    sol = chol_solve(factor, fill_rhs, info)
    call assert_true(info == 0, 'sparse fill-in solve status')
    call assert_vector_close(matmul(fill_dense, sol), fill_rhs, 1.0e-10_dp, 'sparse fill-in solve')

    star_dense = 0.0_dp
    star_dense(1, 1) = 10.0_dp
    do i = 2, 10
        star_dense(i, i) = 2.0_dp
        star_dense(1, i) = 1.0_dp
        star_dense(i, 1) = 1.0_dp
    end do
    a = csr_from_dense(star_dense)
    call minimum_degree_order(a, mmd_perm, mmd_inv_perm, nofsub, info)
    call assert_true(info == 0, 'minimum-degree ordering status')
    call assert_true(all(mmd_perm == [10, 9, 8, 7, 6, 5, 4, 3, 2, 1]), 'minimum-degree star ordering')
    call assert_true(all(mmd_inv_perm(mmd_perm) == [(i, i = 1, 10)]), 'minimum-degree inverse permutation')
    call assert_true(nofsub == 9, 'minimum-degree subscript estimate')
    natural_factor = csr_chol(a, info, ordering='natural')
    call assert_true(info == 0, 'natural-order star cholesky status')
    call assert_true(trim(natural_factor%numeric_kernel) == 'row_sparse', 'natural numerical kernel recorded')
    call assert_true(size(natural_factor%lower%values) == 55, 'natural-order star fill')
    factor = csr_chol(a, info)
    call assert_true(info == 0, 'minimum-degree star cholesky status')
    call assert_true(trim(factor%ordering) == 'minimum_degree', 'minimum-degree ordering recorded')
    call assert_true(trim(factor%numeric_kernel) == 'ng_peyton_level8', 'Ng-Peyton numerical kernel recorded')
    call assert_true(size(factor%lower%values) == 19, 'minimum-degree star fill reduction')
    call assert_true(size(factor%lower%values) < size(natural_factor%lower%values), 'minimum-degree reduces fill')
    call assert_true(factor%nnzl == 19, 'symbolic star factor nonzero count')
    call assert_true(factor%nsub == 18, 'symbolic star compressed subscript count')
    call assert_true(factor%nsuper == 9, 'symbolic star supernode count')
    call assert_true(sum(factor%col_count) == factor%nnzl, 'symbolic column counts total')
    call assert_true(factor%xlnz(11) - 1 == factor%nnzl, 'supernodal numerical pointer extent')
    call assert_true(size(factor%lnz) == factor%nnzl, 'supernodal numerical value extent')
    call assert_true(all(factor%xsuper(1:factor%nsuper + 1) == [1, 2, 3, 4, 5, 6, 7, 8, 9, 11]), &
        'symbolic star supernode partition')
    call assert_true(all(factor%col_count == [2, 2, 2, 2, 2, 2, 2, 2, 2, 1]), &
        'symbolic star column counts')
    call assert_true(factor%level == 8, 'Ng-Peyton default level metadata')
    call assert_true(factor%cache_kb == 1024, 'Ng-Peyton default cache metadata')
    factor_stats = chol_statistics(factor)
    call assert_true(factor_stats%supernodal, 'Ng-Peyton factor statistics available')
    call assert_true(factor_stats%nsuper == 9 .and. factor_stats%nnzl == 19 .and. factor_stats%nsub == 18, &
        'Ng-Peyton factor statistics counts')
    call assert_true(factor_stats%max_supernode_columns == 2, 'Ng-Peyton largest supernode column count')
    call assert_true(factor_stats%max_supernode_nonzeros == 3, 'Ng-Peyton largest supernode nonzero count')
    call assert_close(factor_stats%factorization_flops, 27.0_dp, 0.0_dp, 'Ng-Peyton factor operation count')
    call assert_close(factor_stats%triangular_solve_flops, 56.0_dp, 0.0_dp, 'Ng-Peyton solve operation count')
    call assert_true(factor_stats%tmp_size == factor%tmp_size, 'Ng-Peyton factor workspace statistic')
    call assert_close(factor%tiny, 1.0e-30_dp, 0.0_dp, 'Ng-Peyton default tiny metadata')
    call assert_close(factor%large, 1.0e128_dp, 0.0_dp, 'Ng-Peyton default large metadata')
    star_rhs = [(real(i, dp), i = 1, 10)]
    sol = chol_solve(factor, star_rhs, info)
    call assert_true(info == 0, 'minimum-degree star solve status')
    call assert_vector_close(matmul(star_dense, sol), star_rhs, 1.0e-10_dp, 'minimum-degree star solve')
    forward_sol = chol_forwardsolve(factor, star_rhs, info)
    call assert_true(info == 0, 'Ng-Peyton block forwardsolve status')
    composed_sol = chol_backsolve(factor, forward_sol, info)
    call assert_true(info == 0, 'Ng-Peyton block backsolve status')
    call assert_vector_close(composed_sol, sol, 1.0e-11_dp, 'Ng-Peyton block triangular solve composition')

    do lev = 1, 8
        if (all(lev /= [1, 2, 4, 8])) cycle
        level_factor = csr_chol(a, info, level=lev, cache_kb=64)
        call assert_true(info == 0, 'Ng-Peyton selectable numerical level status')
        call assert_true(level_factor%level == lev, 'Ng-Peyton selectable numerical level metadata')
        call assert_true(level_factor%cache_kb == 64, 'Ng-Peyton selectable cache metadata')
        sol = chol_solve(level_factor, star_rhs, info)
        call assert_true(info == 0, 'Ng-Peyton selectable numerical level solve status')
        call assert_vector_close(matmul(star_dense, sol), star_rhs, 1.0e-10_dp, &
            'Ng-Peyton selectable numerical level solve')
    end do

    level_factor = csr_chol(a, info, level=3)
    call assert_true(info == -4 .and. level_factor%info == -4, 'invalid Ng-Peyton level rejected')
    level_factor = csr_chol(a, info, cache_kb=16)
    call assert_true(info == -5 .and. level_factor%info == -5, 'invalid Ng-Peyton cache size rejected')
    level_factor = csr_chol(a, info, tiny=-1.0_dp)
    call assert_true(info == -6 .and. level_factor%info == -6, 'invalid Ng-Peyton tiny threshold rejected')
    level_factor = csr_chol(a, info, nnzlmax=18)
    call assert_true(info == 4 .and. level_factor%info == 4, 'Ng-Peyton nnzlmax ceiling status')
    level_factor = csr_chol(a, info, nsubmax=17)
    call assert_true(info == 5 .and. level_factor%info == 5, 'Ng-Peyton nsubmax ceiling status')
    level_factor = csr_chol(a, info, tmpmax=0)
    call assert_true(info == 8 .and. level_factor%info == 8, 'Ng-Peyton tmpmax ceiling status')

    b_dense = 0.0_dp
    b_dense(1, 1) = 1.0_dp
    b_dense(2, 2) = -1.0_dp
    level_factor = csr_chol(csr_from_dense(b_dense(1:2, 1:2)), info)
    call assert_true(info == 17 .and. level_factor%info == 0, 'Ng-Peyton indefinite diagonal tiny replacement')

    star_rhs_matrix(:, 1) = star_rhs
    star_rhs_matrix(:, 2) = 2.0_dp * star_rhs
    solm = chol_solve_matrix(factor, star_rhs_matrix, info)
    call assert_true(info == 0, 'Ng-Peyton direct multi-RHS solve status')
    call assert_vector_close(matmul(star_dense, solm(:, 1)), star_rhs, 1.0e-10_dp, &
        'Ng-Peyton direct multi-RHS solve first column')
    call assert_vector_close(matmul(star_dense, solm(:, 2)), 2.0_dp * star_rhs, 1.0e-10_dp, &
        'Ng-Peyton direct multi-RHS solve second column')
    solm = chol_forwardsolve_matrix(factor, star_rhs_matrix, info)
    call assert_true(info == 0, 'Ng-Peyton direct multi-RHS forwardsolve status')
    solm = chol_backsolve_matrix(factor, solm, info)
    call assert_true(info == 0, 'Ng-Peyton direct multi-RHS backsolve status')
    call assert_vector_close(matmul(star_dense, solm(:, 1)), star_rhs, 1.0e-10_dp, &
        'Ng-Peyton direct multi-RHS triangular composition')

    b_dense = 0.0_dp
    b_dense(1, 1) = 2.0_dp
    b_dense(1, 2) = 1.0_dp
    b_dense(2, 2) = 2.0_dp
    level_factor = csr_chol(csr_from_dense(b_dense(1:2, 1:2)), info)
    call assert_true(info == -2 .and. level_factor%info == -2, 'nonsymmetric Cholesky input rejected')

    b_dense = 0.0_dp
    b_dense(1, 1) = 1.0e-40_dp
    b_dense(2, 2) = 1.0_dp
    level_factor = csr_chol(csr_from_dense(b_dense(1:2, 1:2)), info)
    call assert_true(info == 17, 'tiny Cholesky pivot warning status')
    call assert_true(level_factor%info == 0, 'tiny Cholesky pivot factor remains usable')
    call assert_true(level_factor%numeric_info == 17 .and. level_factor%tiny_pivots == 1, &
        'tiny Cholesky pivot metadata')
    sol = chol_solve(level_factor, [1.0_dp, 1.0_dp], info)
    call assert_true(info == 0, 'tiny-pivot Cholesky factor solve status')
    call assert_close(sol(1), 1.0e-128_dp, 1.0e-140_dp, 'tiny-pivot replacement solution')
    call assert_close(sol(2), 1.0_dp, 1.0e-12_dp, 'tiny-pivot unaffected solution')
    sol = csr_solve(csr_from_dense(b_dense(1:2, 1:2)), [1.0_dp, 1.0_dp], info)
    call assert_true(info == 17 .and. size(sol) == 2, 'tiny-pivot warning propagated through csr_solve')
    det = csr_determinant(csr_from_dense(b_dense(1:2, 1:2)), logarithm=.false., info=info)
    call assert_true(info == 17, 'tiny-pivot warning propagated through determinant')
    call assert_close(det, 1.0e128_dp, 1.0e116_dp, 'tiny-pivot replacement determinant')

    x_dense(:, 1) = 1.0_dp
    x_dense(:, 2) = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp]
    y = [1.0_dp, 2.0_dp, 2.0_dp, 4.0_dp]
    fit = slm_fit_csr(csr_from_dense(x_dense), y, info)
    call assert_true(info == 0, 'slm fit status')
    call assert_vector_close(fit%coefficients, [0.0_dp, 0.9_dp], 1.0e-10_dp, 'slm coefficients')
    call assert_close(slm_deviance(fit), 0.7_dp, 1.0e-10_dp, 'slm deviance')
    aic = slm_extract_aic(fit)
    call assert_close(aic(1), 2.0_dp, tol, 'slm aic edf')
    smry = slm_summary(fit, y, intercept=.true., info=info)
    call assert_true(info == 0, 'slm summary status')
    call assert_close(smry%sigma, sqrt(0.35_dp), 1.0e-10_dp, 'slm sigma')
    call assert_close(smry%covariance(1, 1), 0.525_dp, 1.0e-10_dp, 'slm covariance')
    call assert_close(smry%p_value(2), 0.0766194831233613_dp, 1.0e-12_dp, 'slm coefficient p value')
    call assert_true(smry%has_f_statistic, 'slm F statistic present')
    call assert_close(smry%f_statistic, 11.5714285714286_dp, 1.0e-11_dp, 'slm F statistic')
    call assert_close(smry%f_p_value, 0.0766194831233613_dp, 1.0e-12_dp, 'slm F p value')
    call assert_close(smry%correlation(1, 1), 1.0_dp, 1.0e-12_dp, 'slm coefficient correlation diagonal')

    ym(:, 1) = y
    ym(:, 2) = [2.0_dp, 1.0_dp, 3.0_dp, 5.0_dp]
    mfit = mslm_fit_csr(csr_from_dense(x_dense), ym, info)
    call assert_true(info == 0, 'mslm fit status')
    call assert_true(all(shape(mfit%coefficients) == [2, 2]), 'mslm coefficient shape')
    call assert_vector_close(mfit%coefficients(:, 1), [0.0_dp, 0.9_dp], 1.0e-10_dp, 'mslm response 1 coefficients')
    call assert_vector_close(mfit%coefficients(:, 2), [0.0_dp, 1.1_dp], 1.0e-10_dp, 'mslm response 2 coefficients')
    msmry = mslm_summary(mfit, ym, intercept=.true., info=info)
    call assert_true(info == 0, 'mslm summary status')
    call assert_true(size(msmry%responses) == 2, 'mslm summary response count')
    call assert_close(msmry%responses(1)%rss, 0.7_dp, 1.0e-10_dp, 'mslm response 1 rss')
    call assert_close(msmry%responses(2)%rss, 2.7_dp, 1.0e-10_dp, 'mslm response 2 rss')
    call assert_close(msmry%responses(2)%f_p_value, 0.168478159379700_dp, 1.0e-12_dp, 'mslm response 2 F p value')

    call test_ng_numeric_levels()

    hb_path = 'test/data/rua_32_ax.rua'
    inquire(file=trim(hb_path), exist=exists)
    if (.not. exists) then
        hb_path = '../test/data/rua_32_ax.rua'
        inquire(file=trim(hb_path), exist=exists)
    end if
    call assert_true(exists, 'Harwell-Boeing test fixture exists')
    hb = read_matrix_hb(trim(hb_path), info)
    call assert_true(info == 0, 'Harwell-Boeing read status')
    call assert_true(hb%nrow == 32 .and. hb%ncol == 32, 'Harwell-Boeing dimensions')
    call assert_true(size(hb%values) == 126, 'Harwell-Boeing nnz')
    call assert_close(hb%values(1), 101.0_dp, tol, 'Harwell-Boeing first value')
    call assert_close(hb%values(126), 3232.0_dp, tol, 'Harwell-Boeing last value')
    c = hb_model_matrix(hb, info)
    call assert_true(info == 0, 'Harwell-Boeing model matrix status')
    dense = csr_to_dense(c)
    call assert_close(dense(1, 1), 101.0_dp, tol, 'Harwell-Boeing matrix entry 1,1')
    call assert_close(dense(26, 1), 126.0_dp, tol, 'Harwell-Boeing matrix entry 26,1')
    call assert_close(dense(32, 32), 3232.0_dp, tol, 'Harwell-Boeing matrix entry 32,32')
    hb_rhs = hb_model_response(hb, info)
    call assert_true(info == 0 .and. all(shape(hb_rhs) == [32, 2]), 'Harwell-Boeing rhs shape')
    call assert_close(hb_rhs(1, 1), 1001.0_dp, tol, 'Harwell-Boeing rhs first')
    call assert_close(hb_rhs(1, 2), 4606.0_dp, tol, 'Harwell-Boeing rhs second column')
    hb_guess = hb_model_guess(hb, info)
    call assert_true(info == 0 .and. all(shape(hb_guess) == [32, 2]), 'Harwell-Boeing guess shape')
    call assert_close(hb_guess(32, 2), 1.0_dp, tol, 'Harwell-Boeing guess value')
    hb_xexact = hb_model_xexact(hb, info)
    call assert_true(info == 0 .and. all(shape(hb_xexact) == [32, 2]), 'Harwell-Boeing exact shape')
    call assert_close(hb_xexact(10, 1), 1.0_dp, tol, 'Harwell-Boeing exact value')
    call assert_close(hb_xexact(1, 2), 1.0_dp, tol, 'Harwell-Boeing exact second response')

    hb_path = 'test/data/symmetric_3.rsa'
    inquire(file=trim(hb_path), exist=exists)
    if (.not. exists) then
        hb_path = '../test/data/symmetric_3.rsa'
        inquire(file=trim(hb_path), exist=exists)
    end if
    call assert_true(exists, 'symmetric Harwell-Boeing test fixture exists')
    hb = read_matrix_hb(trim(hb_path), info)
    call assert_true(info == 0 .and. hb%symmetric, 'symmetric Harwell-Boeing read')
    c = hb_model_matrix(hb, info)
    call assert_true(info == 0, 'symmetric Harwell-Boeing model matrix status')
    call assert_matrix_close(csr_to_dense(c), reshape([4.0_dp, 1.0_dp, 0.0_dp, 1.0_dp, 3.0_dp, 2.0_dp, &
        0.0_dp, 2.0_dp, 5.0_dp], [3, 3]), tol, 'symmetric Harwell-Boeing matrix')
    hb_rhs = hb_model_response(hb, info)
    call assert_true(info == -1 .and. all(shape(hb_rhs) == [0, 0]), 'Harwell-Boeing absent rhs')

    print '(a)', 'test_sparsem: PASS'

contains


    subroutine test_ng_numeric_levels()
        real(dp), allocatable :: matrix(:,:)
        real(dp), allocatable :: rhs(:)
        real(dp), allocatable :: solution(:)
        type(csr_matrix) :: sparse
        type(csr_cholesky) :: reference_factor
        type(csr_cholesky) :: test_factor
        integer :: case_id
        integer :: column
        integer :: level_index
        integer :: n
        integer :: row
        integer :: status
        integer, parameter :: levels(3) = [2, 4, 8]

        do case_id = 1, 12
            n = 5 + mod(case_id, 9)
            allocate(matrix(n, n), source=0.0_dp)
            do row = 1, n
                matrix(row, row) = real(n + 4, dp)
                do column = row + 1, n
                    if (mod(11 * row + 17 * column + 3 * case_id + row * column, 13) < 4) then
                        matrix(row, column) = 0.02_dp * real(1 + mod(row + 2 * column + case_id, 7), dp)
                        matrix(column, row) = matrix(row, column)
                    end if
                end do
            end do
            sparse = csr_from_dense(matrix)
            reference_factor = csr_chol(sparse, status, level=1, cache_kb=64)
            call assert_true(status == 0, 'Ng-Peyton level suite reference factorization')
            allocate(rhs(n))
            rhs = [(real(row, dp), row = 1, n)]

            do level_index = 1, size(levels)
                test_factor = csr_chol(sparse, status, level=levels(level_index), cache_kb=64)
                call assert_true(status == 0, 'Ng-Peyton level suite factorization')
                call assert_true(size(test_factor%lnz) == size(reference_factor%lnz), &
                    'Ng-Peyton level suite factor storage size')
                if (size(test_factor%lnz) > 0) then
                    call assert_true(maxval(abs(test_factor%lnz - reference_factor%lnz)) < 5.0e-12_dp, &
                        'Ng-Peyton numerical levels agree')
                end if
                solution = chol_solve(test_factor, rhs, status)
                call assert_true(status == 0, 'Ng-Peyton level suite solve')
                call assert_vector_close(matmul(matrix, solution), rhs, 2.0e-10_dp, &
                    'Ng-Peyton level suite residual')
            end do
            deallocate(matrix, rhs)
        end do
    end subroutine test_ng_numeric_levels

    subroutine assert_true(condition, label)
        logical, intent(in) :: condition !! Condition that must evaluate true for the test to pass.
        character(len=*), intent(in) :: label !! Short description printed when the assertion fails.
        if (.not. condition) then
            print '(a)', 'FAIL: ' // trim(label)
            error stop 1
        end if
    end subroutine assert_true

    subroutine assert_close(actual, expected, atol, label)
        real(dp), intent(in) :: actual !! Scalar result produced by the routine under test.
        real(dp), intent(in) :: expected !! Reference scalar value.
        real(dp), intent(in) :: atol !! Maximum allowed absolute error.
        character(len=*), intent(in) :: label !! Short description printed when the assertion fails.
        if (abs(actual - expected) > atol) then
            print '(a,2(1x,es24.15))', 'FAIL: ' // trim(label), actual, expected
            error stop 1
        end if
    end subroutine assert_close

    subroutine assert_vector_close(actual, expected, atol, label)
        real(dp), intent(in) :: actual(:) !! Vector result produced by the routine under test.
        real(dp), intent(in) :: expected(:) !! Reference vector values.
        real(dp), intent(in) :: atol !! Maximum allowed elementwise absolute error.
        character(len=*), intent(in) :: label !! Short description printed when the assertion fails.
        call assert_true(size(actual) == size(expected), trim(label) // ' size')
        if (size(actual) > 0) then
            if (maxval(abs(actual - expected)) > atol) then
                print '(a)', 'FAIL: ' // trim(label)
                error stop 1
            end if
        end if
    end subroutine assert_vector_close

    subroutine assert_matrix_close(actual, expected, atol, label)
        real(dp), intent(in) :: actual(:,:) !! Matrix result produced by the routine under test.
        real(dp), intent(in) :: expected(:,:) !! Reference matrix values.
        real(dp), intent(in) :: atol !! Maximum allowed elementwise absolute error.
        character(len=*), intent(in) :: label !! Short description printed when the assertion fails.
        call assert_true(all(shape(actual) == shape(expected)), trim(label) // ' shape')
        if (size(actual) > 0) then
            if (maxval(abs(actual - expected)) > atol) then
                print '(a)', 'FAIL: ' // trim(label)
                error stop 1
            end if
        end if
    end subroutine assert_matrix_close

end program test_sparsem
