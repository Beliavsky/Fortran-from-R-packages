! SPDX-License-Identifier: GPL-2.0-or-later
! Derived from MCMCglmm 2.36 computational code; see NOTICE.md and upstream/.
module mcmcglmm_matrix
    use r_kinds, only : dp
    use r_linalg, only : cholesky_factor, inverse_matrix, solve_spd, spd_inverse_logdet, symmetric_eigen
    use mcmcglmm_rng, only : rng_state, rng_normal
    implicit none
    private

    real(dp), parameter :: log_two_pi = 1.8378770664093454835606594728112353_dp

    public :: covariance_to_correlation
    public :: block_diagonal
    public :: commutation_matrix
    public :: kronecker_product
    public :: mvn_log_density
    public :: sample_mvn_covariance
    public :: sample_mvn_precision
    public :: symmetric_eigenvalues

contains

    pure subroutine covariance_to_correlation(covariance, correlation, info)
        real(dp), intent(in) :: covariance(:, :) !! Square covariance matrix with positive diagonal entries.
        real(dp), allocatable, intent(out) :: correlation(:, :) !! Allocated correlation matrix of the same shape.
        integer, intent(out) :: info !! Zero on success; nonzero for nonsquare input or a nonpositive diagonal entry.
        real(dp), allocatable :: scale(:)
        integer :: i
        integer :: j
        integer :: n

        info = 0
        n = size(covariance, 1)
        if (size(covariance, 2) /= n) then
            allocate(correlation(0, 0))
            info = 1
            return
        end if
        allocate(correlation(n, n), scale(n))
        do i = 1, n
            if (covariance(i, i) <= 0.0_dp) then
                correlation = 0.0_dp
                info = 2
                return
            end if
            scale(i) = sqrt(covariance(i, i))
        end do
        do j = 1, n
            do i = 1, n
                correlation(i, j) = covariance(i, j) / (scale(i) * scale(j))
            end do
        end do
    end subroutine covariance_to_correlation

    pure subroutine kronecker_product(a, b, result_matrix)
        real(dp), intent(in) :: a(:, :) !! Left matrix factor with arbitrary rectangular shape.
        real(dp), intent(in) :: b(:, :) !! Right matrix factor with arbitrary rectangular shape.
        real(dp), allocatable, intent(out) :: result_matrix(:, :) !! Allocated Kronecker product A tensor B.
        integer :: i
        integer :: j
        integer :: mb
        integer :: nb

        mb = size(b, 1)
        nb = size(b, 2)
        allocate(result_matrix(size(a, 1) * mb, size(a, 2) * nb))
        do j = 1, size(a, 2)
            do i = 1, size(a, 1)
                result_matrix((i - 1) * mb + 1:i * mb, (j - 1) * nb + 1:j * nb) = a(i, j) * b
            end do
        end do
    end subroutine kronecker_product

    pure subroutine commutation_matrix(m, n, matrix_k)
        integer, intent(in) :: m !! Positive row dimension of an m by n matrix before vectorization.
        integer, intent(in) :: n !! Positive column dimension of an m by n matrix before vectorization.
        real(dp), allocatable, intent(out) :: matrix_k(:, :) !! Permutation matrix satisfying K vec(A) = vec(transpose(A)).
        integer :: i
        integer :: j
        integer :: source_index
        integer :: target_index

        allocate(matrix_k(m * n, m * n))
        matrix_k = 0.0_dp
        do j = 1, n
            do i = 1, m
                source_index = i + (j - 1) * m
                target_index = j + (i - 1) * n
                matrix_k(target_index, source_index) = 1.0_dp
            end do
        end do
    end subroutine commutation_matrix

    pure subroutine block_diagonal(blocks, block_rows, block_cols, result_matrix, info)
        real(dp), intent(in) :: blocks(:, :, :) !! Dense blocks in the active leading slice of each third-dimension plane.
        integer, intent(in) :: block_rows(:) !! Number of active rows in each block plane.
        integer, intent(in) :: block_cols(:) !! Number of active columns in each block plane.
        real(dp), allocatable, intent(out) :: result_matrix(:, :) !! Allocated block-diagonal matrix.
        integer, intent(out) :: info !! Zero on success; nonzero for inconsistent block metadata.
        integer :: k
        integer :: row_offset
        integer :: col_offset
        integer :: total_rows
        integer :: total_cols

        info = 0
        if (size(block_rows) /= size(block_cols) .or. size(block_rows) /= size(blocks, 3)) then
            allocate(result_matrix(0, 0))
            info = 1
            return
        end if
        if (any(block_rows < 0) .or. any(block_cols < 0)) then
            allocate(result_matrix(0, 0))
            info = 2
            return
        end if
        if (any(block_rows > size(blocks, 1)) .or. any(block_cols > size(blocks, 2))) then
            allocate(result_matrix(0, 0))
            info = 3
            return
        end if

        total_rows = sum(block_rows)
        total_cols = sum(block_cols)
        allocate(result_matrix(total_rows, total_cols))
        result_matrix = 0.0_dp
        row_offset = 0
        col_offset = 0
        do k = 1, size(block_rows)
            result_matrix(row_offset + 1:row_offset + block_rows(k), &
                          col_offset + 1:col_offset + block_cols(k)) = &
                blocks(1:block_rows(k), 1:block_cols(k), k)
            row_offset = row_offset + block_rows(k)
            col_offset = col_offset + block_cols(k)
        end do
    end subroutine block_diagonal

    pure subroutine mvn_log_density(x, mean_value, covariance, log_density, info)
        real(dp), intent(in) :: x(:) !! Evaluation vector of length n.
        real(dp), intent(in) :: mean_value(:) !! Mean vector of length n.
        real(dp), intent(in) :: covariance(:, :) !! Symmetric positive-definite n by n covariance matrix.
        real(dp), intent(out) :: log_density !! Multivariate-normal log density, including the normalization constant.
        integer, intent(out) :: info !! Zero on success; nonzero for shape or factorization failure.
        real(dp), allocatable :: inverse(:, :)
        real(dp), allocatable :: deviation(:)
        real(dp) :: logdet
        integer :: n

        info = 0
        n = size(x)
        if (size(mean_value) /= n .or. size(covariance, 1) /= n .or. size(covariance, 2) /= n) then
            log_density = -huge(1.0_dp)
            info = 1
            return
        end if
        call spd_inverse_logdet(covariance, inverse, logdet, info)
        if (info /= 0) then
            log_density = -huge(1.0_dp)
            return
        end if
        deviation = x - mean_value
        log_density = -0.5_dp * (real(n, dp) * log_two_pi + logdet + dot_product(deviation, matmul(inverse, deviation)))
    end subroutine mvn_log_density

    pure subroutine sample_mvn_covariance(state, mean_value, covariance, sample, info)
        type(rng_state), intent(inout) :: state !! Generator state used for the standard-normal innovations.
        real(dp), intent(in) :: mean_value(:) !! Mean vector of length n.
        real(dp), intent(in) :: covariance(:, :) !! Symmetric positive-definite n by n covariance matrix.
        real(dp), allocatable, intent(out) :: sample(:) !! Allocated multivariate-normal sample of length n.
        integer, intent(out) :: info !! Zero on success; nonzero for shape or Cholesky failure.
        real(dp), allocatable :: factor(:, :)
        real(dp), allocatable :: z(:)
        integer :: i
        integer :: n

        info = 0
        n = size(mean_value)
        allocate(sample(n))
        if (size(covariance, 1) /= n .or. size(covariance, 2) /= n) then
            sample = 0.0_dp
            info = 1
            return
        end if
        call cholesky_factor(covariance, factor, info)
        if (info /= 0) then
            sample = 0.0_dp
            return
        end if
        allocate(z(n))
        do i = 1, n
            call rng_normal(state, z(i))
        end do
        sample = mean_value + matmul(factor, z)
    end subroutine sample_mvn_covariance

    pure subroutine sample_mvn_precision(state, rhs, precision, sample, mean_value, info)
        type(rng_state), intent(inout) :: state !! Generator state used for the conditional Gaussian innovation.
        real(dp), intent(in) :: rhs(:) !! Precision-weighted mean vector P times mu.
        real(dp), intent(in) :: precision(:, :) !! Symmetric positive-definite precision matrix P.
        real(dp), allocatable, intent(out) :: sample(:) !! Allocated sample from N(P^{-1} rhs, P^{-1}).
        real(dp), allocatable, intent(out) :: mean_value(:) !! Allocated conditional mean P^{-1} rhs.
        integer, intent(out) :: info !! Zero on success; nonzero for shape or linear-algebra failure.
        real(dp), allocatable :: covariance(:, :)
        real(dp) :: logdet

        info = 0
        if (size(precision, 1) /= size(rhs) .or. size(precision, 2) /= size(rhs)) then
            allocate(sample(0), mean_value(0))
            info = 1
            return
        end if
        allocate(mean_value(size(rhs)))
        call solve_spd(precision, rhs, mean_value, info)
        if (info /= 0) then
            allocate(sample(0))
            return
        end if
        call spd_inverse_logdet(precision, covariance, logdet, info)
        if (info /= 0) then
            allocate(sample(0))
            return
        end if
        call sample_mvn_covariance(state, mean_value, covariance, sample, info)
    end subroutine sample_mvn_precision

    subroutine symmetric_eigenvalues(matrix_a, values, info, descending)
        real(dp), intent(in) :: matrix_a(:, :) !! Symmetric square matrix whose eigenvalues are requested.
        real(dp), allocatable, intent(out) :: values(:) !! Allocated eigenvalues ordered according to descending when successful.
        integer, intent(out) :: info !! Zero on success; nonzero on eigensolver failure.
        logical, intent(in), optional :: descending !! If present and true, order eigenvalues from largest to smallest.
        real(dp), allocatable :: vectors(:, :)
        logical :: order_descending

        order_descending = .false.
        if (present(descending)) order_descending = descending
        call symmetric_eigen(matrix_a, values, vectors, info, descending=order_descending)
    end subroutine symmetric_eigenvalues

end module mcmcglmm_matrix
