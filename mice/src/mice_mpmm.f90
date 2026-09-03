! SPDX-License-Identifier: GPL-2.0-or-later
! Derived from R package mice 3.19.0 by Stef van Buuren, Karin Groothuis-Oudshoorn,
! and mice contributors; see NOTICE.md and PROVENANCE.md for attribution.
! Multivariate predictive mean matching derived from mice.impute.mpmm.R.
module mice_mpmm
    use r_kinds, only : dp
    use r_linalg, only : inverse_matrix
    use mice_impute_continuous, only : impute_pmm
    use mice_rng, only : mice_rng_state
    use mice_status, only : mice_ok, mice_invalid_argument, mice_invalid_shape, mice_no_observed, mice_singular
    implicit none
    private

    public :: impute_mpmm

contains

    subroutine impute_mpmm(y, observed, x, where, rng, imputed, info, donors, ridge)
        real(dp), intent(in) :: y(:, :) !! Multivariate incomplete response block whose rows share a missingness pattern.
        logical, intent(in) :: observed(:) !! True for rows where all response components are observed and may donate.
        real(dp), intent(in) :: x(:, :) !! Complete predictor matrix used in PMM for the canonical response score.
        logical, intent(in) :: where(:) !! True for rows requiring imputation of the entire response block.
        type(mice_rng_state), intent(inout) :: rng !! RNG state used by predictive mean matching.
        real(dp), allocatable, intent(out) :: imputed(:, :) !! Donor response rows corresponding to true entries of `where`.
        integer, intent(out) :: info !! `mice_ok` on success or a package status code on failure.
        integer, intent(in), optional :: donors !! PMM donor-pool size; default five.
        real(dp), intent(in), optional :: ridge !! PMM regression ridge; default `1e-5`.

        integer, allocatable :: obs_index(:), mis_index(:)
        real(dp), allocatable :: yobs(:, :), xobs(:, :), syy(:, :), sxx(:, :), syx(:, :), sxy(:, :)
        real(dp), allocatable :: inv_yy(:, :), inv_xx(:, :), canonical(:, :), direction(:), score(:), score_imp(:)
        real(dp) :: ridge_value
        integer :: donor_count, i, j, k, nmis, nobs, ny, status

        if (size(observed) /= size(y, 1) .or. size(where) /= size(y, 1) .or. size(x, 1) /= size(y, 1)) then
            info = mice_invalid_shape
            return
        end if
        nobs = count(observed)
        nmis = count(where)
        ny = size(y, 2)
        if (nobs < 2) then
            info = mice_no_observed
            return
        end if
        if (ny < 1 .or. size(x, 2) < 1) then
            info = mice_invalid_argument
            return
        end if
        donor_count = 5
        if (present(donors)) donor_count = donors
        ridge_value = 1.0e-5_dp
        if (present(ridge)) ridge_value = ridge
        if (donor_count < 1 .or. ridge_value < 0.0_dp) then
            info = mice_invalid_argument
            return
        end if

        allocate(obs_index(nobs), mis_index(nmis), yobs(nobs, ny), xobs(nobs, size(x, 2)))
        j = 0
        do i = 1, size(y, 1)
            if (.not. observed(i)) cycle
            j = j + 1
            obs_index(j) = i
            yobs(j, :) = y(i, :)
            xobs(j, :) = x(i, :)
        end do
        j = 0
        do i = 1, size(y, 1)
            if (.not. where(i)) cycle
            j = j + 1
            mis_index(j) = i
        end do
        call sample_covariance(yobs, yobs, syy)
        call sample_covariance(xobs, xobs, sxx)
        call sample_covariance(yobs, xobs, syx)
        sxy = transpose(syx)
        call inverse_matrix(syy, inv_yy, status)
        if (status /= 0) then
            call add_covariance_ridge(syy, ridge_value)
            call inverse_matrix(syy, inv_yy, status)
        end if
        if (status /= 0) then
            info = mice_singular
            return
        end if
        call inverse_matrix(sxx, inv_xx, status)
        if (status /= 0) then
            call add_covariance_ridge(sxx, ridge_value)
            call inverse_matrix(sxx, inv_xx, status)
        end if
        if (status /= 0) then
            info = mice_singular
            return
        end if
        canonical = matmul(matmul(matmul(inv_yy, syx), inv_xx), sxy)
        call dominant_eigenvector(canonical, direction)
        allocate(score(size(y, 1)))
        score = matmul(y, direction)
        call impute_pmm(score, observed, x, where, rng, score_imp, info, donors=donor_count, ridge=ridge_value)
        if (info /= mice_ok) return

        allocate(imputed(nmis, ny))
        do i = 1, nmis
            k = 1
            do j = 2, nobs
                if (abs(score(obs_index(j)) - score_imp(i)) < abs(score(obs_index(k)) - score_imp(i))) k = j
            end do
            imputed(i, :) = y(obs_index(k), :)
        end do
        info = mice_ok
    end subroutine impute_mpmm

    pure subroutine sample_covariance(a, b, covariance)
        real(dp), intent(in) :: a(:, :) !! First observed data matrix with observations in rows.
        real(dp), intent(in) :: b(:, :) !! Second observed data matrix sharing the rows of `a`.
        real(dp), allocatable, intent(out) :: covariance(:, :) !! Sample cross-covariance using denominator `n-1`.

        real(dp), allocatable :: ca(:, :), cb(:, :)
        real(dp) :: ma, mb
        integer :: j, n

        n = size(a, 1)
        allocate(ca(n, size(a, 2)), cb(n, size(b, 2)), covariance(size(a, 2), size(b, 2)))
        do j = 1, size(a, 2)
            ma = sum(a(:, j)) / real(n, dp)
            ca(:, j) = a(:, j) - ma
        end do
        do j = 1, size(b, 2)
            mb = sum(b(:, j)) / real(n, dp)
            cb(:, j) = b(:, j) - mb
        end do
        covariance = matmul(transpose(ca), cb) / real(max(n - 1, 1), dp)
    end subroutine sample_covariance

    pure subroutine add_covariance_ridge(matrix, ridge)
        real(dp), intent(inout) :: matrix(:, :) !! Symmetric covariance matrix regularized in place.
        real(dp), intent(in), value :: ridge !! Relative ridge magnitude; zero uses a machine-scale fallback.
        real(dp) :: scale
        integer :: j

        scale = max(maxval(abs(matrix)), 1.0_dp)
        do j = 1, min(size(matrix, 1), size(matrix, 2))
            matrix(j, j) = matrix(j, j) + max(ridge * max(abs(matrix(j, j)), scale), &
                                                sqrt(epsilon(1.0_dp)) * scale)
        end do
    end subroutine add_covariance_ridge

    pure subroutine dominant_eigenvector(matrix, vector)
        real(dp), intent(in) :: matrix(:, :) !! Square canonical-direction matrix whose dominant right eigenvector is required.
        real(dp), allocatable, intent(out) :: vector(:) !! Unit-length dominant direction, with deterministic sign.

        real(dp), allocatable :: next(:)
        real(dp) :: norm_value
        integer :: i, iter, n, start

        n = size(matrix, 1)
        allocate(vector(n), next(n))
        if (maxval(abs(matrix)) <= sqrt(epsilon(1.0_dp))) then
            vector = 0.0_dp
            vector(1) = 1.0_dp
            return
        end if
        start = 1
        do i = 2, n
            if (sum(abs(matrix(:, i))) > sum(abs(matrix(:, start)))) start = i
        end do
        vector = matrix(:, start)
        norm_value = sqrt(dot_product(vector, vector))
        if (norm_value <= tiny(1.0_dp)) then
            vector = 1.0_dp / sqrt(real(n, dp))
        else
            vector = vector / norm_value
        end if
        do iter = 1, 500
            next = matmul(matrix, vector)
            norm_value = sqrt(dot_product(next, next))
            if (norm_value <= tiny(1.0_dp)) exit
            next = next / norm_value
            if (maxval(abs(next - vector)) <= 1.0e-10_dp .or. maxval(abs(next + vector)) <= 1.0e-10_dp) then
                vector = next
                exit
            end if
            vector = next
        end do
        do i = 1, n
            if (abs(vector(i)) > sqrt(epsilon(1.0_dp))) then
                if (vector(i) < 0.0_dp) vector = -vector
                exit
            end if
        end do
    end subroutine dominant_eigenvector

end module mice_mpmm
