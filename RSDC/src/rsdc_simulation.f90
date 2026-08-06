! SPDX-License-Identifier: GPL-3.0-only
module rsdc_simulation
    use rsdc_kinds, only: dp
    use rsdc_types, only: rsdc_simulation_result
    use rsdc_parameters, only: transition_at
    use rsdc_random, only: seed_rsdc, random_multivariate_normal, categorical_draw
    use rsdc_linalg, only: is_positive_definite
    implicit none
    private
    public :: rsdc_simulate, rsdc_simulate_fixed

contains

    subroutine rsdc_simulate(nobs, x, beta, mu, covariance, result, seed, ok)
        integer, intent(in) :: nobs
        real(dp), intent(in) :: x(:, :), beta(:, :), mu(:, :), covariance(:, :, :)
        type(rsdc_simulation_result), intent(out) :: result
        integer, intent(in), optional :: seed
        logical, intent(out), optional :: ok
        integer :: t, n, k, state
        real(dp), allocatable :: pmat(:, :), init(:)
        logical :: good
        n = size(mu, 1); k = size(mu, 2)
        good = nobs >= 1 .and. size(x, 1) == nobs .and. size(beta, 1) == n
        good = good .and. size(covariance, 1) == k .and. size(covariance, 2) == k
        good = good .and. size(covariance, 3) == n
        if (.not. good) then
            if (present(ok)) ok = .false.
            return
        end if
        do state = 1, n
            if (.not. is_positive_definite(covariance(:, :, state), 1.0e-10_dp)) then
                if (present(ok)) ok = .false.
                return
            end if
        end do
        if (present(seed)) call seed_rsdc(seed)
        allocate(result%states(nobs), result%observations(nobs, k))
        allocate(result%transition_matrices(n, n, nobs), pmat(n, n), init(n))
        result%transition_matrices = 0.0_dp
        call transition_at(beta, x(1, :), pmat)
        init = matmul(transpose(pmat), [(1.0_dp / real(n, dp), state=1,n)])
        result%states(1) = categorical_draw(init)
        call random_multivariate_normal(mu(result%states(1), :), &
            covariance(:, :, result%states(1)), result%observations(1, :), good)
        if (.not. good) then
            if (present(ok)) ok = .false.
            return
        end if
        do t = 2, nobs
            call transition_at(beta, x(t, :), pmat)
            result%transition_matrices(:, :, t) = pmat
            state = result%states(t - 1)
            result%states(t) = categorical_draw(pmat(state, :))
            call random_multivariate_normal(mu(result%states(t), :), &
                covariance(:, :, result%states(t)), result%observations(t, :), good)
            if (.not. good) exit
        end do
        if (present(ok)) ok = good
    end subroutine rsdc_simulate

    subroutine rsdc_simulate_fixed(nobs, pmat, mu, covariance, result, seed, ok)
        integer, intent(in) :: nobs
        real(dp), intent(in) :: pmat(:, :), mu(:, :), covariance(:, :, :)
        type(rsdc_simulation_result), intent(out) :: result
        integer, intent(in), optional :: seed
        logical, intent(out), optional :: ok
        integer :: t, n, k, state
        real(dp), allocatable :: init(:)
        logical :: good
        n = size(mu, 1); k = size(mu, 2)
        good = nobs >= 1 .and. size(pmat, 1) == n .and. size(pmat, 2) == n
        if (.not. good) then
            if (present(ok)) ok = .false.
            return
        end if
        if (present(seed)) call seed_rsdc(seed)
        allocate(result%states(nobs), result%observations(nobs, k))
        allocate(result%transition_matrices(n, n, nobs), init(n))
        result%transition_matrices = 0.0_dp
        init = matmul(transpose(pmat), [(1.0_dp / real(n, dp), state=1,n)])
        result%states(1) = categorical_draw(init)
        call random_multivariate_normal(mu(result%states(1), :), &
            covariance(:, :, result%states(1)), result%observations(1, :), good)
        do t = 2, nobs
            result%transition_matrices(:, :, t) = pmat
            state = result%states(t - 1)
            result%states(t) = categorical_draw(pmat(state, :))
            call random_multivariate_normal(mu(result%states(t), :), &
                covariance(:, :, result%states(t)), result%observations(t, :), good)
            if (.not. good) exit
        end do
        if (present(ok)) ok = good
    end subroutine rsdc_simulate_fixed
end module rsdc_simulation
