! SPDX-License-Identifier: GPL-2.0-or-later
module ltsa_forecast
    use ltsa_kinds, only : dp
    use ltsa_status, only : ltsa_error, ltsa_success, ltsa_invalid_input, ltsa_not_positive_definite, set_error
    use ltsa_types, only : dl_ar_result, forecast_result
    use ltsa_linalg, only : toeplitz_matrix
    use ltsa_durbin_levinson, only : dl_acf_to_ar
    use ltsa_toeplitz, only : trench_inverse, toeplitz_inverse_update
    use ltsa_arma, only : ar_to_ma
    implicit none
    private

    public :: prediction_variance, trench_forecast

contains

    subroutine prediction_variance(r, max_lead, variances, error, use_durbin_levinson)
        real(dp), intent(in) :: r(:)
        integer, intent(in) :: max_lead
        real(dp), allocatable, intent(out) :: variances(:)
        type(ltsa_error), intent(out) :: error
        logical, intent(in), optional :: use_durbin_levinson
        logical :: use_dl
        type(dl_ar_result) :: ar_result
        real(dp), allocatable :: psi(:), g(:,:), covariance(:,:), inverse(:,:), temp(:)
        integer :: i, j, n
        use_dl = .true.
        if (present(use_durbin_levinson)) use_dl = use_durbin_levinson
        error%code = ltsa_success
        error%message = ''
        if (size(r) < 1 .or. r(1) <= 0.0_dp .or. max_lead < 1) then
            allocate(variances(0))
            call set_error(error, ltsa_invalid_input, 'r(1) and max_lead must be positive')
            return
        end if
        allocate(variances(max_lead))
        if (use_dl) then
            if (size(r) == 1) then
                variances = r(1)
                return
            end if
            ar_result = dl_acf_to_ar(r(2:)/r(1))
            if (.not. ar_result%error%ok()) then
                error = ar_result%error
                deallocate(variances)
                allocate(variances(0))
                return
            end if
            psi = ar_to_ma(ar_result%phi, max_lead-1)
            variances = r(1)*ar_result%prediction_variance(size(ar_result%prediction_variance))*cumsum_square(psi)
        else
            n = size(r)-max_lead
            if (n < 1) then
                deallocate(variances)
                allocate(variances(0))
                call set_error(error, ltsa_invalid_input, 'length(r) must exceed max_lead for exact prediction variances')
                return
            end if
            covariance = toeplitz_matrix(r(1:n))
            call trench_inverse(covariance, inverse, error)
            if (.not. error%ok()) then
                deallocate(variances)
                allocate(variances(0))
                return
            end if
            allocate(g(max_lead,n), temp(n))
            do j = 1, max_lead
                do i = 1, n
                    g(j,i) = r(n+j-i+1)
                end do
                temp = matmul(inverse,g(j,:))
                variances(j) = r(1)-dot_product(g(j,:),temp)
            end do
            if (any(variances <= 0.0_dp)) then
                call set_error(error, ltsa_not_positive_definite, 'autocovariance sequence gives nonpositive forecast variance')
            end if
        end if
    end subroutine prediction_variance

    function cumsum_square(x) result(y)
        real(dp), intent(in) :: x(:)
        real(dp) :: y(size(x))
        integer :: i
        y(1) = x(1)*x(1)
        do i = 2, size(x)
            y(i) = y(i-1)+x(i)*x(i)
        end do
    end function cumsum_square

    function trench_forecast(z, r, mean_value, origin, max_lead, update_algorithm) result(result_value)
        real(dp), intent(in) :: z(:), r(:), mean_value
        integer, intent(in) :: origin, max_lead
        logical, intent(in), optional :: update_algorithm
        type(forecast_result) :: result_value
        logical :: update
        real(dp), allocatable :: centered(:), covariance(:,:), inverse(:,:), g(:), weights(:), updated(:,:)
        integer :: i, j, t, row, nobs, nz
        update = .true.
        if (present(update_algorithm)) update = update_algorithm
        nz = size(z)
        if (origin < 0 .or. origin > nz .or. max_lead < 1 .or. size(r) < nz+max_lead) then
            allocate(result_value%forecasts(0,0), result_value%sd_forecasts(0,0))
            call set_error(result_value%error, ltsa_invalid_input, 'invalid origin, lead, or autocovariance length')
            return
        end if
        allocate(centered(nz))
        centered = z-mean_value
        allocate(result_value%forecasts(nz-origin+1,max_lead))
        allocate(result_value%sd_forecasts(nz-origin+1,max_lead))
        if (origin == 0) then
            result_value%forecasts(1,:) = mean_value
            result_value%sd_forecasts(1,:) = sqrt(r(1))
            if (nz == 0) return
            nobs = 1
            covariance = toeplitz_matrix(r(1:1))
            call trench_inverse(covariance, inverse, result_value%error)
            if (.not. result_value%error%ok()) return
            row = 2
        else
            nobs = origin
            covariance = toeplitz_matrix(r(1:nobs))
            call trench_inverse(covariance, inverse, result_value%error)
            if (.not. result_value%error%ok()) return
            row = 1
        end if
        do t = nobs, nz
            allocate(g(t), weights(t))
            do j = 1, max_lead
                do i = 1, t
                    g(i) = r(t+j-i+1)
                end do
                weights = matmul(inverse,g)
                result_value%forecasts(row,j) = mean_value+dot_product(weights,centered(1:t))
                result_value%sd_forecasts(row,j) = sqrt(max(0.0_dp,r(1)-dot_product(g,weights)))
            end do
            deallocate(g,weights)
            if (t == nz) exit
            if (update) then
                call toeplitz_inverse_update(inverse, r(1:t), r(t+1), updated, result_value%error)
                if (.not. result_value%error%ok()) return
                call move_alloc(updated,inverse)
            else
                covariance = toeplitz_matrix(r(1:t+1))
                call trench_inverse(covariance, inverse, result_value%error)
                if (.not. result_value%error%ok()) return
            end if
            row = row+1
        end do
    end function trench_forecast

end module ltsa_forecast
