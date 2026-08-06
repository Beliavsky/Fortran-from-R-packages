! SPDX-License-Identifier: GPL-2.0-or-later
module ltsa_durbin_levinson
    use ltsa_kinds, only : dp
    use ltsa_status, only : ltsa_error, ltsa_success, ltsa_invalid_input, ltsa_not_positive_definite, set_error
    use ltsa_types, only : dl_ar_result
    use ltsa_random, only : ltsa_normal
    implicit none
    private

    public :: durbin_levinson_table, dl_acf_to_ar, dl_residuals, dl_loglikelihood, dl_simulate

contains

    subroutine durbin_levinson_table(r, phi_table, pacf, variances, error)
        real(dp), intent(in) :: r(:)
        real(dp), allocatable, intent(out) :: phi_table(:,:), pacf(:), variances(:)
        type(ltsa_error), intent(out) :: error
        integer :: j, k, n
        real(dp) :: reflection, numerator, tol
        error%code = ltsa_success
        error%message = ''
        n = size(r)
        if (n < 1 .or. r(1) <= 0.0_dp) then
            allocate(phi_table(0,0), pacf(0), variances(0))
            call set_error(error, ltsa_invalid_input, 'autocovariance sequence must start with a positive variance')
            return
        end if
        allocate(phi_table(max(0,n-1),max(0,n-1)), source=0.0_dp)
        allocate(pacf(max(0,n-1)), source=0.0_dp)
        allocate(variances(n), source=0.0_dp)
        variances(1) = r(1)
        tol = epsilon(1.0_dp)*max(1.0_dp,abs(r(1)))*real(max(1,n),dp)
        do k = 1, n-1
            numerator = r(k+1)
            do j = 1, k-1
                numerator = numerator-phi_table(k-1,j)*r(k-j+1)
            end do
            if (variances(k) <= tol) then
                call set_error(error, ltsa_not_positive_definite, 'autocovariance sequence is not positive definite')
                return
            end if
            reflection = numerator/variances(k)
            pacf(k) = reflection
            do j = 1, k-1
                phi_table(k,j) = phi_table(k-1,j)-reflection*phi_table(k-1,k-j)
            end do
            phi_table(k,k) = reflection
            variances(k+1) = variances(k)*(1.0_dp-reflection*reflection)
            if (variances(k+1) <= tol) then
                call set_error(error, ltsa_not_positive_definite, 'autocovariance sequence is not positive definite')
                return
            end if
        end do
    end subroutine durbin_levinson_table

    function dl_acf_to_ar(r) result(result_value)
        real(dp), intent(in) :: r(:)
        type(dl_ar_result) :: result_value
        real(dp), allocatable :: acvf(:), table(:,:), pacf(:), variances(:)
        integer :: n
        n = size(r)
        allocate(acvf(n+1))
        acvf(1) = 1.0_dp
        if (n > 0) acvf(2:) = r
        if (n > 0 .and. abs(r(1)) >= 1.0_dp) then
            allocate(result_value%phi(0), result_value%pacf(0), result_value%prediction_variance(0))
            call set_error(result_value%error, ltsa_invalid_input, 'lag-one autocorrelation must have absolute value below one')
            return
        end if
        call durbin_levinson_table(acvf, table, pacf, variances, result_value%error)
        if (.not. result_value%error%ok()) then
            allocate(result_value%phi(0), result_value%pacf(0), result_value%prediction_variance(0))
            return
        end if
        allocate(result_value%phi(n), result_value%pacf(n), result_value%prediction_variance(n))
        if (n > 0) then
            result_value%phi = table(n,1:n)
            result_value%pacf = pacf
            result_value%prediction_variance = variances(2:)
        end if
    end function dl_acf_to_ar

    subroutine dl_residuals(r, z, residuals, error, standardized, prediction_variances)
        real(dp), intent(in) :: r(:), z(:)
        real(dp), allocatable, intent(out) :: residuals(:)
        type(ltsa_error), intent(out) :: error
        logical, intent(in), optional :: standardized
        real(dp), allocatable, intent(out), optional :: prediction_variances(:)
        real(dp), allocatable :: table(:,:), pacf(:), variances(:)
        real(dp) :: prediction
        logical :: standardize
        integer :: j, t, n
        standardize = .true.
        if (present(standardized)) standardize = standardized
        if (size(r) /= size(z) .or. size(z) < 1) then
            allocate(residuals(0))
            if (present(prediction_variances)) allocate(prediction_variances(0))
            call set_error(error, ltsa_invalid_input, 'r and z must be nonempty and have equal lengths')
            return
        end if
        call durbin_levinson_table(r, table, pacf, variances, error)
        if (.not. error%ok()) then
            allocate(residuals(0))
            if (present(prediction_variances)) allocate(prediction_variances(0))
            return
        end if
        n = size(z)
        allocate(residuals(n))
        residuals(1) = z(1)
        do t = 2, n
            prediction = 0.0_dp
            do j = 1, t-1
                prediction = prediction+table(t-1,j)*z(t-j)
            end do
            residuals(t) = z(t)-prediction
        end do
        if (standardize) residuals = residuals/sqrt(variances)
        if (present(prediction_variances)) then
            allocate(prediction_variances(n))
            prediction_variances = variances
        end if
    end subroutine dl_residuals

    function dl_loglikelihood(r, z, error) result(value)
        real(dp), intent(in) :: r(:), z(:)
        type(ltsa_error), intent(out), optional :: error
        real(dp) :: value
        type(ltsa_error) :: local_error
        real(dp), allocatable :: residuals(:), variances(:)
        real(dp) :: s
        call dl_residuals(r, z, residuals, local_error, standardized=.false., prediction_variances=variances)
        if (.not. local_error%ok()) then
            value = -huge(1.0_dp)
            if (present(error)) error = local_error
            return
        end if
        s = sum(residuals*residuals/variances)
        value = -0.5_dp*real(size(z),dp)*log(s/real(size(z),dp))-0.5_dp*sum(log(variances))
        if (present(error)) error = local_error
    end function dl_loglikelihood

    subroutine dl_simulate(n, r, z, error, innovations)
        integer, intent(in) :: n
        real(dp), intent(in) :: r(:)
        real(dp), allocatable, intent(out) :: z(:)
        type(ltsa_error), intent(out) :: error
        real(dp), intent(in), optional :: innovations(:)
        real(dp), allocatable :: rwork(:), table(:,:), pacf(:), variances(:), a(:)
        real(dp) :: prediction
        integer :: j, t, nr
        nr = size(r)
        if (n < 1 .or. nr < 1) then
            allocate(z(0))
            call set_error(error, ltsa_invalid_input, 'n and the autocovariance sequence must be positive')
            return
        end if
        if (present(innovations)) then
            if (size(innovations) < n) then
                allocate(z(0))
                call set_error(error, ltsa_invalid_input, 'too few supplied innovations')
                return
            end if
        end if
        allocate(rwork(n), source=0.0_dp)
        rwork(1:min(n,nr)) = r(1:min(n,nr))
        call durbin_levinson_table(rwork, table, pacf, variances, error)
        if (.not. error%ok()) then
            allocate(z(0))
            return
        end if
        allocate(a(n), z(n))
        if (present(innovations)) then
            a = innovations(1:n)
        else
            do t = 1, n
                a(t) = ltsa_normal()
            end do
        end if
        z(1) = a(1)*sqrt(variances(1))
        do t = 2, n
            prediction = 0.0_dp
            do j = 1, t-1
                prediction = prediction+table(t-1,j)*z(t-j)
            end do
            z(t) = prediction+a(t)*sqrt(variances(t))
        end do
    end subroutine dl_simulate

end module ltsa_durbin_levinson
