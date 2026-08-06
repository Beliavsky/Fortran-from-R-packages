! SPDX-License-Identifier: GPL-2.0-or-later
module ltsa_likelihood
    use ltsa_kinds, only : dp, pi
    use ltsa_status, only : ltsa_error, ltsa_success, ltsa_invalid_input, set_error
    use ltsa_types, only : exact_likelihood_result
    use ltsa_toeplitz, only : trench_quadratic_logdet
    implicit none
    private

    public :: trench_loglikelihood, exact_loglikelihood

contains

    function trench_loglikelihood(r, z, error) result(value)
        real(dp), intent(in) :: r(:), z(:)
        type(ltsa_error), intent(out), optional :: error
        real(dp) :: value, quadratic, logdet
        type(ltsa_error) :: local_error
        call trench_quadratic_logdet(r, z, quadratic, logdet, local_error)
        if (.not. local_error%ok()) then
            value = -huge(1.0_dp)
        else
            value = -0.5_dp*real(size(z),dp)*log(quadratic/real(size(z),dp))-0.5_dp*logdet
        end if
        if (present(error)) error = local_error
    end function trench_loglikelihood

    function exact_loglikelihood(r, z, innovation_variance) result(result_value)
        real(dp), intent(in) :: r(:), z(:)
        logical, intent(in), optional :: innovation_variance
        type(exact_likelihood_result) :: result_value
        logical :: estimate_variance
        real(dp) :: quadratic, logdet, r0
        estimate_variance = .true.
        if (present(innovation_variance)) estimate_variance = innovation_variance
        if (size(r) /= size(z) .or. size(r) < 1 .or. r(1) <= 0.0_dp) then
            call set_error(result_value%error, ltsa_invalid_input, 'r and z must be equal-length and r(1) positive')
            return
        end if
        r0 = r(1)
        call trench_quadratic_logdet(r, z, quadratic, logdet, result_value%error)
        if (.not. result_value%error%ok()) return
        if (estimate_variance) then
            result_value%log_likelihood = -0.5_dp*real(size(z),dp)*(1.0_dp+log(2.0_dp*pi)) &
                                          -0.5_dp*(real(size(z),dp)*log(quadratic/real(size(z),dp))+logdet)
            result_value%sigma_sq = quadratic/(real(size(z),dp)*r0)
        else
            result_value%log_likelihood = -0.5_dp*real(size(z),dp)*log(2.0_dp*pi*r0) &
                                          -0.5_dp*logdet-0.5_dp*quadratic/r0
            result_value%sigma_sq = r0
        end if
    end function exact_loglikelihood

end module ltsa_likelihood
