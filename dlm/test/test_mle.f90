! SPDX-License-Identifier: GPL-2.0-or-later
module test_mle_types
    use dlm, only : dp, dlm_model, dlm_success, dlm_invalid_argument, dlm_model_builder
    implicit none
    private

    type, extends(dlm_model_builder), public :: iid_variance_builder
        real(dp) :: variance_scale = 1.0_dp
    contains
        procedure, pass :: build => build_iid_variance
    end type iid_variance_builder

contains

    subroutine build_iid_variance(self, parameters, model, info)
        class(iid_variance_builder), intent(in) :: self !! Builder carrying a fixed positive multiplier for observation variance.
        real(dp), intent(in) :: parameters(:) !! One-element vector containing the logarithm of the scaled observation variance.
        type(dlm_model), intent(out) :: model !! One-state model whose observations are independent zero-mean Gaussian values.
        integer, intent(out) :: info !! Zero for a one-element parameter vector or `dlm_invalid_argument` otherwise.

        real(dp) :: observation_variance

        info = dlm_success
        if (size(parameters) /= 1 .or. self%variance_scale <= 0.0_dp) then
            info = dlm_invalid_argument
            return
        end if
        observation_variance = self%variance_scale * exp(parameters(1))
        allocate(model%m0(1), model%c0(1, 1), model%ff(1, 1), model%v(1, 1))
        allocate(model%gg(1, 1), model%w(1, 1))
        model%m0 = 0.0_dp
        model%c0 = 0.0_dp
        model%ff = 0.0_dp
        model%v = observation_variance
        model%gg = 1.0_dp
        model%w = 0.0_dp
    end subroutine build_iid_variance

end module test_mle_types

program test_mle
    use dlm, only : dp, dlm_success, dlm_mle, dlm_mle_result
    use test_mle_types, only : iid_variance_builder
    implicit none

    type(iid_variance_builder) :: builder
    type(dlm_mle_result) :: fit
    real(dp), parameter :: expected_log_variance = log(2.0_dp)
    real(dp), parameter :: expected_objective = 0.5_dp * (5.0_dp * log(2.0_dp) + 5.0_dp)
    real(dp) :: y(5)
    integer :: info

    y = [-2.0_dp, -1.0_dp, 0.0_dp, 1.0_dp, 2.0_dp]
    call dlm_mle(y, [0.0_dp], builder, fit, info, lower=[-3.0_dp], upper=[3.0_dp], &
                 max_iterations=200, max_evaluations=1000, tolerance=1.0e-10_dp)
    if (info /= dlm_success) error stop "dlmMLE adapter failed"
    if (.not. allocated(fit%par)) error stop "dlmMLE parameter result missing"
    if (size(fit%par) /= 1) error stop "dlmMLE parameter shape failed"
    if (.not. fit%converged) error stop "dlmMLE optimizer did not converge"
    if (abs(fit%par(1) - expected_log_variance) > 2.0e-4_dp) error stop "dlmMLE variance estimate failed"
    if (abs(fit%value - expected_objective) > 2.0e-7_dp) error stop "dlmMLE objective value failed"
    if (fit%function_count < 1) error stop "dlmMLE function count failed"

    print *, "test_mle: PASS"
end program test_mle
