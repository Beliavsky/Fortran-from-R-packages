! SPDX-License-Identifier: GPL-2.0-or-later
module mle_variance_types
    use dlm, only : dp, dlm_model, dlm_success, dlm_invalid_argument, dlm_model_builder
    implicit none
    private

    type, extends(dlm_model_builder), public :: variance_builder
        real(dp) :: scale = 1.0_dp
    contains
        procedure, pass :: build => build_variance_model
    end type variance_builder

contains

    subroutine build_variance_model(self, parameters, model, info)
        class(variance_builder), intent(in) :: self !! Builder holding a fixed positive observation-variance scale.
        real(dp), intent(in) :: parameters(:) !! One-element vector containing the logarithm of scaled observation variance.
        type(dlm_model), intent(out) :: model !! IID zero-mean Gaussian model represented as a one-state DLM.
        integer, intent(out) :: info !! Zero for valid builder input or `dlm_invalid_argument` otherwise.

        info = dlm_success
        if (size(parameters) /= 1 .or. self%scale <= 0.0_dp) then
            info = dlm_invalid_argument
            return
        end if
        allocate(model%m0(1), model%c0(1, 1), model%ff(1, 1), model%v(1, 1))
        allocate(model%gg(1, 1), model%w(1, 1))
        model%m0 = 0.0_dp
        model%c0 = 0.0_dp
        model%ff = 0.0_dp
        model%v = self%scale * exp(parameters(1))
        model%gg = 1.0_dp
        model%w = 0.0_dp
    end subroutine build_variance_model

end module mle_variance_types

program mle_variance
    use dlm, only : dp, dlm_success, dlm_mle, dlm_mle_result
    use mle_variance_types, only : variance_builder
    implicit none

    type(variance_builder) :: builder
    type(dlm_mle_result) :: fit
    real(dp) :: y(5)
    integer :: info

    y = [-2.0_dp, -1.0_dp, 0.0_dp, 1.0_dp, 2.0_dp]
    call dlm_mle(y, [0.0_dp], builder, fit, info, lower=[-3.0_dp], upper=[3.0_dp])
    if (info /= dlm_success .or. .not. fit%converged) error stop "dlmMLE example failed"
    print '(a,f10.5)', 'fitted log variance = ', fit%par(1)
end program mle_variance
