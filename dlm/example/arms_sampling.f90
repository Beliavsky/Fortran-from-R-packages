! SPDX-License-Identifier: GPL-2.0-or-later
module arms_sampling_types
    use dlm, only : dp, dlm_indicator, dlm_log_density
    implicit none
    private

    type, extends(dlm_indicator), public :: interval_support
        real(dp) :: half_width = 4.0_dp
    contains
        procedure, pass :: evaluate => interval_evaluate
    end type interval_support

    type, extends(dlm_log_density), public :: gaussian_density
        real(dp) :: precision = 1.0_dp
    contains
        procedure, pass :: evaluate => gaussian_log_density
    end type gaussian_density

contains

    pure logical function interval_evaluate(self, x) result(inside)
        class(interval_support), intent(in) :: self !! Symmetric interval support object containing its half-width.
        real(dp), intent(in) :: x(:) !! Candidate point tested for membership in the interval support.

        inside = size(x) == 1 .and. abs(x(1)) <= self%half_width
    end function interval_evaluate

    pure real(dp) function gaussian_log_density(self, x) result(log_density)
        class(gaussian_density), intent(in) :: self !! Gaussian target object containing its precision parameter.
        real(dp), intent(in) :: x(:) !! Point whose unnormalized Gaussian log density is evaluated.

        log_density = -0.5_dp * self%precision * sum(x**2)
    end function gaussian_log_density

end module arms_sampling_types

program arms_sampling
    use dlm, only : dp, dlm_success, arms
    use arms_sampling_types, only : interval_support, gaussian_density
    implicit none

    type(interval_support) :: support
    type(gaussian_density) :: density
    real(dp), allocatable :: draws(:, :)
    integer :: info

    call arms([0.0_dp], density, support, 1000, draws, info, seed=12345)
    if (info /= dlm_success) error stop "ARMS example failed"
    print '(a,f10.5)', 'sample mean = ', sum(draws(:, 1)) / real(size(draws, 1), dp)
end program arms_sampling
