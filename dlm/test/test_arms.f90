! SPDX-License-Identifier: GPL-2.0-or-later
module test_arms_types
    use dlm, only : dp, dlm_indicator, dlm_log_density
    implicit none
    private

    type, extends(dlm_indicator), public :: box_indicator
        real(dp) :: half_width = 4.0_dp
    contains
        procedure, pass :: evaluate => box_evaluate
    end type box_indicator

    type, extends(dlm_log_density), public :: standard_normal_density
        real(dp) :: precision = 1.0_dp
    contains
        procedure, pass :: evaluate => normal_log_density
    end type standard_normal_density

    type, extends(dlm_log_density), public :: bimodal_normal_density
        real(dp) :: mode = 2.0_dp
        real(dp) :: precision = 4.0_dp
    contains
        procedure, pass :: evaluate => bimodal_log_density
    end type bimodal_normal_density

contains

    pure logical function box_evaluate(self, x) result(inside)
        class(box_indicator), intent(in) :: self !! Hypercube support object containing the common component half-width.
        real(dp), intent(in) :: x(:) !! Candidate point tested for membership in the closed support hypercube.

        inside = all(abs(x) <= self%half_width)
    end function box_evaluate

    pure real(dp) function normal_log_density(self, x) result(log_density)
        class(standard_normal_density), intent(in) :: self !! Independent Gaussian target object containing common precision.
        real(dp), intent(in) :: x(:) !! Point whose unnormalized independent Gaussian log density is evaluated.

        log_density = -0.5_dp * self%precision * sum(x**2)
    end function normal_log_density

    pure real(dp) function bimodal_log_density(self, x) result(log_density)
        class(bimodal_normal_density), intent(in) :: self !! Symmetric two-mode Gaussian-mixture target parameters.
        real(dp), intent(in) :: x(:) !! Point whose unnormalized symmetric Gaussian-mixture log density is evaluated.

        real(dp) :: left_log, maximum_log, right_log

        left_log = -0.5_dp * self%precision * (x(1) + self%mode)**2
        right_log = -0.5_dp * self%precision * (x(1) - self%mode)**2
        maximum_log = max(left_log, right_log)
        log_density = maximum_log + log(exp(left_log - maximum_log) + exp(right_log - maximum_log))
    end function bimodal_log_density

end module test_arms_types

program test_arms
    use dlm, only : dp, dlm_success, arms
    use test_arms_types, only : box_indicator, standard_normal_density, bimodal_normal_density
    implicit none

    type(box_indicator) :: support
    type(standard_normal_density) :: density
    type(bimodal_normal_density) :: bimodal_density
    real(dp), allocatable :: draws_a(:, :), draws_b(:, :), draws_mv(:, :), draws_bimodal(:, :)
    integer :: info

    call arms([0.0_dp], density, support, 2000, draws_a, info, seed=24681357)
    if (info /= dlm_success) error stop "univariate ARMS sampling failed"
    if (any(shape(draws_a) /= [2000, 1])) error stop "univariate ARMS shape failed"
    if (any(abs(draws_a(:, 1)) > support%half_width + 1.0e-8_dp)) error stop "univariate ARMS support failed"
    if (abs(sum(draws_a(:, 1)) / real(size(draws_a, 1), dp)) > 0.12_dp) error stop "univariate ARMS mean failed"
    if (sum(draws_a(:, 1)**2) / real(size(draws_a, 1), dp) < 0.70_dp) error stop "univariate ARMS variance low"
    if (sum(draws_a(:, 1)**2) / real(size(draws_a, 1), dp) > 1.30_dp) error stop "univariate ARMS variance high"

    call arms([0.0_dp], density, support, 2000, draws_b, info, seed=24681357)
    if (info /= dlm_success) error stop "repeat ARMS sampling failed"
    if (maxval(abs(draws_a - draws_b)) > 0.0_dp) error stop "ARMS deterministic seeding failed"

    call arms([0.0_dp], bimodal_density, support, 2000, draws_bimodal, info, seed=86420)
    if (info /= dlm_success) error stop "non-log-concave ARMS sampling failed"
    if (any(abs(draws_bimodal(:, 1)) > support%half_width + 1.0e-8_dp)) then
        error stop "non-log-concave ARMS support failed"
    end if
    if (count(draws_bimodal(:, 1) < -0.5_dp) < 250) error stop "non-log-concave ARMS left mode failed"
    if (count(draws_bimodal(:, 1) > 0.5_dp) < 250) error stop "non-log-concave ARMS right mode failed"

    call arms([0.0_dp, 0.0_dp], density, support, 400, draws_mv, info, seed=97531)
    if (info /= dlm_success) error stop "multivariate hit-and-run ARMS failed"
    if (any(shape(draws_mv) /= [400, 2])) error stop "multivariate ARMS shape failed"
    if (any(abs(draws_mv) > support%half_width + 1.0e-8_dp)) error stop "multivariate ARMS support failed"
    if (maxval(abs(sum(draws_mv, dim=1) / real(size(draws_mv, 1), dp))) > 0.35_dp) then
        error stop "multivariate ARMS mean failed"
    end if

    print *, "test_arms: PASS"
end program test_arms
