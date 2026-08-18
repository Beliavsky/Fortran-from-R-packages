! SPDX-License-Identifier: GPL-2.0-or-later
! Modern Fortran translation of computational code from R package good 1.0.2.

program test_glm
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
    use good, only : dp, good_glm_fit, good_prediction, good_glm_summary, glm_good, predict_good, summary_good
    implicit none

    integer, allocatable :: y(:)
    real(dp), allocatable :: x(:, :), x2(:, :)
    type(good_glm_fit) :: fit, fit2, fit3
    type(good_prediction) :: pred
    type(good_glm_summary) :: smry
    integer :: n, i, k
    real(dp) :: sample_mean

    ! Dataset from the original R package example/tests.
    n = 46 + 76 + 24 + 9 + 1
    allocate(y(n), x(n, 1))
    k = 0
    call put_repeated(0, 46)
    call put_repeated(1, 76)
    call put_repeated(2, 24)
    call put_repeated(3, 9)
    call put_repeated(4, 1)
    x = 1.0_dp

    call glm_good(y, x, 'log', fit, max_iter=300)
    if (.not. fit%converged) error stop 'intercept-only glm_good did not converge'
    if (size(fit%coefficients) /= 2) error stop 'wrong coefficient count'
    if (fit%coefficients(2) >= 0.0_dp) error stop 'log-link z is not below one'
    sample_mean = sum(real(y, dp)) / real(n, dp)
    if (abs(sum(fit%fitted_values) / real(n, dp) - sample_mean) > 5.0e-4_dp) &
        error stop 'intercept-only fitted mean mismatch'
    if (any(ieee_is_nan(fit%vcov))) error stop 'vcov contains NaN'
    if (abs(fit%coefficients(1) + 4.77857735_dp) > 2.0e-5_dp) error stop 's estimate reference mismatch'
    if (abs(fit%coefficients(2) + 2.86615814_dp) > 2.0e-5_dp) error stop 'beta estimate reference mismatch'
    if (abs(fit%loglik + 187.5738649834_dp) > 2.0e-7_dp) error stop 'loglik reference mismatch'

    call predict_good(fit, x, pred, with_se=.true.)
    if (maxval(abs(pred%fit - fit%fitted_values)) > 1.0e-9_dp) error stop 'prediction mismatch'
    if (any(pred%se_fit <= 0.0_dp)) error stop 'prediction SE invalid'

    call summary_good(fit, y, x, smry)
    if (abs(smry%aic - (4.0_dp - 2.0_dp * fit%loglik)) > 1.0e-10_dp) error stop 'AIC mismatch'
    if (size(smry%lrt) /= 2) error stop 'intercept summary LRT count mismatch'

    ! Add a simple covariate and exercise logit regression/LRT path.
    allocate(x2(n, 2))
    x2(:, 1) = 1.0_dp
    do i = 1, n
        x2(i, 2) = (real(i, dp) - 0.5_dp * real(n + 1, dp)) / real(n, dp)
    end do
    call glm_good(y, x2, 'logit', fit2, max_iter=250)
    if (.not. fit2%converged) error stop 'logit glm_good did not converge'
    call summary_good(fit2, y, x2, smry)
    if (size(smry%lrt) /= 1 .or. smry%df(1) /= 1) error stop 'covariate LRT summary failed'

    call glm_good(y, x, 'identity', fit3, max_iter=300)
    if (.not. fit3%converged) error stop 'identity glm_good did not converge'
    if (fit3%coefficients(2) <= 0.0_dp .or. fit3%coefficients(2) >= 1.0_dp) &
        error stop 'identity link fitted z outside (0,1)'

    print '(a)', 'test_glm: PASS'

contains

    subroutine put_repeated(value, count)
        integer, intent(in) :: value, count
        integer :: j
        do j = 1, count
            k = k + 1
            y(k) = value
        end do
    end subroutine put_repeated

end program test_glm
