! Computational translation of R package changepoint 2.3.
! Upstream license declaration: GPL (unversioned).
! See NOTICE.md and UPSTREAM.md for authorship and provenance.
program test_changepoint
use r_kinds, only : dp
use changepoint, only : cp_cost_mean_normal, cp_cost_meanvar_normal, cp_cost_exponential, cp_cost_gamma, cp_cost_poisson
use changepoint, only : cp_pelt, cp_binseg, cp_segneigh, cp_amoc, cp_amoc_css, cp_amoc_cusum
use changepoint, only : cp_binseg_css, cp_binseg_cusum, cp_crops
use changepoint, only : cp_regression_amoc, cp_regression_pelt
use changepoint, only : cp_segment_means, cp_segment_variances_mle, cp_segment_regression_fits
use changepoint, only : cp_segment_trend_fits
use changepoint, only : cp_penalty_value, cp_amoc_asymptotic_value, cp_decision
use changepoint, only : changepoint_result, binseg_result, segneigh_result, amoc_result, crops_solution
implicit none

integer :: failures
failures = 0
call test_pelt_mean(failures)
call test_amoc(failures)
call test_binseg_segneigh(failures)
call test_pelt_vs_exact_dp(failures)
call test_distributions(failures)
call test_nonparametric(failures)
call test_crops(failures)
call test_regression(failures)
call test_fits_penalties(failures)
if (failures /= 0) then
    error stop 1
end if
print '(a)', 'all changepoint tests passed'

contains

subroutine test_pelt_mean(failures)
    integer, intent(inout) :: failures
    real(dp) :: x(12)
    type(changepoint_result) :: r
    x = [0.0_dp, 0.1_dp, -0.1_dp, 0.05_dp, 5.0_dp, 5.1_dp, 4.9_dp, 5.05_dp, &
        -2.0_dp, -2.1_dp, -1.9_dp, -2.05_dp]
    call cp_pelt(x, cp_cost_mean_normal, 1.0_dp, 2, r)
    call expect(r%status == 0, 'PELT status', failures)
    call expect(r%ncpts == 2, 'PELT number of changes', failures)
    if (r%ncpts == 2) then
        call expect(all(r%cpts == [4, 8]), 'PELT locations', failures)
    end if
end subroutine test_pelt_mean

subroutine test_amoc(failures)
    integer, intent(inout) :: failures
    real(dp) :: x(10)
    type(amoc_result) :: r
    x = [1.0_dp, 1.1_dp, 0.9_dp, 1.0_dp, 1.05_dp, 4.0_dp, 4.1_dp, 3.9_dp, 4.0_dp, 4.05_dp]
    call cp_amoc(x, cp_cost_mean_normal, 0.5_dp, 2, r)
    call expect(r%status == 0 .and. r%changed .and. r%cpt == 5, 'AMOC normal mean', failures)
end subroutine test_amoc

subroutine test_binseg_segneigh(failures)
    integer, intent(inout) :: failures
    real(dp) :: x(12)
    type(binseg_result) :: b
    type(segneigh_result) :: s
    x = [0.0_dp, 0.1_dp, -0.1_dp, 0.05_dp, 5.0_dp, 5.1_dp, 4.9_dp, 5.05_dp, &
        -2.0_dp, -2.1_dp, -1.9_dp, -2.05_dp]
    call cp_binseg(x, cp_cost_mean_normal, 1.0_dp, 2, 4, b)
    call expect(b%status == 0, 'BinSeg status', failures)
    call expect(b%ncpts >= 2, 'BinSeg count', failures)
    if (b%ncpts >= 2) call expect(all(b%cpts(1:2) == [4, 8]), 'BinSeg first locations', failures)
    call cp_segneigh(x, cp_cost_mean_normal, 1.0_dp, 2, 4, s)
    call expect(s%status == 0 .and. s%ncpts == 2, 'SegNeigh count', failures)
    if (s%ncpts == 2) call expect(all(s%cpts == [4, 8]), 'SegNeigh locations', failures)
end subroutine test_binseg_segneigh

subroutine test_pelt_vs_exact_dp(failures)
    integer, intent(inout) :: failures
    real(dp) :: x(15)
    type(changepoint_result) :: p
    type(segneigh_result) :: s

    x = [0.2_dp, -0.4_dp, 0.1_dp, 0.3_dp, 2.0_dp, 2.4_dp, 1.8_dp, 2.1_dp, &
        -1.2_dp, -0.8_dp, -1.0_dp, -1.3_dp, 0.5_dp, 0.4_dp, 0.7_dp]
    call cp_pelt(x, cp_cost_mean_normal, 1.7_dp, 2, p)
    call cp_segneigh(x, cp_cost_mean_normal, 1.7_dp, 2, 6, s)
    call expect(p%status == 0 .and. s%status == 0, 'PELT/exact DP status', failures)
    call expect(p%ncpts == s%ncpts, 'PELT/exact DP count', failures)
    if (p%ncpts == s%ncpts .and. p%ncpts > 0) then
        call expect(all(p%cpts == s%cpts), 'PELT/exact DP locations', failures)
    end if
end subroutine test_pelt_vs_exact_dp

subroutine test_distributions(failures)
    integer, intent(inout) :: failures
    real(dp) :: e(10), g(10), p(10)
    type(changepoint_result) :: r
    e = [1.0_dp, 1.2_dp, 0.8_dp, 1.1_dp, 1.0_dp, 5.0_dp, 5.2_dp, 4.8_dp, 5.1_dp, 4.9_dp]
    g = e
    p = [1.0_dp, 2.0_dp, 1.0_dp, 2.0_dp, 1.0_dp, 8.0_dp, 9.0_dp, 8.0_dp, 9.0_dp, 8.0_dp]
    call cp_pelt(e, cp_cost_exponential, 1.0_dp, 2, r)
    call expect(r%status == 0 .and. r%ncpts >= 1, 'Exponential PELT', failures)
    call cp_pelt(g, cp_cost_gamma, 1.0_dp, 2, r, shape=2.0_dp)
    call expect(r%status == 0 .and. r%ncpts >= 1, 'Gamma PELT', failures)
    call cp_pelt(p, cp_cost_poisson, 1.0_dp, 2, r)
    call expect(r%status == 0 .and. r%ncpts >= 1, 'Poisson PELT', failures)
    call cp_pelt(e, cp_cost_meanvar_normal, 2.0_dp, 2, r)
    call expect(r%status == 0, 'Mean-variance normal PELT', failures)
    e = [0.0_dp, 0.0_dp, 0.2_dp, 0.1_dp, 0.0_dp, 3.0_dp, 3.2_dp, 2.8_dp, 3.1_dp, 2.9_dp]
    call cp_pelt(e, cp_cost_exponential, 1.0_dp, 2, r)
    call expect(r%status == 0, 'Exponential accepts zero data', failures)
    p = [1.0_dp, 2.0_dp, 1.5_dp, 2.0_dp, 1.0_dp, 8.0_dp, 9.0_dp, 8.0_dp, 9.0_dp, 8.0_dp]
    call cp_pelt(p, cp_cost_poisson, 1.0_dp, 2, r)
    call expect(r%status /= 0, 'Poisson rejects noninteger data', failures)
end subroutine test_distributions

subroutine test_nonparametric(failures)
    integer, intent(inout) :: failures
    real(dp) :: x(10), v(10)
    type(amoc_result) :: a
    type(binseg_result) :: b
    x = [0.0_dp, 0.1_dp, -0.1_dp, 0.05_dp, 0.0_dp, 4.0_dp, 4.1_dp, 3.9_dp, 4.0_dp, 4.1_dp]
    v = [0.1_dp, -0.1_dp, 0.05_dp, -0.05_dp, 0.1_dp, 3.0_dp, -3.0_dp, 2.5_dp, -2.5_dp, 3.2_dp]
    call cp_amoc_cusum(x, 0.2_dp, 2, a)
    call expect(a%status == 0 .and. a%changed .and. a%cpt == 5, 'CUSUM AMOC', failures)
    call cp_amoc_css(v, 0.5_dp, 2, a)
    call expect(a%status == 0 .and. a%changed, 'CSS AMOC', failures)
    call cp_binseg_cusum(x, 0.1_dp, 2, 3, b)
    call expect(b%status == 0 .and. b%ncpts >= 1, 'CUSUM BinSeg', failures)
    call cp_binseg_css(v, 0.5_dp, 2, 3, b)
    call expect(b%status == 0 .and. b%ncpts >= 1, 'CSS BinSeg', failures)
end subroutine test_nonparametric

subroutine test_crops(failures)
    integer, intent(inout) :: failures
    real(dp) :: x(12)
    type(crops_solution), allocatable :: sol(:)
    integer :: status
    x = [0.0_dp, 0.1_dp, -0.1_dp, 0.05_dp, 5.0_dp, 5.1_dp, 4.9_dp, 5.05_dp, &
        -2.0_dp, -2.1_dp, -1.9_dp, -2.05_dp]
    allocate(sol(0))
    call cp_crops(x, cp_cost_mean_normal, 0.2_dp, 200.0_dp, 2, sol, status)
    call expect(status == 0 .and. size(sol) >= 2, 'CROPS solutions', failures)
    if (size(sol) >= 2) call expect(sol(1)%ncpts >= sol(size(sol))%ncpts, 'CROPS monotone counts', failures)
end subroutine test_crops

subroutine test_regression(failures)
    integer, intent(inout) :: failures
    integer, parameter :: n = 12
    real(dp) :: y(n), x(n, 2)
    type(amoc_result) :: a
    type(changepoint_result) :: r
    integer :: i
    do i = 1, n
        x(i, 1) = 1.0_dp
        x(i, 2) = real(i, dp)
        if (i <= 6) then
            y(i) = 1.0_dp + 0.5_dp * x(i, 2)
        else
            y(i) = -2.0_dp + 1.5_dp * x(i, 2)
        end if
    end do
    call cp_regression_amoc(y, x, 1.0_dp, 3, -1.0_dp, a)
    call expect(a%status == 0 .and. a%changed .and. a%cpt == 6, 'Regression AMOC', failures)
    call cp_regression_pelt(y, x, 1.0_dp, 3, -1.0_dp, r)
    call expect(r%status == 0 .and. r%ncpts == 1, 'Regression PELT count', failures)
    if (r%ncpts == 1) call expect(r%cpts(1) == 6, 'Regression PELT location', failures)
end subroutine test_regression

subroutine test_fits_penalties(failures)
    integer, intent(inout) :: failures
    real(dp) :: x(8), value
    real(dp), allocatable :: means(:), vars(:), beta(:, :), sig2(:), theta0(:), theta1(:)
    real(dp) :: y(8), design(8, 1)
    integer, allocatable :: decided(:)
    integer :: status, i
    x = [1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 3.0_dp, 3.0_dp, 3.0_dp, 3.0_dp]
    call cp_segment_means(x, [4], means, status)
    call expect(status == 0 .and. maxval(abs(means - [1.0_dp, 3.0_dp])) < 1.0e-12_dp, 'Segment means', failures)
    call cp_segment_variances_mle(x, [4], vars, status)
    call expect(status == 0 .and. maxval(abs(vars)) < 1.0e-12_dp, 'Segment variances', failures)
    x = [5.0_dp, 8.0_dp, 11.0_dp, 14.0_dp, 1.0_dp, 3.0_dp, 5.0_dp, 7.0_dp]
    call cp_segment_trend_fits(x, [4], theta0, theta1, status)
    call expect(status == 0 .and. maxval(abs(theta0 - [2.0_dp, -1.0_dp])) < 1.0e-12_dp, &
        'Segment trend starts', failures)
    call expect(maxval(abs(theta1 - [14.0_dp, 7.0_dp])) < 1.0e-12_dp, 'Segment trend ends', failures)
    design(:, 1) = 1.0_dp
    do i = 1, 8
        y(i) = merge(2.0_dp, 5.0_dp, i <= 4)
    end do
    call cp_segment_regression_fits(y, design, [4], beta, sig2, status)
    call expect(status == 0 .and. abs(beta(1, 1) - 2.0_dp) < 1.0e-10_dp .and. &
        abs(beta(2, 1) - 5.0_dp) < 1.0e-10_dp, 'Segment regression fit', failures)
    call cp_penalty_value('BIC', 100, 1, value, status)
    call expect(status == 0 .and. abs(value - 2.0_dp * log(100.0_dp)) < 1.0e-12_dp, 'BIC penalty', failures)
    call cp_amoc_asymptotic_value(cp_cost_mean_normal, 100, 20.0_dp, 5.0_dp, value, status)
    call expect(status == 0 .and. value >= 0.0_dp .and. value <= 1.0_dp, 'AMOC asymptotic value', failures)
    call cp_decision([5, 7], [10.0_dp, 3.0_dp], 4.0_dp, 10, decided, status, [2.0_dp, 1.0_dp])
    call expect(status == 0 .and. all(decided == [5, 10]), 'Decision helper', failures)
end subroutine test_fits_penalties

subroutine expect(condition, label, failures)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: label
    integer, intent(inout) :: failures
    if (.not. condition) then
        failures = failures + 1
        print '(a)', 'FAIL: ' // trim(label)
    end if
end subroutine expect

end program test_changepoint
