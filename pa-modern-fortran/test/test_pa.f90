! SPDX-License-Identifier: GPL-2.0-only
! Copyright (C) 2010-2023 Yang Lu and David Kane
! Copyright (C) 2026 Modern Fortran translation contributors
! This program is free software under GNU GPL version 2 only.
program test_pa
  use pa
  implicit none

  call test_utils_and_exposure()
  call test_brinson()
  call test_regression()
  print '(a)', 'All performance-attribution tests passed.'

contains

  subroutine test_utils_and_exposure()
    real(dp) :: v(10), wp(10), wb(10), numeric2(6,1)
    integer :: c(10), periods(10), cats(6,2), levels(2), status
    real(dp), allocatable :: ranks(:), x(:, :)
    integer, allocatable :: quint(:), gs(:), ge(:)
    type(exposure_result) :: e
    type(exposure_multi_result) :: em
    integer :: i

    v = [(real(i,dp), i=1,10)]
    wp = 0.1_dp
    wb = [0.05_dp,0.15_dp,0.05_dp,0.15_dp,0.05_dp,0.15_dp,0.05_dp,0.15_dp,0.05_dp,0.15_dp]
    call average_ranks([1.0_dp,2.0_dp,2.0_dp,4.0_dp], ranks)
    call assert_close_vec(ranks, [1.0_dp,2.5_dp,2.5_dp,4.0_dp], 1.0e-14_dp, 'average ranks')
    call quintile_groups(v, quint)
    call assert_int_vec(quint, [1,1,2,2,3,3,4,4,5,5], 'quintile groups')

    call continuous_exposure(v, wp, wb, e, status)
    call assert_true(status == 0, 'continuous exposure status')
    call assert_close_vec(e%portfolio, [0.2_dp,0.2_dp,0.2_dp,0.2_dp,0.2_dp], 1.0e-14_dp, 'continuous portfolio')
    call assert_close_vec(e%benchmark, [0.2_dp,0.2_dp,0.2_dp,0.2_dp,0.2_dp], 1.0e-14_dp, 'continuous benchmark')

    c = [1,1,2,2,3,3,1,2,3,3]
    call categorical_exposure(c, wp, wb, e, status)
    call assert_true(status == 0, 'categorical exposure status')
    call assert_int_vec(e%group, [1,2,3], 'categorical group')
    call assert_close_vec(e%portfolio, [0.3_dp,0.3_dp,0.4_dp], 1.0e-14_dp, 'categorical portfolio')

    periods = [1,1,1,1,1,2,2,2,2,2]
    call continuous_exposure_multi(periods, v, wp, wb, em, status)
    call assert_true(status == 0, 'continuous multi status')
    call assert_close_vec(sum(em%portfolio,dim=1), [0.5_dp,0.5_dp], 1.0e-14_dp, 'multi exposure sums')

    cats(:,1) = [1,2,3,1,2,3]
    cats(:,2) = [1,1,2,2,3,3]
    levels = [3,3]
    numeric2(:,1) = [1.0_dp,2.0_dp,3.0_dp,4.0_dp,5.0_dp,6.0_dp]
    call build_design_matrix(cats, levels, numeric2, x, gs, ge, status)
    call assert_true(status == 0, 'design status')
    call assert_true(all(shape(x) == [6,6]), 'design shape')
    call assert_int_vec(gs, [1,4,6], 'design group starts')
    call assert_int_vec(ge, [3,5,6], 'design group ends')
    call assert_close_vec(x(1,:), [1.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,1.0_dp], 1.0e-14_dp, 'design first row')
    call assert_close_vec(x(6,:), [0.0_dp,0.0_dp,1.0_dp,0.0_dp,1.0_dp,6.0_dp], 1.0e-14_dp, 'design last row')
  end subroutine test_utils_and_exposure

  subroutine test_brinson()
    integer :: category(12), period(12)
    real(dp) :: wb(12), wp(12), r(12)
    type(brinson_period_result) :: one
    type(brinson_multi_result) :: multi
    type(attribution_summary) :: summary

    category = [1,1,2,2,3,3, 1,1,2,2,3,3]
    period = [1,1,1,1,1,1, 2,2,2,2,2,2]
    wb = [0.15_dp,0.15_dp,0.20_dp,0.20_dp,0.15_dp,0.15_dp, &
          0.10_dp,0.20_dp,0.15_dp,0.25_dp,0.10_dp,0.20_dp]
    wp = [0.20_dp,0.10_dp,0.25_dp,0.15_dp,0.10_dp,0.20_dp, &
          0.15_dp,0.20_dp,0.20_dp,0.20_dp,0.15_dp,0.10_dp]
    r = [0.02_dp,0.01_dp,-0.01_dp,0.03_dp,0.04_dp,0.00_dp, &
         0.01_dp,0.025_dp,0.02_dp,-0.005_dp,0.03_dp,0.015_dp]

    call fit_brinson_period(category(1:6), wb(1:6), wp(1:6), r(1:6), one)
    call assert_true(one%status == 0, 'Brinson single status')
    call assert_close_vec(one%weight_portfolio, [0.3_dp,0.4_dp,0.3_dp], 1.0e-14_dp, 'Brinson weights')
    call assert_close_vec(one%return_portfolio, [0.0166666666666667_dp,0.005_dp,0.0133333333333333_dp], &
                          1.0e-13_dp, 'Brinson portfolio category returns')
    call assert_close_vec(one%return_benchmark, [0.015_dp,0.01_dp,0.02_dp], 1.0e-14_dp, &
                          'Brinson benchmark category returns')
    call assert_close_vec(one%q, [0.011_dp,0.011_dp,0.0145_dp,0.0145_dp], 1.0e-14_dp, 'Brinson q')
    call assert_close_vec(one%aggregate, [0.0_dp,-0.0035_dp,0.0_dp,-0.0035_dp], 1.0e-14_dp, 'Brinson aggregate')
    call assert_close(sum(one%category_effect(:,1)), one%aggregate(1), 1.0e-14_dp, 'category allocation sum')
    call assert_close(sum(one%category_effect(:,2)), one%aggregate(2), 1.0e-14_dp, 'category selection sum')
    call assert_close(sum(one%category_effect(:,3)), one%aggregate(3), 1.0e-14_dp, 'category interaction sum')

    call fit_brinson_multi(period, category, wb, wp, r, multi)
    call assert_true(multi%status == 0, 'Brinson multi status')
    call assert_close_vec(multi%raw(:,2), [0.0_dp,0.00202142857142857_dp,-0.000271428571428571_dp,0.00175_dp], &
                          1.0e-13_dp, 'Brinson period 2')

    call summarize_brinson_multi(multi, 'arithmetic', summary)
    call assert_close_vec(summary%aggregate, [0.0_dp,-0.00147857142857143_dp,-0.000271428571428571_dp,-0.00175_dp], &
                          1.0e-13_dp, 'Brinson arithmetic')
    call summarize_brinson_multi(multi, 'geometric', summary)
    call assert_close_vec(summary%aggregate, [0.0_dp,-0.001504460714285714_dp,-0.000274414285714286_dp,-0.001778875_dp], &
                          2.0e-13_dp, 'Brinson geometric')
    call summarize_brinson_multi(multi, 'linking', summary)
    call assert_close_vec(summary%linking_coefficient, [1.0148117279_dp,1.0131234579_dp], 2.0e-9_dp, &
                          'Brinson linking coefficients')
    call assert_close(summary%aggregate(4), -0.001778875_dp, 2.0e-13_dp, 'Brinson linked active')
    call assert_close(sum(summary%aggregate(1:3)), summary%aggregate(4), 2.0e-13_dp, 'Brinson linked additivity')
  end subroutine test_brinson

  subroutine test_regression()
    integer, parameter :: n = 9
    integer :: categorical(n,1), levels(1), status, i
    integer :: period(2*n)
    real(dp) :: numeric(n,2), wb(n), wp(n), beta1(5), beta2(5), y1(n), y2(n)
    real(dp), allocatable :: x(:, :), xmulti(:, :), y(:), wbm(:), wpm(:), one_summary(:)
    integer, allocatable :: gs(:), ge(:)
    type(regression_period_result) :: one
    type(regression_multi_result) :: multi
    type(regression_summary) :: summary

    categorical(:,1) = [1,2,3,1,2,3,1,2,3]
    levels = [3]
    numeric(:,1) = [-1.0_dp,-0.5_dp,0.0_dp,0.5_dp,1.0_dp,-0.8_dp,0.2_dp,0.7_dp,-0.3_dp]
    numeric(:,2) = [0.4_dp,-0.2_dp,0.1_dp,0.6_dp,-0.5_dp,0.3_dp,-0.7_dp,0.2_dp,0.8_dp]
    wb = [0.10_dp,0.12_dp,0.08_dp,0.13_dp,0.09_dp,0.11_dp,0.12_dp,0.15_dp,0.10_dp]
    wp = [0.12_dp,0.10_dp,0.08_dp,0.10_dp,0.12_dp,0.10_dp,0.14_dp,0.12_dp,0.12_dp]
    beta1 = [0.01_dp,0.02_dp,-0.005_dp,0.003_dp,-0.002_dp]
    beta2 = [0.008_dp,0.018_dp,-0.002_dp,-0.001_dp,0.004_dp]

    call build_design_matrix(categorical, levels, numeric, x, gs, ge, status)
    call assert_true(status == 0, 'regression design status')
    y1 = matmul(x,beta1)
    y2 = matmul(x,beta2)
    call fit_regression_period(y1, x, wb, wp, one)
    call assert_true(one%status == 0, 'regression period status')
    call assert_true(one%rank == 5, 'regression rank')
    call assert_close_vec(one%coefficients, beta1, 2.0e-13_dp, 'regression coefficients')
    call assert_close_vec(one%active_exposure, [0.01_dp,-0.02_dp,0.01_dp,-0.01_dp,-0.028_dp], &
                          1.0e-14_dp, 'regression active exposure')
    call assert_close_vec(one%contribution, [0.0001_dp,-0.0004_dp,-0.00005_dp,-0.00003_dp,0.000056_dp], &
                          2.0e-14_dp, 'regression contribution')
    call summarize_regression_period(one, gs, ge, one_summary, status)
    call assert_close_vec(one_summary, [-0.00035_dp,-0.00003_dp,0.000056_dp,0.0_dp,0.008712_dp,0.009036_dp,-0.000324_dp], &
                          2.0e-13_dp, 'regression return summary')

    allocate(xmulti(2*n,5), y(2*n), wbm(2*n), wpm(2*n))
    xmulti(1:n,:) = x
    xmulti(n+1:2*n,:) = x
    y(1:n) = y1
    y(n+1:2*n) = y2
    wbm(1:n) = wb
    wbm(n+1:2*n) = wb
    wpm(1:n) = wp
    wpm(n+1:2*n) = wp
    period = [(1,i=1,n),(2,i=1,n)]
    call fit_regression_multi(period, y, xmulti, wbm, wpm, multi)
    call assert_true(multi%status == 0, 'regression multi status')
    call assert_close_vec(multi%coefficients(:,1), beta1, 2.0e-13_dp, 'regression multi beta1')
    call assert_close_vec(multi%coefficients(:,2), beta2, 2.0e-13_dp, 'regression multi beta2')

    call summarize_regression_multi(multi, gs, ge, 'arithmetic', summary)
    call assert_close_vec(summary%aggregate, [-0.00065_dp,-0.00002_dp,-0.000056_dp,0.0_dp,0.017468_dp,0.018194_dp,-0.000726_dp], &
                          3.0e-13_dp, 'regression arithmetic')
    call summarize_regression_multi(multi, gs, ge, 'geometric', summary)
    call assert_close_vec(summary%aggregate, [-0.000649895_dp,-0.0000200003_dp,-0.000056006272_dp,0.0_dp, &
                          0.017544282272_dp,0.018276751688_dp,-0.000732469416_dp], 5.0e-13_dp, 'regression geometric')
    call summarize_regression_multi(multi, gs, ge, 'linking', summary)
    call assert_close_vec(summary%linking_coefficient, [1.008911566_dp,1.008910619_dp], 2.0e-9_dp, &
                          'regression linking coefficient')
    call assert_close(summary%aggregate(size(summary%aggregate)), -0.000732469416_dp, 5.0e-13_dp, &
                      'regression linked active')
  end subroutine test_regression

  subroutine assert_close(actual, expected, tol, label)
    real(dp), intent(in) :: actual, expected, tol
    character(len=*), intent(in) :: label
    if (abs(actual-expected) > tol) then
      write(*,'(a,2es24.15)') trim(label)//' failed: ', actual, expected
      error stop 1
    end if
  end subroutine assert_close

  subroutine assert_close_vec(actual, expected, tol, label)
    real(dp), intent(in) :: actual(:), expected(:), tol
    character(len=*), intent(in) :: label
    if (size(actual) /= size(expected)) then
      write(*,'(a)') trim(label)//' size failed'
      error stop 1
    end if
    if (size(actual) > 0 .and. maxval(abs(actual-expected)) > tol) then
      write(*,'(a)') trim(label)//' failed'
      write(*,'(*(es18.9,1x))') actual
      write(*,'(*(es18.9,1x))') expected
      error stop 1
    end if
  end subroutine assert_close_vec

  subroutine assert_int_vec(actual, expected, label)
    integer, intent(in) :: actual(:), expected(:)
    character(len=*), intent(in) :: label
    if (size(actual) /= size(expected) .or. any(actual /= expected)) then
      write(*,'(a)') trim(label)//' failed'
      error stop 1
    end if
  end subroutine assert_int_vec

  subroutine assert_true(condition, label)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: label
    if (.not. condition) then
      write(*,'(a)') trim(label)//' failed'
      error stop 1
    end if
  end subroutine assert_true

end program test_pa
