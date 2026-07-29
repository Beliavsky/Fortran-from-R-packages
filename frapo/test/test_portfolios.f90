! SPDX-License-Identifier: GPL-3.0-or-later
program test_portfolios
  use frapo
  implicit none

  real(dp) :: returns(20, 3), covariance(3, 3), sigma_w(3), total_risk
  real(dp), parameter :: expected_gmv(3) = [ &
    0.14621071595666_dp, 0.62650118429361_dp, 0.22728809974973_dp]
  real(dp), parameter :: expected_md(3) = [ &
    0.24110697_dp, 0.45974267_dp, 0.29915036_dp]
  type(portfolio_result) :: gmv, md, mtd, erc
  integer :: i

  do i = 1, 20
    returns(i, 1) = 0.01_dp * sin(real(i, dp))
    returns(i, 2) = 0.005_dp * cos(real(i, dp))
    returns(i, 3) = 0.008_dp * sin(0.7_dp * real(i, dp))
  end do

  gmv = pgmv(returns, percentage=.false.)
  call assert_equal_int(gmv%status, frapo_ok, 'PGMV status')
  call assert_vector(gmv%weights, expected_gmv, 5.0e-9_dp, 'PGMV reference')

  md = pmd(returns, percentage=.false.)
  call assert_equal_int(md%status, frapo_ok, 'PMD status')
  call assert_vector(md%weights, expected_md, 2.0e-7_dp, 'PMD reference')

  mtd = pmtd(returns, method=tdc_empirical, k=4, percentage=.false.)
  call assert_equal_int(mtd%status, frapo_ok, 'PMTD status')
  call assert_close(sum(mtd%weights), 1.0_dp, 1.0e-10_dp, 'PMTD budget')
  call assert_true(minval(mtd%weights) >= -1.0e-10_dp, 'PMTD nonnegative')

  covariance = reshape([0.04_dp, 0.006_dp, 0.004_dp, &
                        0.006_dp, 0.09_dp, 0.01_dp, &
                        0.004_dp, 0.01_dp, 0.16_dp], [3, 3])
  erc = perc(covariance, percentage=.false.)
  call assert_equal_int(erc%status, frapo_ok, 'PERC status')
  call assert_close(sum(erc%weights), 1.0_dp, 1.0e-12_dp, 'PERC budget')
  sigma_w = matmul(covariance, erc%weights)
  total_risk = sqrt(dot_product(erc%weights, sigma_w))
  sigma_w = erc%weights * sigma_w / total_risk
  call assert_close(maxval(sigma_w) - minval(sigma_w), 0.0_dp, 2.0e-10_dp, 'equal risk contributions')

  print '(a)', 'test_portfolios: PASS'

contains
  subroutine assert_close(actual, expected, tolerance, label)
    real(dp), intent(in) :: actual, expected, tolerance
    character(len=*), intent(in) :: label
    if (abs(actual - expected) > tolerance) then
      write(*, '(a,2es24.15)') trim(label)//': ', actual, expected
      error stop 1
    end if
  end subroutine assert_close

  subroutine assert_vector(actual, expected, tolerance, label)
    real(dp), intent(in) :: actual(:), expected(:), tolerance
    character(len=*), intent(in) :: label
    if (maxval(abs(actual - expected)) > tolerance) then
      write(*, '(a)') trim(label)//': values differ'
      error stop 1
    end if
  end subroutine assert_vector

  subroutine assert_true(condition, label)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: label
    if (.not. condition) then
      write(*, '(a)') trim(label)//': failed'
      error stop 1
    end if
  end subroutine assert_true

  subroutine assert_equal_int(actual, expected, label)
    integer, intent(in) :: actual, expected
    character(len=*), intent(in) :: label
    if (actual /= expected) then
      write(*, '(a,2i8)') trim(label)//': ', actual, expected
      error stop 1
    end if
  end subroutine assert_equal_int
end program test_portfolios
