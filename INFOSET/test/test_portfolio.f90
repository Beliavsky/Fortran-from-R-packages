! SPDX-License-Identifier: GPL-2.0-or-later
program test_portfolio
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use infoset, only : dp, ptf_construction, summary_ptf
  use infoset, only : portfolio_result, portfolio_summary, infoset_success
  implicit none
  integer, parameter :: n = 180, p = 4, nw = 6
  real(dp) :: prices(n,p), left_risk(p,nw), r
  type(portfolio_result) :: fit
  type(portfolio_summary) :: statistics
  integer :: i, j, t

  prices(1,:) = [100.0_dp, 75.0_dp, 130.0_dp, 90.0_dp]
  do i = 2, n
    do j = 1, p
      r = 0.0002_dp*real(j,dp) + 0.004_dp*sin(0.07_dp*real(i*j,dp)) &
        + 0.001_dp*cos(0.19_dp*real(i+j,dp))
      if (mod(i+3*j,53) == 0) r = r - 0.025_dp*real(j,dp)/real(p,dp)
      prices(i,j) = prices(i-1,j)*exp(r)
    end do
  end do
  do t = 1, nw
    do j = 1, p
      left_risk(j,t) = 0.01_dp*real(j,dp) + 0.001_dp*real(t,dp)
    end do
  end do

  call ptf_construction(prices, 80, 20, 'M', fit)
  call validate(fit, 'M')
  call summary_ptf(fit%oos_returns, statistics)
  call check(statistics%status == infoset_success, 'summary status')
  call check(statistics%count == size(fit%oos_returns), 'summary count')
  call check(statistics%minimum <= statistics%median .and. &
    statistics%median <= statistics%maximum, 'summary ordering')

  call ptf_construction(prices, 80, 20, 'C_M', fit, left_risk)
  call validate(fit, 'C_M')
  call ptf_construction(prices, 80, 20, 'EDC', fit)
  call validate(fit, 'EDC')
  call ptf_construction(prices, 80, 20, 'C_EDC', fit, left_risk)
  call validate(fit, 'C_EDC')
  print '(a)', 'test_portfolio: PASS'
contains
  subroutine validate(result, label)
    type(portfolio_result), intent(in) :: result
    character(len=*), intent(in) :: label
    integer :: k
    call check(result%status == infoset_success, trim(label)//' status')
    call check(all(shape(result%weights) == [p,nw-1]), trim(label)//' weights shape')
    call check(all(shape(result%oos_returns) == [21,nw-1]), trim(label)//' oos shape')
    call check(all(ieee_is_finite(result%oos_returns)), trim(label)//' finite')
    do k = 1, size(result%weights,2)
      call check(abs(sum(result%weights(:,k))-1.0_dp) < 1.0e-7_dp, &
        trim(label)//' full investment')
      call check(minval(result%weights(:,k)) >= -1.0e-8_dp, &
        trim(label)//' lower bound')
      call check(maxval(result%weights(:,k)) <= 1.0_dp+1.0e-8_dp, &
        trim(label)//' upper bound')
    end do
  end subroutine validate

  subroutine check(condition, message)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message
    if (.not. condition) then
      write(*,'(a)') 'FAIL: '//message
      error stop 1
    end if
  end subroutine check
end program test_portfolio
