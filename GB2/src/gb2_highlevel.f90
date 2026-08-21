! SPDX-License-Identifier: GPL-2.0-or-later
module gb2_highlevel
  use gb2_kinds, only : dp
  use gb2_optimizer, only : optimization_result
  use gb2_likelihood, only : fit_gb2_full, fit_gb2_profile
  use gb2_indicators, only : main_gb2
  use gb2_empirical, only : main_emp
  implicit none
  private
  public :: mlfit_gb2
contains
  subroutine mlfit_gb2(x,w,full_fit,profile_fit,empirical,full_indicators,profile_indicators)
    real(dp), intent(in) :: x(:),w(:)
    type(optimization_result), intent(out) :: full_fit,profile_fit
    real(dp), intent(out) :: empirical(6),full_indicators(6),profile_indicators(6)
    real(dp), allocatable :: xp(:),wp(:)
    integer :: n
    n=count(x>0.0_dp)
    if(n<=0) error stop 'mlfit_gb2: no positive observations'
    allocate(xp(n),wp(n))
    xp=pack(x,x>0.0_dp)
    wp=pack(w,x>0.0_dp)
    call fit_gb2_full(xp,full_fit,wp)
    call fit_gb2_profile(xp,profile_fit,wp)
    call main_emp(xp,wp,empirical)
    call main_gb2(0.6_dp,full_fit%par(1),full_fit%par(2),full_fit%par(3),full_fit%par(4),full_indicators)
    call main_gb2(0.6_dp,profile_fit%par(1),profile_fit%par(2),profile_fit%par(3),profile_fit%par(4),profile_indicators)
  end subroutine mlfit_gb2
end module gb2_highlevel
