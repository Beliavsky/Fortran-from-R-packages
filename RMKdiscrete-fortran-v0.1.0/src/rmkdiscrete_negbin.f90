! SPDX-License-Identifier: GPL-2.0-or-later
! Modern Fortran translation of RMKdiscrete 0.1 by Robert M. Kirkpatrick.
! See LICENSE, COPYING, and TRANSLATION_NOTES.md for provenance.
module rmkdiscrete_negbin
  use rmkdiscrete_kinds, only : dp
  use rmkdiscrete_math, only : qnan, real_equal, dnbinom_prob, rnbinom_prob
  implicit none
  private
  public :: dnegbin, negbin_from_nu_p, negbin_from_mu_sigma2, negbin_from_mu_nu, negbin_from_mu_p
  public :: negbin_from_sigma2_p
  public :: rnegbin, rnegbin_sample
contains

  real(dp) function dnegbin(x, nu, p, mu, give_log) result(v)
    integer, intent(in) :: x
    real(dp), intent(in), optional :: nu, p, mu
    logical, intent(in), optional :: give_log
    real(dp) :: pp, nnu
    integer :: nprovided
    nprovided=0
    if(present(nu)) nprovided=nprovided+1
    if(present(p)) nprovided=nprovided+1
    if(present(mu)) nprovided=nprovided+1
    if(nprovided/=2) then
      v=qnan()
      return
    end if
    if(present(nu)) nnu=nu
    if(present(p)) pp=p
    if(.not.present(nu)) then
      if(p<=0.0_dp .or. p>=1.0_dp .or. mu<0.0_dp) then
        v=qnan()
        return
      end if
      nnu=mu*p/(1.0_dp-p)
    else if(.not.present(p)) then
      if(nnu<0.0_dp .or. mu<0.0_dp) then
        v=qnan()
        return
      end if
      if(real_equal(nnu,0.0_dp) .and. real_equal(mu,0.0_dp)) then
        pp=1.0_dp
      else
        pp=nnu/(nnu+mu)
      end if
    end if
    v=dnbinom_prob(x,nnu,pp,give_log)
  end function dnegbin

  pure subroutine negbin_from_nu_p(nu,p,mu,sigma2,status)
    real(dp), intent(in) :: nu,p
    real(dp), intent(out) :: mu,sigma2
    integer, intent(out), optional :: status
    if (nu < 0.0_dp .or. p <= 0.0_dp .or. p > 1.0_dp) then
      mu=qnan()
      sigma2=qnan()
      if(present(status)) status=1
      return
    end if
    mu=nu*(1.0_dp-p)/p
    sigma2=nu*(1.0_dp-p)/p**2
    if(present(status)) status=0
  end subroutine negbin_from_nu_p

  pure subroutine negbin_from_mu_sigma2(mu,sigma2,nu,p,status)
    real(dp), intent(in) :: mu,sigma2
    real(dp), intent(out) :: nu,p
    integer, intent(out), optional :: status
    if (mu < 0.0_dp .or. sigma2 < mu .or. sigma2 <= 0.0_dp) then
      nu=qnan()
      p=qnan()
      if(present(status)) status=1
      return
    end if
    p=mu/sigma2
    if (real_equal(p,1.0_dp)) then
      nu=huge(1.0_dp)
    else
      nu=mu*p/(1.0_dp-p)
    end if
    if(present(status)) status=0
  end subroutine negbin_from_mu_sigma2

  pure subroutine negbin_from_mu_nu(mu,nu,sigma2,p,status)
    real(dp), intent(in) :: mu,nu
    real(dp), intent(out) :: sigma2,p
    integer, intent(out), optional :: status
    if (mu < 0.0_dp .or. nu <= 0.0_dp) then
      sigma2=qnan()
      p=qnan()
      if(present(status)) status=1
      return
    end if
    sigma2=mu+mu**2/nu
    p=nu/(nu+mu)
    if(present(status)) status=0
  end subroutine negbin_from_mu_nu

  pure subroutine negbin_from_sigma2_p(sigma2,p,mu,nu,status)
    real(dp), intent(in) :: sigma2,p
    real(dp), intent(out) :: mu,nu
    integer, intent(out), optional :: status
    if(sigma2<0.0_dp .or. p<=0.0_dp .or. p>=1.0_dp) then
      mu=qnan()
      nu=qnan()
      if(present(status)) status=1
      return
    end if
    mu=p*sigma2
    nu=mu*p/(1.0_dp-p)
    if(present(status)) status=0
  end subroutine negbin_from_sigma2_p

  pure subroutine negbin_from_mu_p(mu,p,sigma2,nu,status)
    real(dp), intent(in) :: mu,p
    real(dp), intent(out) :: sigma2,nu
    integer, intent(out), optional :: status
    if (mu < 0.0_dp .or. p <= 0.0_dp .or. p >= 1.0_dp) then
      sigma2=qnan()
      nu=qnan()
      if(present(status)) status=1
      return
    end if
    sigma2=mu/p
    nu=mu*p/(1.0_dp-p)
    if(present(status)) status=0
  end subroutine negbin_from_mu_p

  integer function rnegbin(nu,p) result(x)
    real(dp), intent(in) :: nu,p
    x=rnbinom_prob(nu,p)
  end function rnegbin

  function rnegbin_sample(n,nu,p) result(out)
    integer, intent(in) :: n
    real(dp), intent(in) :: nu,p
    integer, allocatable :: out(:)
    integer :: i
    allocate(out(max(0,n)))
    do i=1,n
      out(i)=rnegbin(nu,p)
    end do
  end function rnegbin_sample
end module rmkdiscrete_negbin
