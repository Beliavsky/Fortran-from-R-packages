! SPDX-License-Identifier: GPL-2.0-or-later
program test_extended_api
  use runuran
  implicit none
  type(continuous_distribution)::d
  type(discrete_distribution)::q
  integer::fails
  fails=0

  d=ud_continuous_cdf(ecdf,0.0_dp,10.0_dp,name='cdf-only exponential')
  call chk(d%cdf(log(2.0_dp)),0.5_dp,2e-12_dp,'cdf-only CDF')
  call chk(d%pdf(1.0_dp),exp(-1.0_dp),2e-7_dp,'cdf-only numerical PDF')

  d=ud_continuous_logpdf(elogpdf,0.0_dp,20.0_dp,name='logpdf-only exponential')
  call chk(d%pdf(1.0_dp),exp(-1.0_dp)/(1.0_dp-exp(-20.0_dp)),2e-9_dp,'logpdf-only PDF')
  call chk(d%dlogpdf(1.0_dp),-1.0_dp,2e-6_dp,'dlogpdf')

  q=ud_probability_vector([1.0_dp,2.0_dp,1.0_dp],lb=-1)
  call chk(q%pmf(0),0.5_dp,2e-15_dp,'probability vector PMF')
  if(q%quantile(0.7_dp)/=0)then
    print *,'probability vector quantile failed'
    fails=fails+1
  end if

  if(fails==0)then
    print *,'test_extended_api: PASS'
  else
    error stop 1
  end if
contains
  real(dp) function ecdf(x,p) result(y)
    real(dp),intent(in)::x,p(:)
    if(size(p)<0) stop
    if(x<=0.0_dp)then
      y=0.0_dp
    else
      y=1.0_dp-exp(-x)
    end if
  end function ecdf
  real(dp) function elogpdf(x,p) result(y)
    real(dp),intent(in)::x,p(:)
    if(size(p)<0) stop
    if(x<0.0_dp)then
      y=-huge(1.0_dp)
    else
      y=-x
    end if
  end function elogpdf
  subroutine chk(x,e,t,n)
    real(dp),intent(in)::x,e,t
    character(len=*),intent(in)::n
    if(abs(x-e)>t*max(1.0_dp,abs(e)))then
      print *,trim(n),' got=',x,' expected=',e
      fails=fails+1
    end if
  end subroutine chk
end program test_extended_api
