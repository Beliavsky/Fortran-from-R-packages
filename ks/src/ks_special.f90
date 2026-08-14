! SPDX-License-Identifier: GPL-2.0-only
module ks_special
  use ks_kinds, only: dp
  implicit none
  private
  public :: regularized_gamma_p, chi_square_cdf
contains
  pure real(dp) function regularized_gamma_p(a,x) result(p)
    real(dp), intent(in) :: a,x
    real(dp) :: ap,del,sumg,b,c,d,h,an,eps,fpmin,gln
    integer :: n
    if(a<=0.0_dp .or. x<0.0_dp) then;p=0.0_dp;return;end if
    if(x<=0.0_dp) then;p=0.0_dp;return;end if
    eps=5.0e-15_dp;fpmin=tiny(1.0_dp)/eps;gln=log_gamma(a)
    if(x<a+1.0_dp) then
      ap=a;sumg=1.0_dp/a;del=sumg
      do n=1,10000
        ap=ap+1.0_dp;del=del*x/ap;sumg=sumg+del
        if(abs(del)<=abs(sumg)*eps)exit
      end do
      p=sumg*exp(-x+a*log(x)-gln)
    else
      b=x+1.0_dp-a;c=1.0_dp/fpmin;d=1.0_dp/b;h=d
      do n=1,10000
        an=-real(n,dp)*(real(n,dp)-a);b=b+2.0_dp
        d=an*d+b;if(abs(d)<fpmin)d=fpmin
        c=b+an/c;if(abs(c)<fpmin)c=fpmin
        d=1.0_dp/d;del=d*c;h=h*del
        if(abs(del-1.0_dp)<=eps)exit
      end do
      p=1.0_dp-exp(-x+a*log(x)-gln)*h
    end if
    p=min(1.0_dp,max(0.0_dp,p))
  end function

  pure real(dp) function chi_square_cdf(x,df) result(p)
    real(dp),intent(in)::x
    integer,intent(in)::df
    if(df<=0)then;p=0.0_dp;else;p=regularized_gamma_p(0.5_dp*real(df,dp),0.5_dp*x);end if
  end function
end module ks_special
