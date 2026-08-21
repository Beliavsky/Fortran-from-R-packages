module relsurv_stats
  use relsurv_kinds, only : dp
  implicit none
  private
  public :: chi_square_sf, normal_cdf
contains
  pure function normal_cdf(x) result(p)
    real(dp),intent(in)::x
    real(dp)::p
    p=0.5_dp*erfc(-x/sqrt(2.0_dp))
  end function normal_cdf

  function chi_square_sf(x,df) result(p)
    real(dp),intent(in)::x
    integer,intent(in)::df
    real(dp)::p
    if(x<=0.0_dp)then; p=1.0_dp; return; end if
    p=gammq(0.5_dp*real(df,dp),0.5_dp*x)
  end function chi_square_sf

  function gammq(a,x) result(q)
    real(dp),intent(in)::a,x
    real(dp)::q,gln,gamser,gammcf
    if(x<0.0_dp .or. a<=0.0_dp)then; q=0.0_dp; return; end if
    gln=log_gamma(a)
    if(x<a+1.0_dp)then
      call gser(gamser,a,x,gln); q=1.0_dp-gamser
    else
      call gcf(gammcf,a,x,gln); q=gammcf
    end if
    q=max(0.0_dp,min(1.0_dp,q))
  end function gammq

  subroutine gser(gamser,a,x,gln)
    real(dp),intent(out)::gamser
    real(dp),intent(in)::a,x,gln
    integer,parameter::itmax=1000
    real(dp),parameter::eps=3.0e-14_dp
    real(dp)::ap,del,summ
    integer::n
    if(x<=0.0_dp)then; gamser=0.0_dp; return; end if
    ap=a; summ=1.0_dp/a; del=summ
    do n=1,itmax
      ap=ap+1.0_dp; del=del*x/ap; summ=summ+del
      if(abs(del)<abs(summ)*eps)exit
    end do
    gamser=summ*exp(-x+a*log(x)-gln)
  end subroutine gser

  subroutine gcf(gammcf,a,x,gln)
    real(dp),intent(out)::gammcf
    real(dp),intent(in)::a,x,gln
    integer,parameter::itmax=1000
    real(dp),parameter::eps=3.0e-14_dp,fpmin=1.0e-300_dp
    real(dp)::an,b,c,d,del,h
    integer::i
    b=x+1.0_dp-a; c=1.0_dp/fpmin; d=1.0_dp/max(abs(b),fpmin); if(b<0.0_dp)d=-d
    h=d
    do i=1,itmax
      an=-real(i,dp)*(real(i,dp)-a); b=b+2.0_dp
      d=an*d+b; if(abs(d)<fpmin)d=fpmin
      c=b+an/c; if(abs(c)<fpmin)c=fpmin
      d=1.0_dp/d; del=d*c; h=h*del
      if(abs(del-1.0_dp)<eps)exit
    end do
    gammcf=exp(-x+a*log(x)-gln)*h
  end subroutine gcf
end module relsurv_stats
