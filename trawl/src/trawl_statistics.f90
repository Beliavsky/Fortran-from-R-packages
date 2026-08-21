module trawl_statistics
  use trawl_kinds, only : dp
  implicit none
  private
  public :: sample_mean,sample_variance,empirical_acf
contains
  pure real(dp) function sample_mean(x) result(m)
    real(dp),intent(in)::x(:)
    if(size(x)==0) then; m=0.0_dp; else; m=sum(x)/real(size(x),dp); end if
  end function
  pure real(dp) function sample_variance(x) result(v)
    real(dp),intent(in)::x(:)
    real(dp)::m
    if(size(x)<2) then; v=0.0_dp; return; end if
    m=sample_mean(x); v=sum((x-m)**2)/real(size(x)-1,dp)
  end function
  subroutine empirical_acf(x,lag_max,acf)
    real(dp),intent(in)::x(:)
    integer,intent(in)::lag_max
    real(dp),allocatable,intent(out)::acf(:)
    real(dp)::m,den
    integer::k,n,lm
    n=size(x); lm=max(0,min(lag_max,n-1)); allocate(acf(lm))
    if(lm==0) return
    m=sample_mean(x); den=sum((x-m)**2)
    if(den<=tiny(1.0_dp)) then; acf=0.0_dp; return; end if
    do k=1,lm
      acf(k)=sum((x(1:n-k)-m)*(x(1+k:n)-m))/den
    end do
  end subroutine
end module trawl_statistics
