module trawl_intersection
  use trawl_kinds, only : dp
  use trawl_types, only : trawl_spec
  use trawl_functions, only : eval_trawl,trawl_tail_area
  implicit none
  private
  public :: fit_trawl_intersection,fit_trawl_intersection_exp,fit_trawl_intersection_lm
contains
  real(dp) function fit_trawl_intersection(spec1,spec2,lm1,lm2) result(r12)
    type(trawl_spec),intent(in)::spec1,spec2
    real(dp),intent(in)::lm1,lm2
    real(dp)::roots(8),lb,up,mid
    integer::nr
    lb=-1000.0_dp; up=-1.0e-6_dp
    call all_roots(spec1,spec2,lb,up,roots,nr)
    if(nr==0 .or. nr>=3) then
      r12=min(lm1,lm2); return
    end if
    mid=0.5_dp*(lb+roots(1))
    if(nr==1) then
      if(diff(mid,spec1,spec2)<=0.0_dp) then
        r12=trawl_tail_area(roots(1),spec1) + &
          trawl_tail_area(0.0_dp,spec2)-trawl_tail_area(roots(1),spec2)
      else
        r12=trawl_tail_area(roots(1),spec2) + &
          trawl_tail_area(0.0_dp,spec1)-trawl_tail_area(roots(1),spec1)
      end if
    else
      if(diff(mid,spec1,spec2)<=0.0_dp) then
        r12=trawl_tail_area(roots(1),spec1) + &
          trawl_tail_area(roots(2),spec2)-trawl_tail_area(roots(1),spec2) + &
          trawl_tail_area(0.0_dp,spec1)-trawl_tail_area(roots(2),spec1)
      else
        r12=trawl_tail_area(roots(1),spec2) + &
          trawl_tail_area(roots(2),spec1)-trawl_tail_area(roots(1),spec1) + &
          trawl_tail_area(0.0_dp,spec2)-trawl_tail_area(roots(2),spec2)
      end if
    end if
  end function

  real(dp) function fit_trawl_intersection_exp(lambda1,lambda2,lm1,lm2) result(r12)
    real(dp),intent(in)::lambda1,lambda2,lm1,lm2
    type(trawl_spec)::s1,s2
    s1%kind='Exp';s1%lambda1=lambda1
    s2%kind='Exp';s2%lambda1=lambda2
    r12=fit_trawl_intersection(s1,s2,lm1,lm2)
  end function

  real(dp) function fit_trawl_intersection_lm(alpha1,h1,alpha2,h2,lm1,lm2) result(r12)
    real(dp),intent(in)::alpha1,h1,alpha2,h2,lm1,lm2
    type(trawl_spec)::s1,s2
    s1%kind='LM';s1%alpha=alpha1;s1%h=h1
    s2%kind='LM';s2%alpha=alpha2;s2%h=h2
    r12=fit_trawl_intersection(s1,s2,lm1,lm2)
  end function

  real(dp) function diff(x,s1,s2) result(v)
    real(dp),intent(in)::x
    type(trawl_spec),intent(in)::s1,s2
    v=eval_trawl(x,s1)-eval_trawl(x,s2)
  end function

  subroutine all_roots(s1,s2,lo,hi,roots,nr)
    type(trawl_spec),intent(in)::s1,s2
    real(dp),intent(in)::lo,hi
    real(dp),intent(out)::roots(:)
    integer,intent(out)::nr
    integer,parameter::ngrid=5000
    integer::i
    real(dp)::x0,x1,f0,f1,dx,r,tol
    nr=0; roots=0.0_dp; dx=(hi-lo)/real(ngrid,dp); tol=1.0e-10_dp
    x0=lo;f0=diff(x0,s1,s2)
    do i=1,ngrid
      x1=lo+real(i,dp)*dx;f1=diff(x1,s1,s2)
      if(abs(f0)<=tol .and. abs(f1)<=tol) then
        x0=x1;f0=f1;cycle
      end if
      if(f0*f1<0.0_dp) then
        r=bisect_root(s1,s2,x0,x1)
        if(nr==0) then
          nr=1; roots(nr)=r
        else if(abs(r-roots(nr))>1.0e-6_dp) then
          nr=nr+1; if(nr<=size(roots)) roots(nr)=r
        end if
        if(nr>=size(roots)) return
      else if(abs(f1)<=tol .and. abs(f0)>tol) then
        nr=nr+1; if(nr<=size(roots)) roots(nr)=x1
        if(nr>=size(roots)) return
      end if
      x0=x1;f0=f1
    end do
  end subroutine

  real(dp) function bisect_root(s1,s2,a0,b0) result(r)
    type(trawl_spec),intent(in)::s1,s2
    real(dp),intent(in)::a0,b0
    real(dp)::a,b,c,fa,fc
    integer::i
    a=a0;b=b0;fa=diff(a,s1,s2)
    do i=1,100
      c=0.5_dp*(a+b);fc=diff(c,s1,s2)
      if(abs(fc)<1.0e-13_dp .or. abs(b-a)<1.0e-12_dp) exit
      if(fa*fc<=0.0_dp) then
        b=c
      else
        a=c;fa=fc
      end if
    end do
    r=0.5_dp*(a+b)
  end function
end module trawl_intersection
