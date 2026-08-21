! SPDX-License-Identifier: GPL-2.0-or-later
module gb2_indicators
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use gb2_kinds, only : dp
  use gb2_special, only : log_beta, beta_fn, hypergeo3f2_1, quiet_nan
  use gb2_distribution, only : pgb2, qgb2
  use gb2_moments, only : moment_gb2, incomplete_moment_gb2
  implicit none
  private
  public :: arpt_gb2, arpr_gb2, rmpg_gb2, qsr_gb2, gb2_gini, gini_gb2
  public :: gini_b2, gini_dag, gini_sm, main_gb2, main2_gb2, thomae_gb2
contains
  real(dp) function arpt_gb2(prop,shape1,scale,shape2,shape3) result(v)
    real(dp), intent(in) :: prop,shape1,scale,shape2,shape3
    v=prop*qgb2(0.5_dp,shape1,scale,shape2,shape3)
  end function arpt_gb2

  real(dp) function arpr_gb2(prop,shape1,shape2,shape3) result(v)
    real(dp), intent(in) :: prop,shape1,shape2,shape3
    v=pgb2(arpt_gb2(prop,shape1,1.0_dp,shape2,shape3),shape1,1.0_dp,shape2,shape3)
  end function arpr_gb2

  real(dp) function rmpg_gb2(arpr,shape1,shape2,shape3) result(v)
    real(dp), intent(in) :: arpr,shape1,shape2,shape3
    real(dp) :: qlo,qhi
    qlo=qgb2(arpr/2.0_dp,shape1,1.0_dp,shape2,shape3)
    qhi=qgb2(arpr,shape1,1.0_dp,shape2,shape3)
    v=1.0_dp-qlo/qhi
  end function rmpg_gb2

  real(dp) function qsr_gb2(shape1,shape2,shape3) result(v)
    real(dp), intent(in) :: shape1,shape2,shape3
    real(dp) :: q20,q80,lo,hi
    q20=qgb2(0.2_dp,shape1,1.0_dp,shape2,shape3)
    q80=qgb2(0.8_dp,shape1,1.0_dp,shape2,shape3)
    lo=incomplete_moment_gb2(q20,1.0_dp,shape1,1.0_dp,shape2,shape3)
    hi=incomplete_moment_gb2(q80,1.0_dp,shape1,1.0_dp,shape2,shape3)
    v=(1.0_dp-hi)/lo
  end function qsr_gb2

  real(dp) function thomae_term(u,l,lb,tol,maxiter,ok) result(g1)
    real(dp), intent(in) :: u(3),l(2),lb,tol
    integer, intent(in) :: maxiter
    logical, intent(out) :: ok
    real(dp) :: g(5),ex(4),excess,exopt,lopt(2),uopt(3),f32,logg
    integer, parameter :: i1(4)=[1,1,2,3], i2(4)=[2,4,3,4]
    integer :: k,best,j,n
    excess=l(1)+l(2)-sum(u)
    g(1)=l(1)+l(2)-u(2)-u(3)
    g(2)=l(1)+l(2)-u(1)-u(3)
    g(3)=l(1)+l(2)-u(1)-u(2)
    g(4)=l(1)
    g(5)=l(2)
    if(excess<=0.0_dp) then
    g1=quiet_nan()
    ok=.false.
    return
    end if
    do k=1,4
    ex(k)=0.5_dp*sum(g)-g(i1(k))-g(i2(k))
    end do
    best=maxloc(ex,dim=1)
    exopt=ex(best)
    lopt=[g(i1(best)),g(i2(best))]
    n=0
    do j=1,5
      if(j/=i1(best) .and. j/=i2(best)) then
        n=n+1
        uopt(n)=g(j)-exopt
      end if
    end do
    f32=hypergeo3f2_1(uopt,lopt,tol=tol,maxiter=maxiter,converged=ok)
    if(.not.ok .or. f32<=0.0_dp) then
    g1=quiet_nan()
    return
    end if
    logg=sum(log_gamma(l)-log_gamma(lopt))+log_gamma(excess)-log_gamma(exopt)+log(f32)+lb
    if(logg>log(huge(1.0_dp))) then
    g1=huge(1.0_dp)
    else
    g1=exp(logg)
    end if
  end function thomae_term

  real(dp) function thomae_gb2(u,l,lb,tol,maxiter,converged) result(value)
    real(dp), intent(in) :: u(3),l(2),lb
    real(dp), intent(in), optional :: tol
    integer, intent(in), optional :: maxiter
    logical, intent(out), optional :: converged
    real(dp) :: t
    integer :: mi
    logical :: ok
    t=1.0e-8_dp
    if(present(tol)) t=tol
    mi=10000
    if(present(maxiter)) mi=maxiter
    value=thomae_term(u,l,lb,t,mi,ok)
    if(present(converged)) converged=ok
  end function thomae_gb2

  real(dp) function gb2_gini(shape1,shape2,shape3,tol,maxiter) result(gini)
    real(dp), intent(in) :: shape1,shape2,shape3
    real(dp), intent(in), optional :: tol
    integer, intent(in), optional :: maxiter
    real(dp) :: t,excess,u(3),l1(2),l2(2),lb,g1,g2
    integer :: mi
    logical :: ok1,ok2
    t=1.0e-8_dp
    if(present(tol)) t=tol
    mi=10000
    if(present(maxiter)) mi=maxiter
    if(shape1<=0.0_dp .or. shape2<=0.0_dp .or. shape3<=0.0_dp) then
    gini=quiet_nan()
    return
    end if
    excess=shape3-1.0_dp/shape1
    if(excess<=0.0_dp) then
    gini=quiet_nan()
    return
    end if
    if(excess<1.0e-10_dp) then
    gini=1.0_dp
    return
    end if
    u=[1.0_dp,shape2+shape3,2.0_dp*shape2+1.0_dp/shape1]
    l1=[shape2+1.0_dp,2.0_dp*(shape2+shape3)]
    l2=[shape2+1.0_dp+1.0_dp/shape1,2.0_dp*(shape2+shape3)]
    lb=log_beta(2.0_dp*shape3-1.0_dp/shape1,2.0_dp*shape2+1.0_dp/shape1)- &
      log_beta(shape2,shape3)-log_beta(shape2+1.0_dp/shape1,shape3-1.0_dp/shape1)
    g1=thomae_term(u,l1,lb,t,mi,ok1)
    g2=thomae_term(u,l2,lb,t,mi,ok2)
    if(.not.ok1 .or. .not.ok2 .or. .not.ieee_is_finite(g1)) then
    gini=1.0_dp
    return
    end if
    gini=g1/shape2-g2/(shape2+1.0_dp/shape1)
    if(gini>1.0_dp+1.0e-10_dp) gini=1.0_dp
    gini=max(0.0_dp,min(1.0_dp,gini))
  end function gb2_gini

  real(dp) function gini_gb2(shape1,shape2,shape3) result(v)
    real(dp), intent(in) :: shape1,shape2,shape3
    v=gb2_gini(shape1,shape2,shape3)
  end function gini_gb2

  pure real(dp) function gini_b2(shape2,shape3) result(v)
    real(dp), intent(in) :: shape2,shape3
    v=beta_fn(2.0_dp*shape2,2.0_dp*shape3-1.0_dp)/(beta_fn(shape2,shape3)**2)*(2.0_dp/shape2)
  end function gini_b2

  pure real(dp) function gini_dag(shape1,shape2) result(v)
    real(dp), intent(in) :: shape1,shape2
    v=exp(log_gamma(shape2)+log_gamma(2.0_dp*shape2+1.0_dp/shape1)-log_gamma(2.0_dp*shape2)-log_gamma(shape2+1.0_dp/shape1))-1.0_dp
  end function gini_dag

  pure real(dp) function gini_sm(shape1,shape3) result(v)
    real(dp), intent(in) :: shape1,shape3
    v=1.0_dp-exp(log_gamma(shape3)+log_gamma(2.0_dp*shape3-1.0_dp/shape1)-log_gamma(2.0_dp*shape3)-log_gamma(shape3-1.0_dp/shape1))
  end function gini_sm

  subroutine main_gb2(prop,shape1,scale,shape2,shape3,values)
    real(dp), intent(in) :: prop,shape1,scale,shape2,shape3
    real(dp), intent(out) :: values(6)
    real(dp) :: ar
    ar=arpr_gb2(prop,shape1,shape2,shape3)
    values(1)=qgb2(0.5_dp,shape1,scale,shape2,shape3)
    values(2)=moment_gb2(1.0_dp,shape1,scale,shape2,shape3)
    values(3)=100.0_dp*ar
    values(4)=100.0_dp*rmpg_gb2(ar,shape1,shape2,shape3)
    values(5)=qsr_gb2(shape1,shape2,shape3)
    values(6)=gini_gb2(shape1,shape2,shape3)
  end subroutine main_gb2

  subroutine main2_gb2(prop,shape1,scale,shape12,shape13,values)
    real(dp), intent(in) :: prop,shape1,scale,shape12,shape13
    real(dp), intent(out) :: values(6)
    call main_gb2(prop,shape1,scale,shape12/shape1,shape13/shape1,values)
  end subroutine main2_gb2
end module gb2_indicators
