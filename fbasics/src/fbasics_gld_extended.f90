! SPDX-License-Identifier: GPL-2.0-or-later
! Numerical translation derived from the GPL-2-or-later fBasics package.
module fbasics_gld_extended
  use fbasics_kinds, only: dp, clamp
  use fbasics_rng, only: runif_lcg
  use fbasics_stats, only: sample_mean, sample_sd, sample_quantile
  use fbasics_optimize, only: nelder_mead_bounded, numerical_hessian
  use fbasics_linalg, only: matrix_inverse
  use fbasics_distributions, only: distribution_fit
  implicit none
  private
  public :: qgld_fmkl, pgld_fmkl, dgld_fmkl, rgld_fmkl
  public :: qgld_fm5, pgld_fm5, dgld_fm5, rgld_fm5
  public :: fit_gld_extended
  real(dp), allocatable, save :: active_data(:)
  character(len=8), save :: active_family='fmkl'
  character(len=12), save :: active_method='rob'
contains
  pure real(dp) function power_term(u,a) result(v)
    real(dp), intent(in) :: u,a
    if (abs(a)<1.0e-10_dp) then
      v=log(u)
    else
      v=(u**a-1.0_dp)/a
    end if
  end function power_term

  pure real(dp) function qgld_fmkl(p,l1,l2,l3,l4) result(v)
    real(dp), intent(in) :: p,l1,l2,l3,l4
    if (p<0.0_dp .or. p>1.0_dp .or. l2<=0.0_dp) then
      v=huge(1.0_dp); return
    end if
    if (p==0.0_dp) then
      if (l3>0.0_dp) then; v=l1+(-1.0_dp/l3-power_term(1.0_dp,l4))/l2; else; v=-huge(1.0_dp); end if
      return
    else if (p==1.0_dp) then
      if (l4>0.0_dp) then; v=l1+(power_term(1.0_dp,l3)+1.0_dp/l4)/l2; else; v=huge(1.0_dp); end if
      return
    end if
    v=l1+(power_term(p,l3)-power_term(1.0_dp-p,l4))/l2
  end function qgld_fmkl

  pure real(dp) function qgld_fm5(p,l1,l2,l3,l4,l5) result(v)
    real(dp), intent(in) :: p,l1,l2,l3,l4,l5
    real(dp) :: u
    if (p<0.0_dp .or. p>1.0_dp .or. l2<=0.0_dp .or. abs(l5)>1.0_dp) then
      v=huge(1.0_dp); return
    end if
    u=clamp(p,1.0e-14_dp,1.0_dp-1.0e-14_dp)
    v=l1+((1.0_dp-l5)*power_term(u,l3)-(1.0_dp+l5)*power_term(1.0_dp-u,l4))/l2
  end function qgld_fm5

  real(dp) function pgld_fmkl(x,l1,l2,l3,l4) result(v)
    real(dp), intent(in) :: x,l1,l2,l3,l4
    real(dp) :: lo,hi,mid
    integer :: i
    if (l2<=0.0_dp) then; v=0.0_dp; return; end if
    lo=1.0e-12_dp; hi=1.0_dp-1.0e-12_dp
    if (x<=qgld_fmkl(lo,l1,l2,l3,l4)) then; v=0.0_dp; return; end if
    if (x>=qgld_fmkl(hi,l1,l2,l3,l4)) then; v=1.0_dp; return; end if
    do i=1,70
      mid=0.5_dp*(lo+hi)
      if (qgld_fmkl(mid,l1,l2,l3,l4)<x) then; lo=mid; else; hi=mid; end if
    end do
    v=0.5_dp*(lo+hi)
  end function pgld_fmkl

  real(dp) function pgld_fm5(x,l1,l2,l3,l4,l5) result(v)
    real(dp), intent(in) :: x,l1,l2,l3,l4,l5
    real(dp) :: lo,hi,mid
    integer :: i
    if (l2<=0.0_dp .or. abs(l5)>1.0_dp) then; v=0.0_dp; return; end if
    lo=1.0e-12_dp; hi=1.0_dp-1.0e-12_dp
    if (x<=qgld_fm5(lo,l1,l2,l3,l4,l5)) then; v=0.0_dp; return; end if
    if (x>=qgld_fm5(hi,l1,l2,l3,l4,l5)) then; v=1.0_dp; return; end if
    do i=1,70
      mid=0.5_dp*(lo+hi)
      if (qgld_fm5(mid,l1,l2,l3,l4,l5)<x) then; lo=mid; else; hi=mid; end if
    end do
    v=0.5_dp*(lo+hi)
  end function pgld_fm5

  real(dp) function dgld_fmkl(x,l1,l2,l3,l4) result(v)
    real(dp), intent(in) :: x,l1,l2,l3,l4
    real(dp) :: u,den
    if (l2<=0.0_dp) then; v=0.0_dp; return; end if
    u=pgld_fmkl(x,l1,l2,l3,l4)
    if (u<=0.0_dp .or. u>=1.0_dp) then; v=0.0_dp; return; end if
    den=u**(l3-1.0_dp)+(1.0_dp-u)**(l4-1.0_dp)
    if (den>0.0_dp) then; v=l2/den; else; v=0.0_dp; end if
  end function dgld_fmkl

  real(dp) function dgld_fm5(x,l1,l2,l3,l4,l5) result(v)
    real(dp), intent(in) :: x,l1,l2,l3,l4,l5
    real(dp) :: u,den
    if (l2<=0.0_dp .or. abs(l5)>1.0_dp) then; v=0.0_dp; return; end if
    u=pgld_fm5(x,l1,l2,l3,l4,l5)
    if (u<=0.0_dp .or. u>=1.0_dp) then; v=0.0_dp; return; end if
    den=(1.0_dp-l5)*u**(l3-1.0_dp)+(1.0_dp+l5)*(1.0_dp-u)**(l4-1.0_dp)
    if (den>0.0_dp) then; v=l2/den; else; v=0.0_dp; end if
  end function dgld_fm5

  real(dp) function rgld_fmkl(l1,l2,l3,l4) result(v)
    real(dp), intent(in) :: l1,l2,l3,l4
    v=qgld_fmkl(runif_lcg(),l1,l2,l3,l4)
  end function rgld_fmkl

  real(dp) function rgld_fm5(l1,l2,l3,l4,l5) result(v)
    real(dp), intent(in) :: l1,l2,l3,l4,l5
    v=qgld_fm5(runif_lcg(),l1,l2,l3,l4,l5)
  end function rgld_fm5

  subroutine fit_gld_extended(x,family,method,fit,max_iter)
    real(dp), intent(in) :: x(:)
    character(len=*), intent(in) :: family,method
    type(distribution_fit), intent(out) :: fit
    integer, intent(in), optional :: max_iter
    real(dp), allocatable :: start(:),lower(:),upper(:),best(:),h(:,:),cov(:,:)
    real(dp) :: s,m,fbest
    logical :: conv
    integer :: info
    integer :: npar,imax
    active_family=adjustl(family); active_method=adjustl(method)
    if (allocated(active_data)) deallocate(active_data)
    allocate(active_data(size(x))); active_data=x
    s=max(sample_sd(x),1.0e-5_dp); m=sample_quantile(x,0.5_dp)
    if (trim(active_family)=='fm5') then
      npar=5
      allocate(start(npar),lower(npar),upper(npar))
      start=[m,1.0_dp/s,0.1_dp,0.1_dp,0.0_dp]
      lower=[minval(x)-5*s,1.0e-5_dp,-5.0_dp,-5.0_dp,-0.99_dp]
      upper=[maxval(x)+5*s,100.0_dp/s,5.0_dp,5.0_dp,0.99_dp]
    else
      npar=4
      allocate(start(npar),lower(npar),upper(npar))
      start=[m,1.0_dp/s,0.1_dp,0.1_dp]
      lower=[minval(x)-5*s,1.0e-5_dp,-5.0_dp,-5.0_dp]
      upper=[maxval(x)+5*s,100.0_dp/s,5.0_dp,5.0_dp]
    end if
    imax=1800; if (present(max_iter)) imax=max_iter
    call nelder_mead_bounded(gld_extended_objective,start,lower,upper,best,fbest,conv,max_iter=imax,tol=1.0e-8_dp)
    fit%family='gld-'//trim(active_family)
    fit%parameters=best; fit%loglik=-gld_loglik(best); fit%converged=conv
    fit%aic=2.0_dp*npar-2.0_dp*fit%loglik
    fit%bic=log(real(size(x),dp))*npar-2.0_dp*fit%loglik
    call numerical_hessian(gld_extended_objective,best,h)
    call matrix_inverse(h,cov,info)
    fit%hessian=h
    if (info==0) fit%covariance=cov
  end subroutine fit_gld_extended

  real(dp) function gld_extended_objective(p) result(v)
    real(dp), intent(in) :: p(:)
    real(dp), allocatable :: sx(:),cdf(:),sp(:)
    real(dp) :: qobs(5),qmod(5),probs(5),hwidth,lo,hi,center,expected,dens
    integer :: i,j,nb,n
    if (p(2)<=0.0_dp) then
      v=huge(1.0_dp); return
    end if
    if (size(p)==5) then
      if (abs(p(5))>=1.0_dp) then
        v=huge(1.0_dp); return
      end if
    end if
    n=size(active_data)
    select case(trim(active_method))
    case('mle')
      v=gld_loglik(p)
    case('mps')
      allocate(sx(n),cdf(n),sp(n+1)); sx=active_data; call sort_values(sx)
      do i=1,n; cdf(i)=gld_cdf(sx(i),p); end do
      sp(1)=cdf(1); do i=2,n; sp(i)=cdf(i)-cdf(i-1); end do; sp(n+1)=1.0_dp-cdf(n)
      if (minval(sp)<=0.0_dp) then; v=huge(1.0_dp); else; v=-sum(log(sp)); end if
    case('gof')
      allocate(sx(n)); sx=active_data; call sort_values(sx); v=0.0_dp
      do i=1,n
        v=max(v,max(abs(gld_cdf(sx(i),p)-real(i,dp)/n),abs(gld_cdf(sx(i),p)-real(i-1,dp)/n)))
      end do
    case('hist')
      hwidth=2.0_dp*(sample_quantile(active_data,0.75_dp)-sample_quantile(active_data,0.25_dp))/real(n,dp)**(1.0_dp/3.0_dp)
      if (hwidth<=0.0_dp) hwidth=3.5_dp*sample_sd(active_data)/real(n,dp)**(1.0_dp/3.0_dp)
      nb=max(5,ceiling((maxval(active_data)-minval(active_data))/max(hwidth,1.0e-8_dp)))
      lo=minval(active_data); hi=maxval(active_data); hwidth=(hi-lo)/real(nb,dp); v=0.0_dp
      do j=1,nb
        center=lo+(real(j,dp)-0.5_dp)*hwidth; expected=0.0_dp
        do i=1,n
          if (active_data(i)>=lo+(j-1)*hwidth .and. (active_data(i)<lo+j*hwidth .or. j==nb)) expected=expected+1.0_dp
        end do
        expected=expected/(real(n,dp)*hwidth); dens=gld_density(center,p); v=v+(expected-dens)**2
      end do
    case default
      probs=[0.05_dp,0.25_dp,0.50_dp,0.75_dp,0.95_dp]
      do i=1,5; qobs(i)=sample_quantile(active_data,probs(i)); qmod(i)=gld_quantile(probs(i),p); end do
      v=sum((qobs-qmod)**2)/(sample_sd(active_data)**2+1.0e-12_dp)
    end select
  end function gld_extended_objective

  real(dp) function gld_loglik(p) result(v)
    real(dp), intent(in) :: p(:)
    integer :: i
    real(dp) :: d
    v=0.0_dp
    do i=1,size(active_data)
      d=gld_density(active_data(i),p)
      if (d<=tiny(1.0_dp) .or. d/=d) then; v=huge(1.0_dp); return; end if
      v=v-log(d)
    end do
  end function gld_loglik

  real(dp) function gld_quantile(prob,p) result(v)
    real(dp), intent(in) :: prob,p(:)
    if (size(p)==5) then; v=qgld_fm5(prob,p(1),p(2),p(3),p(4),p(5)); else; v=qgld_fmkl(prob,p(1),p(2),p(3),p(4)); end if
  end function gld_quantile
  real(dp) function gld_cdf(x,p) result(v)
    real(dp), intent(in) :: x,p(:)
    if (size(p)==5) then; v=pgld_fm5(x,p(1),p(2),p(3),p(4),p(5)); else; v=pgld_fmkl(x,p(1),p(2),p(3),p(4)); end if
  end function gld_cdf
  real(dp) function gld_density(x,p) result(v)
    real(dp), intent(in) :: x,p(:)
    if (size(p)==5) then; v=dgld_fm5(x,p(1),p(2),p(3),p(4),p(5)); else; v=dgld_fmkl(x,p(1),p(2),p(3),p(4)); end if
  end function gld_density
  subroutine sort_values(a)
    real(dp), intent(inout) :: a(:)
    integer :: i,j
    real(dp) :: key
    do i=2,size(a); key=a(i); j=i-1; do while(j>=1); if (a(j)<=key) exit; a(j+1)=a(j); j=j-1; end do; a(j+1)=key; end do
  end subroutine sort_values
end module fbasics_gld_extended
