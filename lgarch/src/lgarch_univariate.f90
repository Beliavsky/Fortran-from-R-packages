! SPDX-License-Identifier: GPL-2.0-only
! Copyright (C) 2014-2025 Genaro Sucarrat (original R/C package)
! Copyright (C) 2026 Modern Fortran translation contributors
! This file is distributed under the GNU General Public License version 2 only.
module lgarch_univariate
  use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan, ieee_is_finite
  use lgarch_kinds, only : dp, pi
  use lgarch_rng, only : random_normal_vector
  use lgarch_utils, only : mean_value, sample_variance, safe_log_mean_exp
  use lgarch_linalg, only : spectral_radius, inverse_matrix
  use lgarch_optimizer, only : minimize_nelder_mead, numerical_hessian
  implicit none
  private
  integer, parameter, public :: LGARCH_LS=1, LGARCH_ML=2, LGARCH_CEX2=3

  type, public :: lgarch_fit_result
    integer :: arch_order=0, garch_order=0, method=LGARCH_LS, iterations=0
    logical :: mean_correction=.false., converged=.false.
    real(dp) :: objective_arma=0.0_dp, loglik_model=0.0_dp, loglik_arma=0.0_dp
    real(dp) :: rss=0.0_dp, elnz2=0.0_dp, sigma2u=0.0_dp
    real(dp), allocatable :: arma_par(:), lgarch_par(:)
    real(dp), allocatable :: hessian_arma(:,:), vcov_arma(:,:), vcov_lgarch(:,:)
    real(dp), allocatable :: fitted_sd(:), log_sigma2(:), residuals(:), arma_residuals(:)
  end type lgarch_fit_result

  public :: lgarch_simulate, lgarch_arma_recursion, lgarch_objective, fit_lgarch
  public :: lgarch_is_stable
contains
  logical function lgarch_is_stable(arch,garch) result(ok)
    real(dp), intent(in) :: arch(:),garch(:)
    real(dp), allocatable :: phi(:)
    integer :: p
    p=max(size(arch),size(garch)); allocate(phi(p)); phi=0.0_dp
    if(size(arch)>0) phi(:size(arch))=phi(:size(arch))+arch
    if(size(garch)>0) phi(:size(garch))=phi(:size(garch))+garch
    ok=spectral_radius(phi)<1.0_dp
  end function lgarch_is_stable

  subroutine lgarch_simulate(n,constant,arch,garch,y,xreg,innovations,backcast_lnsigma2,backcast_lnz2, &
      log_sigma2,sigma,z,lnz2,stable)
    integer,intent(in)::n
    real(dp),intent(in)::constant,arch(:),garch(:)
    real(dp),intent(out)::y(n)
    real(dp),intent(in),optional::xreg(n),innovations(n),backcast_lnsigma2(:),backcast_lnz2(:)
    real(dp),intent(out),optional::log_sigma2(n),sigma(n),z(n),lnz2(n)
    logical,intent(out),optional::stable
    integer::p,q,m,t,j
    real(dp),allocatable::zz(:),lz(:),ls(:),phi(:),xr(:)
    real(dp)::elnz,els,xmean,innov
    p=size(arch); q=size(garch); m=max(p,q)
    if(m<1) m=1
    allocate(zz(n),lz(n+m),ls(n+m),phi(m),xr(n+m)); phi=0.0_dp; xr=0.0_dp
    if(size(arch)>0) phi(:size(arch))=phi(:size(arch))+arch
    if(size(garch)>0) phi(:size(garch))=phi(:size(garch))+garch
    if(present(innovations)) then; zz=innovations; else; call random_normal_vector(zz); end if
    elnz=mean_value(log(max(zz*zz,tiny(1.0_dp))))
    if(present(backcast_lnz2)) then
      if(size(backcast_lnz2)<m) error stop "lgarch_simulate: insufficient lnz2 backcasts"
      lz(:m)=backcast_lnz2(size(backcast_lnz2)-m+1:)
    else
      lz(:m)=elnz
    end if
    lz(m+1:)=log(max(zz*zz,tiny(1.0_dp)))
    if(present(xreg)) then; xr(m+1:)=xreg; xmean=mean_value(xreg); xr(:m)=xmean; else; xmean=0.0_dp; end if
    if(present(backcast_lnsigma2)) then
      if(size(backcast_lnsigma2)<m) error stop "lgarch_simulate: insufficient lnsigma2 backcasts"
      ls(:m)=backcast_lnsigma2(size(backcast_lnsigma2)-m+1:)
    else
      if(abs(1.0_dp-sum(phi))<=sqrt(epsilon(1.0_dp))) then
        els=0.0_dp
      else
        els=(constant+sum(arch)*elnz+xmean)/(1.0_dp-sum(phi))
      end if
      ls(:m)=els
    end if
    do t=m+1,m+n
      innov=constant+xr(t)
      do j=1,p; innov=innov+arch(j)*lz(t-j); end do
      ls(t)=innov
      do j=1,max(size(arch),size(garch)); ls(t)=ls(t)+phi(j)*ls(t-j); end do
    end do
    y=exp(0.5_dp*ls(m+1:))*zz
    if(present(log_sigma2)) log_sigma2=ls(m+1:)
    if(present(sigma)) sigma=exp(0.5_dp*ls(m+1:))
    if(present(z)) z=zz
    if(present(lnz2)) lnz2=lz(m+1:)
    if(present(stable)) stable=spectral_radius(phi)<1.0_dp
  end subroutine lgarch_simulate

  subroutine lgarch_arma_recursion(y,intercept,phi,theta,u,lny2_adjusted,xreg,xcoef,mean_correction)
    real(dp),intent(in)::y(:),intercept,phi,theta
    real(dp),intent(out)::u(size(y)),lny2_adjusted(size(y))
    real(dp),intent(in),optional::xreg(:,:),xcoef(:)
    logical,intent(in),optional::mean_correction
    integer::n,t
    real(dp)::eln,fit,imean
    real(dp),allocatable::ly(:),innov(:)
    logical::mc
    logical,allocatable::nonzero(:)
    real(dp) :: prev_u
    n=size(y); mc=.false.; if(present(mean_correction)) mc=mean_correction
    allocate(ly(n+1),innov(n+1),nonzero(n)); nonzero=(abs(y)>tiny(1.0_dp))
    if(.not.any(nonzero)) error stop "lgarch_arma_recursion: all observations are zero"
    ly(2:)=log(max(y*y,tiny(1.0_dp)))
    eln=sum(pack(ly(2:),nonzero))/real(count(nonzero),dp)
    where(.not.nonzero) ly(2:)=eln
    if(mc) ly(2:)=ly(2:)-eln
    ly(1)=merge(0.0_dp,eln,mc)
    innov(2:)=intercept
    if(present(xreg)) then
      if(.not.present(xcoef)) error stop "lgarch_arma_recursion: xcoef missing"
      if(size(xreg,1)/=n .or. size(xreg,2)/=size(xcoef)) error stop "lgarch_arma_recursion: xreg mismatch"
      innov(2:)=innov(2:)+matmul(xreg,xcoef)
    end if
    imean=mean_value(innov(2:)); innov(1)=imean
    u=0.0_dp; prev_u=0.0_dp
    do t=2,n+1
      fit=innov(t)+phi*ly(t-1)+theta*prev_u
      if(.not.nonzero(t-1)) then
        ly(t)=fit; u(t-1)=0.0_dp
      else
        u(t-1)=ly(t)-fit
      end if
      prev_u=u(t-1)
    end do
    lny2_adjusted=ly(2:)
  end subroutine lgarch_arma_recursion

  real(dp) function lgarch_objective(y,pars,arch_order,garch_order,method,xreg,mean_correction) result(value)
    real(dp),intent(in)::y(:),pars(:)
    integer,intent(in)::arch_order,garch_order,method
    real(dp),intent(in),optional::xreg(:,:)
    logical,intent(in),optional::mean_correction
    real(dp),allocatable::u(:),ly(:),xcoef(:),uu(:)
    real(dp)::intercept,phi,theta,sigma2
    logical::mc
    logical,allocatable::nz(:)
    integer::idx,k
    mc=.false.; if(present(mean_correction)) mc=mean_correction
    k=0; if(present(xreg)) k=size(xreg,2)
    idx=1; intercept=pars(idx); idx=idx+1; phi=0.0_dp; theta=0.0_dp
    if(arch_order>0) then; phi=pars(idx); idx=idx+1; end if
    if(garch_order>0) then; theta=pars(idx); idx=idx+1; end if
    allocate(xcoef(k)); if(k>0) then; xcoef=pars(idx:idx+k-1); idx=idx+k; end if
    allocate(u(size(y)),ly(size(y)),nz(size(y))); nz=(abs(y)>tiny(1.0_dp))
    if(present(xreg)) then
      call lgarch_arma_recursion(y,intercept,phi,theta,u,ly,xreg,xcoef,mc)
    else
      call lgarch_arma_recursion(y,intercept,phi,theta,u,ly,mean_correction=mc)
    end if
    uu=pack(u,nz)
    select case(method)
    case(LGARCH_LS)
      value=sum(uu*uu)
    case(LGARCH_ML)
      sigma2=pars(idx)
      if(sigma2<=0.0_dp) then; value=-huge(1.0_dp); else
        value=-0.5_dp*real(size(uu),dp)*(log(sigma2)+log(2.0_dp*pi))-0.5_dp*sum(uu*uu)/sigma2
      end if
    case(LGARCH_CEX2)
      sigma2=pars(idx)
      value=-0.5_dp*real(size(uu),dp)*log(2.0_dp*pi)+0.5_dp*sum(uu+sigma2-exp(uu+sigma2))
    case default
      error stop "lgarch_objective: unknown method"
    end select
  end function lgarch_objective

  subroutine fit_lgarch(y,arch_order,garch_order,method,result,xreg,mean_correction,compute_vcov,max_iter,tol, &
      initial_values,lower_bounds,upper_bounds,objective_penalty)
    real(dp),intent(in)::y(:)
    integer,intent(in)::arch_order,garch_order,method
    type(lgarch_fit_result),intent(out)::result
    real(dp),intent(in),optional::xreg(:,:)
    logical,intent(in),optional::mean_correction,compute_vcov
    integer,intent(in),optional::max_iter
    real(dp),intent(in),optional::tol,initial_values(:),lower_bounds(:),upper_bounds(:),objective_penalty
    integer::n,k,npar_full,npar_opt,idx,info,itmax,i
    real(dp)::eln,phi0,theta0,xmean,fmin,opttol,nanv,sdarma,penalty
    real(dp),allocatable::xmat(:,:),x0(:),lo(:),hi(:),full(:),u(:),ly(:),uu(:),h(:,:),hinv(:,:),jac(:,:)
    logical::mc,do_vcov,conv
    logical,allocatable::nz(:)
    if(arch_order<garch_order) error stop "fit_lgarch: garch order cannot exceed arch order"
    if(arch_order<0 .or. arch_order>1 .or. garch_order<0 .or. garch_order>1) &
      error stop "fit_lgarch: estimation supports orders 0 or 1, matching lgarch 0.7"
    if(method<LGARCH_LS .or. method>LGARCH_CEX2) error stop "fit_lgarch: invalid method"
    n=size(y); k=0; if(present(xreg)) k=size(xreg,2)
    allocate(xmat(n,k)); if(k>0) then
      if(size(xreg,1)/=n) error stop "fit_lgarch: xreg row mismatch"
      xmat=xreg
    end if
    allocate(nz(n)); nz=(abs(y)>tiny(1.0_dp)); if(.not.any(nz)) error stop "fit_lgarch: all observations are zero"
    eln=sum(log(pack(y*y,nz)))/real(count(nz),dp)
    mc=.false.; if(present(mean_correction)) mc=mean_correction
    if(method==LGARCH_CEX2) mc=.true.
    do_vcov=.true.; if(present(compute_vcov)) do_vcov=compute_vcov
    npar_full=1+arch_order+garch_order+k+merge(0,1,method==LGARCH_LS)
    npar_opt=npar_full-merge(1,0,mc)
    allocate(x0(npar_opt),lo(npar_opt),hi(npar_opt),full(npar_full))
    full=0.0_dp; idx=1
    phi0=merge(0.9_dp,0.0_dp,arch_order>0); theta0=merge(-0.8_dp,0.0_dp,garch_order>0)
    xmean=0.0_dp
    full(idx)=(1.0_dp-phi0)*eln-xmean; idx=idx+1
    if(arch_order>0) then; full(idx)=phi0; idx=idx+1; end if
    if(garch_order>0) then; full(idx)=theta0; idx=idx+1; end if
    if(k>0) then; full(idx:idx+k-1)=0.01_dp; xmean=mean_value(matmul(xmat,full(idx:idx+k-1))); full(1)=(1.0_dp-phi0)*eln-xmean; idx=idx+k; end if
    if(method==LGARCH_ML) full(idx)=4.94_dp
    if(method==LGARCH_CEX2) full(idx)=-1.27_dp
    if(mc) then; x0=full(2:); else; x0=full; end if
    lo=-1.0e6_dp; hi=1.0e6_dp; idx=merge(0,1,mc)+1
    if(arch_order>0) then; lo(idx)=-0.999999_dp; hi(idx)=0.999999_dp; idx=idx+1; end if
    if(garch_order>0) then; lo(idx)=-0.999999_dp; hi(idx)=0.999999_dp; idx=idx+1; end if
    idx=idx+k
    if(method==LGARCH_ML) then; lo(idx)=1.0e-8_dp; hi(idx)=1.0e4_dp; end if
    if(method==LGARCH_CEX2) then; lo(idx)=-50.0_dp; hi(idx)=50.0_dp; end if
    if(present(initial_values)) then
      if(size(initial_values)==npar_full) then
        full=initial_values
        if(mc) then; x0=full(2:); else; x0=full; end if
      else if(size(initial_values)==npar_opt) then
        x0=initial_values
      else
        error stop "fit_lgarch: initial_values has wrong size"
      end if
    end if
    if(present(lower_bounds)) call override_bounds(lower_bounds,lo,"lower_bounds")
    if(present(upper_bounds)) call override_bounds(upper_bounds,hi,"upper_bounds")
    if(any(lo>=hi)) error stop "fit_lgarch: lower bound must be below upper bound"
    x0=min(max(x0,lo),hi)
    penalty=1.0e30_dp; if(present(objective_penalty)) penalty=objective_penalty
    itmax=3000; if(present(max_iter)) itmax=max_iter
    opttol=1.0e-8_dp; if(present(tol)) opttol=tol
    call minimize_nelder_mead(obj,x0,lo,hi,fmin,conv,result%iterations,itmax,opttol)
    if(mc) then
      full(2:)=x0; full(1)=0.0_dp
      if(arch_order>0) full(1)=(1.0_dp-full(2))*eln
    else
      full=x0
    end if
    result%arch_order=arch_order; result%garch_order=garch_order; result%method=method
    result%mean_correction=mc; result%converged=conv
    allocate(result%arma_par(npar_full)); result%arma_par=full
    allocate(u(n),ly(n)); call recurse_full(full,u,ly)
    uu=pack(u,nz); result%arma_residuals=u; result%rss=sum(uu*uu); result%sigma2u=sample_variance(uu)
    select case(method)
    case(LGARCH_LS); result%objective_arma=fmin; result%elnz2=-safe_log_mean_exp(uu)
    case(LGARCH_ML); result%objective_arma=-fmin; result%elnz2=-safe_log_mean_exp(uu-mean_value(uu)); result%sigma2u=full(npar_full)
    case(LGARCH_CEX2); result%objective_arma=-fmin; result%elnz2=full(npar_full)
    end select
    call make_lgarch_parameters(full,result%elnz2,result%lgarch_par)
    allocate(result%fitted_sd(n),result%log_sigma2(n),result%residuals(n))
    result%log_sigma2=ly-(u+result%elnz2)
    if(mc) result%log_sigma2=result%log_sigma2+eln
    result%fitted_sd=exp(0.5_dp*result%log_sigma2); result%residuals=y/result%fitted_sd
    result%loglik_model=normal_model_loglik(y,result%fitted_sd,nz)
    if(method==LGARCH_LS) then
      sdarma=sqrt(max(sample_variance(uu),tiny(1.0_dp)))
      result%loglik_arma=-0.5_dp*real(size(uu),dp)*log(2.0_dp*pi*sdarma*sdarma)-0.5_dp*sum((uu/sdarma)**2)
    else
      result%loglik_arma=result%objective_arma
    end if
    nanv=ieee_value(0.0_dp,ieee_quiet_nan)
    if(do_vcov) then
      allocate(h(npar_opt,npar_opt),hinv(npar_opt,npar_opt)); call numerical_hessian(obj,x0,h)
      call inverse_matrix(h,hinv,info)
      if(info/=0) hinv=nanv
      if(method==LGARCH_LS .and. info==0) hinv=2.0_dp*result%sigma2u*hinv
      allocate(result%hessian_arma(npar_full,npar_full),result%vcov_arma(npar_full,npar_full)); result%hessian_arma=nanv; result%vcov_arma=nanv
      if(mc) then; result%hessian_arma(2:,2:)=h; result%vcov_arma(2:,2:)=hinv; else; result%hessian_arma=h; result%vcov_arma=hinv; end if
      allocate(result%vcov_lgarch(size(result%lgarch_par),size(result%lgarch_par))); result%vcov_lgarch=nanv
      allocate(jac(size(result%lgarch_par)-1,npar_full)); jac=0.0_dp
      call dynamic_jacobian(jac)
      if(all(ieee_is_finite(result%vcov_arma))) then
        result%vcov_lgarch(:size(jac,1),:size(jac,1))=matmul(jac,matmul(result%vcov_arma,transpose(jac)))
      end if
      if(method==LGARCH_CEX2 .and. ieee_is_finite(result%vcov_arma(npar_full,npar_full))) then
        result%vcov_lgarch(size(result%lgarch_par),size(result%lgarch_par))=result%vcov_arma(npar_full,npar_full)
      else
        result%vcov_lgarch(size(result%lgarch_par),size(result%lgarch_par))=sample_variance(pack(result%residuals**2-log(max(result%residuals**2,tiny(1.0_dp))),nz))/real(count(nz),dp)
      end if
    end if
  contains
    real(dp) function obj(v) result(f)
      real(dp),intent(in)::v(:)
      real(dp)::p(npar_full),raw
      if(mc) then; p(1)=0.0_dp; p(2:)=v; else; p=v; end if
      raw=lgarch_objective(y,p,arch_order,garch_order,method,xmat,mc)
      if(.not.ieee_is_finite(raw)) then; f=penalty
      else if(method==LGARCH_LS) then; f=raw
      else; f=-raw
      end if
    end function obj
    subroutine recurse_full(pv,ur,lr)
      real(dp),intent(in)::pv(:); real(dp),intent(out)::ur(:),lr(:)
      real(dp)::ph,th
      real(dp),allocatable::xc(:)
      integer::ii
      ii=2; ph=0.0_dp; th=0.0_dp
      if(arch_order>0) then; ph=pv(ii); ii=ii+1; end if
      if(garch_order>0) then; th=pv(ii); ii=ii+1; end if
      allocate(xc(k)); if(k>0) xc=pv(ii:ii+k-1)
      if(k>0) then; call lgarch_arma_recursion(y,merge(0.0_dp,pv(1),mc),ph,th,ur,lr,xmat,xc,mc)
      else; call lgarch_arma_recursion(y,merge(0.0_dp,pv(1),mc),ph,th,ur,lr,mean_correction=mc); end if
    end subroutine recurse_full
    subroutine override_bounds(source,target,label)
      real(dp),intent(in)::source(:)
      real(dp),intent(inout)::target(:)
      character(len=*),intent(in)::label
      if(size(source)==npar_full) then
        if(mc) then; target=source(2:); else; target=source; end if
      else if(size(source)==npar_opt) then
        target=source
      else
        error stop "fit_lgarch: "//label//" has wrong size"
      end if
    end subroutine override_bounds
    subroutine make_lgarch_parameters(pv,e,out)
      real(dp),intent(in)::pv(:),e; real(dp),allocatable,intent(out)::out(:)
      real(dp)::ph,th
      integer::ii,oo
      allocate(out(1+arch_order+garch_order+k+1)); ii=2; ph=0.0_dp; th=0.0_dp
      if(arch_order>0) then; ph=pv(ii); ii=ii+1; end if
      if(garch_order>0) then; th=pv(ii); ii=ii+1; end if
      out(1)=pv(1)-(1.0_dp+th)*e; oo=2
      if(arch_order>0) then; out(oo)=ph+th; oo=oo+1; end if
      if(garch_order>0) then; out(oo)=-th; oo=oo+1; end if
      if(k>0) then; out(oo:oo+k-1)=pv(ii:ii+k-1); oo=oo+k; end if
      out(oo)=e
    end subroutine make_lgarch_parameters
    subroutine dynamic_jacobian(jc)
      real(dp),intent(out)::jc(:,:)
      integer::ia,im,ix,ro
      jc=0.0_dp; jc(1,1)=1.0_dp; ia=0; im=0
      if(arch_order>0) ia=2
      if(garch_order>0) im=2+arch_order
      ro=2
      if(arch_order>0) then; jc(ro,ia)=1.0_dp; if(garch_order>0) jc(ro,im)=1.0_dp; ro=ro+1; end if
      if(garch_order>0) then; jc(ro,im)=-1.0_dp; ro=ro+1; end if
      ix=2+arch_order+garch_order
      do i=1,k; jc(ro,ix+i-1)=1.0_dp; ro=ro+1; end do
    end subroutine dynamic_jacobian
  end subroutine fit_lgarch

  pure real(dp) function normal_model_loglik(y,sd,nz) result(ll)
    real(dp),intent(in)::y(:),sd(:); logical,intent(in)::nz(:)
    integer::i
    ll=0.0_dp
    do i=1,size(y)
      if(nz(i)) ll=ll-0.5_dp*(log(2.0_dp*pi)+2.0_dp*log(sd(i))+(y(i)/sd(i))**2)
    end do
  end function normal_model_loglik
end module lgarch_univariate
