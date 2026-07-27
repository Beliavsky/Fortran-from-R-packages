! SPDX-License-Identifier: GPL-2.0-only
! Copyright (C) 2014-2025 Genaro Sucarrat (original R/C package)
! Copyright (C) 2026 Modern Fortran translation contributors
! This file is distributed under the GNU General Public License version 2 only.
module lgarch_multivariate
  use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan, ieee_is_finite
  use lgarch_kinds, only : dp, pi
  use lgarch_rng, only : random_normal_vector
  use lgarch_utils, only : mean_value, sample_variance, safe_log_mean_exp
  use lgarch_linalg, only : solve_linear, inverse_matrix, cholesky_lower, logdet_spd, spectral_radius, correlation_matrix, covariance_matrix
  use lgarch_optimizer, only : minimize_nelder_mead, numerical_hessian
  implicit none
  private
  type, public :: mlgarch_fit_result
    integer :: dimension=0, arch_order=0, garch_order=0, iterations=0
    logical :: converged=.false.
    real(dp) :: objective_varma=0.0_dp, loglik_model=0.0_dp
    real(dp), allocatable :: varma_par(:), mlgarch_par(:), elnz2(:)
    real(dp), allocatable :: innovation_cov(:,:), hessian_varma(:,:), vcov_varma(:,:), vcov_mlgarch(:,:)
    real(dp), allocatable :: fitted_sd(:,:), log_sigma2(:,:), residuals(:,:), varma_residuals(:,:)
  end type mlgarch_fit_result
  public :: rmnorm, mlgarch_simulate, mlgarch_varma_recursion, mlgarch_objective, fit_mlgarch
  public :: mlgarch_is_stable
contains
  subroutine rmnorm(n,mean,vcov,x,info)
    integer,intent(in)::n
    real(dp),intent(in)::mean(:),vcov(:,:)
    real(dp),intent(out)::x(n,size(mean))
    integer,intent(out),optional::info
    real(dp),allocatable::l(:,:),z(:)
    integer::j,istat,m
    m=size(mean); allocate(l(m,m),z(m)); call cholesky_lower(vcov,l,istat)
    if(present(info)) info=istat
    if(istat/=0) then; x=0.0_dp; return; end if
    do j=1,n
      call random_normal_vector(z); x(j,:)=mean+matmul(l,z)
    end do
  end subroutine rmnorm

  logical function mlgarch_is_stable(arch,garch) result(ok)
    real(dp),intent(in)::arch(:,:),garch(:,:)
    real(dp),allocatable::a(:,:),wr(:),wi(:),vl(:,:),vr(:,:),work(:)
    integer::m,lwork,info
    interface
      subroutine dgeev(jobvl,jobvr,n,a,lda,wr,wi,vl,ldvl,vr,ldvr,work,lwork,info)
        import dp
        character(len=1),intent(in)::jobvl,jobvr
        integer,intent(in)::n,lda,ldvl,ldvr,lwork
        integer,intent(out)::info
        real(dp),intent(inout)::a(lda,*),work(*)
        real(dp),intent(out)::wr(*),wi(*),vl(ldvl,*),vr(ldvr,*)
      end subroutine dgeev
    end interface
    m=size(arch,1); if(any(shape(arch)/=[m,m]) .or. any(shape(garch)/=[m,m])) error stop "mlgarch_is_stable: size mismatch"
    allocate(a(m,m),wr(m),wi(m),vl(1,1),vr(1,1)); a=arch+garch; lwork=max(1,4*m); allocate(work(lwork))
    call dgeev('N','N',m,a,m,wr,wi,vl,1,vr,1,work,lwork,info)
    ok=info==0 .and. maxval(sqrt(wr*wr+wi*wi))<1.0_dp
  end function mlgarch_is_stable

  subroutine mlgarch_simulate(n,constant,arch,garch,y,xreg,innovations,innovations_vcov,backcast_lnsigma2,backcast_lnz2, &
      log_sigma2,sigma,z,stable)
    integer,intent(in)::n
    real(dp),intent(in)::constant(:),arch(:,:),garch(:,:)
    real(dp),intent(out)::y(n,size(constant))
    real(dp),intent(in),optional::xreg(:,:),innovations(:,:),innovations_vcov(:,:),backcast_lnsigma2(:),backcast_lnz2(:)
    real(dp),intent(out),optional::log_sigma2(n,size(constant)),sigma(n,size(constant)),z(n,size(constant))
    logical,intent(out),optional::stable
    integer::m,t,info
    real(dp),allocatable::zz(:,:),lz(:,:),ls(:,:),phi(:,:),cov(:,:),rhs(:),init(:),xr(:,:)
    real(dp)::eln(size(constant)),xmean(size(constant))
    m=size(constant)
    if(any(shape(arch)/=[m,m]) .or. any(shape(garch)/=[m,m])) error stop "mlgarch_simulate: coefficient size mismatch"
    allocate(zz(n,m),lz(n+1,m),ls(n+1,m),phi(m,m),cov(m,m),rhs(m),init(m),xr(n+1,m)); phi=arch+garch; xr=0.0_dp
    if(present(innovations)) then
      if(any(shape(innovations)/=[n,m])) error stop "mlgarch_simulate: innovation size mismatch"
      zz=innovations
    else
      cov=0.0_dp
      if(present(innovations_vcov)) then; cov=innovations_vcov; else; do t=1,m; cov(t,t)=1.0_dp; end do; end if
      call rmnorm(n,spread(0.0_dp,1,m),cov,zz,info); if(info/=0) error stop "mlgarch_simulate: covariance is not positive definite"
    end if
    do t=1,m; eln(t)=mean_value(log(max(zz(:,t)*zz(:,t),tiny(1.0_dp)))); end do
    if(present(backcast_lnz2)) then; if(size(backcast_lnz2)/=m) error stop "mlgarch_simulate: lnz2 backcast mismatch"; lz(1,:)=backcast_lnz2
    else; lz(1,:)=eln; end if
    lz(2:,:)=log(max(zz*zz,tiny(1.0_dp)))
    xmean=0.0_dp
    if(present(xreg)) then
      if(any(shape(xreg)/=[n,m])) error stop "mlgarch_simulate: xreg must be n by dimension"
      xr(2:,:)=xreg
      do t=1,m; xmean(t)=mean_value(xreg(:,t)); end do
      xr(1,:)=xmean
    end if
    if(present(backcast_lnsigma2)) then
      if(size(backcast_lnsigma2)/=m) error stop "mlgarch_simulate: lnsigma2 backcast mismatch"
      ls(1,:)=backcast_lnsigma2
    else
      rhs=constant+xmean+matmul(arch,eln)
      call solve_linear(identity(m)-phi,rhs,init,info); if(info/=0) error stop "mlgarch_simulate: unstable or singular initialization"
      ls(1,:)=init
    end if
    do t=2,n+1
      ls(t,:)=constant+xr(t,:)+matmul(phi,ls(t-1,:))+matmul(arch,lz(t-1,:))
    end do
    y=exp(0.5_dp*ls(2:,:))*zz
    if(present(log_sigma2)) log_sigma2=ls(2:,:)
    if(present(sigma)) sigma=exp(0.5_dp*ls(2:,:))
    if(present(z)) z=zz
    if(present(stable)) stable=mlgarch_is_stable(arch,garch)
  contains
    pure function identity(nm) result(a)
      integer,intent(in)::nm; real(dp)::a(nm,nm); integer::ii
      a=0.0_dp; do ii=1,nm; a(ii,ii)=1.0_dp; end do
    end function identity
  end subroutine mlgarch_simulate

  subroutine mlgarch_varma_recursion(y,constant,phi,theta,u,lny2_adjusted,xreg,xcoef)
    real(dp),intent(in)::y(:,:),constant(:),phi(:,:),theta(:,:)
    real(dp),intent(out)::u(size(y,1),size(y,2)),lny2_adjusted(size(y,1),size(y,2))
    real(dp),intent(in),optional::xreg(:,:),xcoef(:,:)
    integer::n,m,t,j
    real(dp),allocatable::ly(:,:),innov(:,:)
    logical,allocatable::nz(:,:)
    real(dp)::means(size(y,2)),fit(size(y,2)),prev_u(size(y,2))
    n=size(y,1); m=size(y,2)
    if(size(constant)/=m .or. any(shape(phi)/=[m,m]) .or. any(shape(theta)/=[m,m])) error stop "mlgarch_varma_recursion: size mismatch"
    allocate(ly(n+1,m),innov(n+1,m),nz(n,m)); nz=(abs(y)>tiny(1.0_dp)); ly(2:,:)=log(max(y*y,tiny(1.0_dp)))
    do j=1,m
      if(.not.any(nz(:,j))) error stop "mlgarch_varma_recursion: an entire series is zero"
      means(j)=sum(pack(ly(2:,j),nz(:,j)))/real(count(nz(:,j)),dp)
      where(.not.nz(:,j)) ly(2:,j)=means(j)
    end do
    ly(1,:)=means; innov(2:,:)=spread(constant,1,n)
    if(present(xreg)) then
      if(.not.present(xcoef)) error stop "mlgarch_varma_recursion: xcoef missing"
      if(size(xreg,1)/=n .or. size(xcoef,1)/=m .or. size(xreg,2)/=size(xcoef,2)) error stop "mlgarch_varma_recursion: xreg mismatch"
      innov(2:,:)=innov(2:,:)+matmul(xreg,transpose(xcoef))
    end if
    innov(1,:)=sum(innov(2:,:),dim=1)/real(n,dp); u=0.0_dp; prev_u=0.0_dp
    do t=2,n+1
      fit=innov(t,:)+matmul(phi,ly(t-1,:))+matmul(theta,prev_u)
      do j=1,m
        if(.not.nz(t-1,j)) then; ly(t,j)=fit(j); u(t-1,j)=0.0_dp; else; u(t-1,j)=ly(t,j)-fit(j); end if
      end do
      prev_u=u(t-1,:)
    end do
    lny2_adjusted=ly(2:,:)
  end subroutine mlgarch_varma_recursion

  real(dp) function mlgarch_objective(y,pars,arch_order,garch_order,xreg) result(ll)
    real(dp),intent(in)::y(:,:),pars(:)
    integer,intent(in)::arch_order,garch_order
    real(dp),intent(in),optional::xreg(:,:)
    integer::n,m,k,idx,i,j,info,nvalid
    real(dp),allocatable::constant(:),phi(:,:),theta(:,:),xcoef(:,:),s(:,:),sinv(:,:),u(:,:),ly(:,:),uv(:,:)
    logical,allocatable::valid(:)
    real(dp)::logdet,quad
    n=size(y,1); m=size(y,2); k=0; if(present(xreg)) k=size(xreg,2)
    allocate(constant(m),phi(m,m),theta(m,m),xcoef(m,k),s(m,m),sinv(m,m),u(n,m),ly(n,m),valid(n))
    idx=1; constant=pars(idx:idx+m-1); idx=idx+m; phi=0.0_dp; theta=0.0_dp
    if(arch_order>0) then; phi=reshape(pars(idx:idx+m*m-1),[m,m]); idx=idx+m*m; end if
    if(garch_order>0) then; theta=reshape(pars(idx:idx+m*m-1),[m,m]); idx=idx+m*m; end if
    if(k>0) then; xcoef=reshape(pars(idx:idx+m*k-1),[m,k]); idx=idx+m*k; end if
    s=0.0_dp
    do i=1,m; s(i,i)=pars(idx); idx=idx+1; end do
    do j=1,m-1; do i=j+1,m; s(i,j)=pars(idx); s(j,i)=pars(idx); idx=idx+1; end do; end do
    call logdet_spd(s,logdet,info); if(info/=0) then; ll=-huge(1.0_dp); return; end if
    call inverse_matrix(s,sinv,info); if(info/=0) then; ll=-huge(1.0_dp); return; end if
    if(k>0) then; call mlgarch_varma_recursion(y,constant,phi,theta,u,ly,xreg,xcoef)
    else; call mlgarch_varma_recursion(y,constant,phi,theta,u,ly); end if
    valid=.not.any(abs(y)<=tiny(1.0_dp),dim=2); nvalid=count(valid); if(nvalid==0) then; ll=-huge(1.0_dp); return; end if
    uv=pack_rows(u,valid); quad=0.0_dp
    do i=1,nvalid; quad=quad+dot_product(uv(i,:),matmul(sinv,uv(i,:))); end do
    ll=-0.5_dp*real(nvalid*m,dp)*log(2.0_dp*pi)-0.5_dp*real(nvalid,dp)*logdet-0.5_dp*quad
  contains
    function pack_rows(a,mask) result(b)
      real(dp),intent(in)::a(:,:); logical,intent(in)::mask(:)
      real(dp),allocatable::b(:,:); integer::ii,jj
      allocate(b(count(mask),size(a,2))); jj=0
      do ii=1,size(mask); if(mask(ii)) then; jj=jj+1; b(jj,:)=a(ii,:); end if; end do
    end function pack_rows
  end function mlgarch_objective

  subroutine fit_mlgarch(y,arch_order,garch_order,result,xreg,compute_vcov,max_iter,tol, &
      initial_values,lower_bounds,upper_bounds,objective_penalty)
    real(dp),intent(in)::y(:,:)
    integer,intent(in)::arch_order,garch_order
    type(mlgarch_fit_result),intent(out)::result
    real(dp),intent(in),optional::xreg(:,:)
    logical,intent(in),optional::compute_vcov
    integer,intent(in),optional::max_iter
    real(dp),intent(in),optional::tol,initial_values(:),lower_bounds(:),upper_bounds(:),objective_penalty
    integer::n,m,k,npar,idx,i,j,info,itmax,nmlg
    real(dp)::fmin,opttol,nanv,epsj,penalty
    real(dp),allocatable::xmat(:,:),x0(:),lo(:),hi(:),u(:,:),ly(:,:),h(:,:),covp(:,:),mapplus(:),mapminus(:),jac(:,:)
    real(dp),allocatable::phi(:,:),theta(:,:),xcoef(:,:),s(:,:),means(:),fullcov(:,:)
    logical::conv,do_vcov
    logical,allocatable::nz(:,:),valid(:)
    real(dp)::term,exmean
    if(arch_order<garch_order) error stop "fit_mlgarch: garch order cannot exceed arch order"
    if(arch_order<0 .or. arch_order>1 .or. garch_order<0 .or. garch_order>1) error stop "fit_mlgarch: orders must be 0 or 1"
    n=size(y,1); m=size(y,2); k=0; if(present(xreg)) k=size(xreg,2)
    allocate(xmat(n,k)); if(k>0) then; if(size(xreg,1)/=n) error stop "fit_mlgarch: xreg row mismatch"; xmat=xreg; end if
    allocate(nz(n,m),valid(n),means(m)); nz=(abs(y)>tiny(1.0_dp)); valid=.not.any(.not.nz,dim=2)
    if(count(valid)==0) error stop "fit_mlgarch: no rows without zeros"
    do j=1,m; if(.not.any(nz(:,j))) error stop "fit_mlgarch: an entire series is zero"; means(j)=sum(log(pack(y(:,j)**2,nz(:,j))))/real(count(nz(:,j)),dp); end do
    npar=m+(arch_order+garch_order)*m*m+k*m+m+m*(m-1)/2
    allocate(x0(npar),lo(npar),hi(npar)); x0=0.0_dp; lo=-1.0e6_dp; hi=1.0e6_dp; idx=1
    x0(idx:idx+m-1)=0.1_dp*means; idx=idx+m
    if(arch_order>0) then
      do i=1,m; x0(idx+(i-1)*m+i-1)=0.9_dp; end do
      lo(idx:idx+m*m-1)=-0.999999_dp; hi(idx:idx+m*m-1)=0.999999_dp; idx=idx+m*m
    end if
    if(garch_order>0) then
      do i=1,m; x0(idx+(i-1)*m+i-1)=-0.8_dp; end do
      lo(idx:idx+m*m-1)=-0.999999_dp; hi(idx:idx+m*m-1)=0.999999_dp; idx=idx+m*m
    end if
    if(k>0) then; x0(idx:idx+m*k-1)=0.01_dp; idx=idx+m*k; end if
    do i=1,m; x0(idx)=4.94_dp; lo(idx)=1.0e-8_dp; hi(idx)=1.0e4_dp; idx=idx+1; end do
    fullcov=covariance_matrix(log(max(y*y,tiny(1.0_dp))))
    do j=1,m-1; do i=j+1,m; x0(idx)=fullcov(i,j); idx=idx+1; end do; end do
    if(present(initial_values)) then
      if(size(initial_values)/=npar) error stop "fit_mlgarch: initial_values has wrong size"
      x0=initial_values
    end if
    if(present(lower_bounds)) then
      if(size(lower_bounds)/=npar) error stop "fit_mlgarch: lower_bounds has wrong size"
      lo=lower_bounds
    end if
    if(present(upper_bounds)) then
      if(size(upper_bounds)/=npar) error stop "fit_mlgarch: upper_bounds has wrong size"
      hi=upper_bounds
    end if
    if(any(lo>=hi)) error stop "fit_mlgarch: lower bound must be below upper bound"
    x0=min(max(x0,lo),hi)
    penalty=1.0e30_dp; if(present(objective_penalty)) penalty=objective_penalty
    itmax=5000; if(present(max_iter)) itmax=max_iter; opttol=1.0e-7_dp; if(present(tol)) opttol=tol
    call minimize_nelder_mead(obj,x0,lo,hi,fmin,conv,result%iterations,itmax,opttol)
    result%dimension=m; result%arch_order=arch_order; result%garch_order=garch_order; result%converged=conv
    result%objective_varma=-fmin; result%varma_par=x0
    allocate(phi(m,m),theta(m,m),xcoef(m,k),s(m,m),u(n,m),ly(n,m)); call unpack(x0,phi,theta,xcoef,s)
    if(k>0) then; call mlgarch_varma_recursion(y,x0(:m),phi,theta,u,ly,xmat,xcoef)
    else; call mlgarch_varma_recursion(y,x0(:m),phi,theta,u,ly); end if
    result%varma_residuals=u; result%innovation_cov=s; allocate(result%elnz2(m))
    do j=1,m
      result%elnz2(j)=-safe_log_mean_exp(pack(u(:,j)-sum(pack(u(:,j),nz(:,j)))/real(count(nz(:,j)),dp),nz(:,j)))
    end do
    call map_parameters(x0,result%elnz2,result%mlgarch_par)
    allocate(result%fitted_sd(n,m),result%log_sigma2(n,m),result%residuals(n,m))
    do j=1,m
      result%log_sigma2(:,j)=ly(:,j)-(u(:,j)+merge(result%elnz2(j),0.0_dp,nz(:,j)))
    end do
    result%fitted_sd=exp(0.5_dp*result%log_sigma2); result%residuals=y/result%fitted_sd
    call model_loglik(result%loglik_model)
    do_vcov=.true.; if(present(compute_vcov)) do_vcov=compute_vcov; nanv=ieee_value(0.0_dp,ieee_quiet_nan)
    if(do_vcov) then
      allocate(h(npar,npar),covp(npar,npar)); call numerical_hessian(obj,x0,h); call inverse_matrix(h,covp,info); if(info/=0) covp=nanv
      result%hessian_varma=h; result%vcov_varma=covp; nmlg=size(result%mlgarch_par)
      allocate(result%vcov_mlgarch(nmlg,nmlg)); result%vcov_mlgarch=nanv
      allocate(jac(nmlg-m,npar),mapplus(nmlg),mapminus(nmlg)); jac=0.0_dp
      do i=1,npar
        epsj=epsilon(1.0_dp)**0.25_dp*max(1.0_dp,abs(x0(i)))
        x0(i)=x0(i)+epsj; call map_parameters(x0,result%elnz2,mapplus)
        x0(i)=x0(i)-2.0_dp*epsj; call map_parameters(x0,result%elnz2,mapminus)
        x0(i)=x0(i)+epsj; jac(:,i)=(mapplus(:nmlg-m)-mapminus(:nmlg-m))/(2.0_dp*epsj)
      end do
      if(all(ieee_is_finite(covp))) result%vcov_mlgarch(:nmlg-m,:nmlg-m)=matmul(jac,matmul(covp,transpose(jac)))
      do j=1,m
        exmean=mean_value(pack(exp(u(:,j)),nz(:,j)))
        term=sample_variance(pack(exp(u(:,j)),nz(:,j)))/(exmean*exmean)+sample_variance(pack(u(:,j),nz(:,j))) &
          -2.0_dp*mean_value(pack(u(:,j)*exp(u(:,j)),nz(:,j)))/exmean
        result%vcov_mlgarch(nmlg-m+j,nmlg-m+j)=term/real(count(nz(:,j)),dp)
      end do
    end if
  contains
    real(dp) function obj(v) result(f)
      real(dp),intent(in)::v(:); real(dp)::raw
      raw=mlgarch_objective(y,v,arch_order,garch_order,xmat)
      if(ieee_is_finite(raw)) then; f=-raw; else; f=penalty; end if
    end function obj
    subroutine unpack(v,ph,th,xc,sc)
      real(dp),intent(in)::v(:); real(dp),intent(out)::ph(:,:),th(:,:),xc(:,:),sc(:,:)
      integer::ii,a,b
      ii=m+1; ph=0.0_dp; th=0.0_dp
      if(arch_order>0) then; ph=reshape(v(ii:ii+m*m-1),[m,m]); ii=ii+m*m; end if
      if(garch_order>0) then; th=reshape(v(ii:ii+m*m-1),[m,m]); ii=ii+m*m; end if
      if(k>0) then; xc=reshape(v(ii:ii+m*k-1),[m,k]); ii=ii+m*k; end if
      sc=0.0_dp; do a=1,m; sc(a,a)=v(ii); ii=ii+1; end do
      do b=1,m-1; do a=b+1,m; sc(a,b)=v(ii); sc(b,a)=v(ii); ii=ii+1; end do; end do
    end subroutine unpack
    subroutine map_parameters(v,e,out)
      real(dp),intent(in)::v(:),e(:); real(dp),allocatable,intent(out)::out(:)
      real(dp),allocatable::ph(:,:),th(:,:),xc(:,:),sc(:,:)
      integer::ii,oo
      allocate(ph(m,m),th(m,m),xc(m,k),sc(m,m)); call unpack(v,ph,th,xc,sc)
      allocate(out(m+(arch_order+garch_order)*m*m+k*m+m)); oo=1
      out(oo:oo+m-1)=v(:m)-matmul(identity(m)+th,e); oo=oo+m
      if(arch_order>0) then; out(oo:oo+m*m-1)=reshape(ph+th,[m*m]); oo=oo+m*m; end if
      if(garch_order>0) then; out(oo:oo+m*m-1)=reshape(-th,[m*m]); oo=oo+m*m; end if
      if(k>0) then; ii=m+(arch_order+garch_order)*m*m+1; out(oo:oo+m*k-1)=v(ii:ii+m*k-1); oo=oo+m*k; end if
      out(oo:)=e
    end subroutine map_parameters
    subroutine model_loglik(ll)
      real(dp),intent(out)::ll
      real(dp),allocatable::zv(:,:),sv(:,:),rr(:,:),ri(:,:)
      integer::nv,ii,istat
      real(dp)::ld,q
      zv=pack_rows(result%residuals,valid); sv=pack_rows(result%fitted_sd,valid); nv=size(zv,1); rr=correlation_matrix(zv)
      allocate(ri(m,m)); call logdet_spd(rr,ld,istat); if(istat/=0) then; ll=-huge(1.0_dp); return; end if
      call inverse_matrix(rr,ri,istat); if(istat/=0) then; ll=-huge(1.0_dp); return; end if
      q=0.0_dp; do ii=1,nv; q=q+dot_product(zv(ii,:),matmul(ri,zv(ii,:))); end do
      ll=-0.5_dp*real(nv,dp)*(real(m,dp)*log(2.0_dp*pi)+ld)-sum(log(sv))-0.5_dp*q
    end subroutine model_loglik
    function pack_rows(a,mask) result(b)
      real(dp),intent(in)::a(:,:); logical,intent(in)::mask(:)
      real(dp),allocatable::b(:,:); integer::ii,jj
      allocate(b(count(mask),size(a,2))); jj=0
      do ii=1,size(mask); if(mask(ii)) then; jj=jj+1; b(jj,:)=a(ii,:); end if; end do
    end function pack_rows
    pure function identity(nm) result(a)
      integer,intent(in)::nm; real(dp)::a(nm,nm); integer::ii
      a=0.0_dp; do ii=1,nm; a(ii,ii)=1.0_dp; end do
    end function identity
  end subroutine fit_mlgarch
end module lgarch_multivariate
