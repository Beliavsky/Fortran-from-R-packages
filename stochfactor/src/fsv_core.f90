! SPDX-License-Identifier: GPL-2.0-or-later
! Copyright (C) 2013-2026 the original stochvol and factorstochvol authors
! and the modern Fortran translation contributors.
! This file is free software: you may redistribute it and/or modify it
! under the terms of the GNU General Public License as published by the
! Free Software Foundation, either version 2, or (at your option) any later version.
! This file is distributed without any warranty; see LICENSE for details.
module fsv_core
  use sv_kinds, only : dp,tiny_var,log2pi
  use sv_rng, only : randn,rand_student_t,rand_gamma,rand_gig_slice
  use sv_linalg, only : covariance_matrix,correlation_from_covariance,symmetric_eigen,solve_spd,inverse_spd,mvn_draw,mvn_logpdf
  use sv_types
  use sv_core, only : sv_single_sweep
  use fsv_types
  implicit none
  private
  public :: simulate_fsv,fsv_covariance_path,fsv_correlation_path,covelement,corelement
  public :: make_dense_loadings,make_sparse_loadings,default_fsv_parameters
  public :: ledermann,expweightcov,static_factor_initialize,preorder,findrestrict
  public :: fit_fsv,predict_fsv,predloglik_fsv,predloglik_fsv_woodbury,aggregate_loglik_draws,dmvnorm_columns
  public :: predict_conditional_fsv,woodbury_precision,eigen_loading_diagnostics
  public :: sign_identify,order_identify,running_covariance,running_correlation
contains
  integer function ledermann(m) result(r)
    integer,intent(in)::m
    r=int(floor((2.0_dp*m+1.0_dp)/2.0_dp-sqrt((2.0_dp*m+1.0_dp)**2/4.0_dp-m*m+m)))
  end function ledermann


  subroutine make_dense_loadings(m,r,loadings)
    integer,intent(in)::m,r
    real(dp),intent(out)::loadings(m,r)
    integer::i,j
    loadings=0.0_dp
    if(r>=1)then
      loadings(1,1)=1.0_dp
      do i=2,m;loadings(i,1)=0.9_dp-0.8_dp*real(i-2,dp)/real(max(1,m-2),dp);end do
    end if
    if(r>=2.and.m>=2)then
      loadings(2,2)=1.0_dp
      do i=3,m;loadings(i,2)=0.1_dp+0.7_dp*real(i-3,dp)/real(max(1,m-3),dp);end do
    end if
    if(r>=3.and.m>=3)then
      loadings(3,3)=1.0_dp
      do i=4,m;loadings(i,3)=0.7_dp-0.3_dp*real(i-4,dp)/real(max(1,m-4),dp);end do
    end if
    do j=4,r
      if(j<=m)loadings(j,j)=1.0_dp
      do i=j+1,m;loadings(i,j)=-0.2_dp+randu_local();end do
    end do
  end subroutine make_dense_loadings

  subroutine make_sparse_loadings(m,r,cutoff,loadings)
    integer,intent(in)::m,r
    real(dp),intent(in)::cutoff
    real(dp),intent(out)::loadings(m,r)
    integer::i,j,target,nactive
    loadings=0.0_dp
    do j=1,r
      if(j<=m)loadings(j,j)=cutoff+0.5_dp
      do i=j+1,m
        if(mod(i+j,3)==0)loadings(i,j)=sign(cutoff+0.1_dp*real(mod(i,5)+1,dp),real(mod(i+j,2)*2-1,dp))
      end do
      target=min(3,m-j+1);nactive=count(abs(loadings(:,j))>cutoff)
      i=j
      do while(nactive<target.and.i<=m)
        if(abs(loadings(i,j))<=cutoff)then;loadings(i,j)=cutoff+0.2_dp;nactive=nactive+1;end if
        i=i+1
      end do
    end do
  end subroutine make_sparse_loadings

  real(dp) function randu_local() result(x)
    use sv_rng, only : randu
    x=randu()
  end function randu_local

  subroutine default_fsv_parameters(m,r,idio_params,factor_params)
    integer,intent(in)::m,r
    real(dp),intent(out)::idio_params(m,3),factor_params(r,2)
    integer::i
    do i=1,m
      idio_params(i,1)=-2.0_dp+0.9_dp*real(i-1,dp)/real(max(1,m-1),dp)
      idio_params(i,2)=0.8_dp+0.18_dp*real(i-1,dp)/real(max(1,m-1),dp)
      idio_params(i,3)=0.6_dp-0.45_dp*real(i-1,dp)/real(max(1,m-1),dp)
    end do
    do i=1,r
      select case(i)
      case(1);factor_params(i,:)=[0.99_dp,0.1_dp]
      case(2);factor_params(i,:)=[0.95_dp,0.3_dp]
      case(3);factor_params(i,:)=[0.97_dp,0.1_dp]
      case default;factor_params(i,:)=[0.97_dp,0.2_dp]
      end select
    end do
  end subroutine default_fsv_parameters

  subroutine simulate_logvol(mu,phi,sigma,n,h0,h,hetero)
    real(dp),intent(in)::mu,phi,sigma
    integer,intent(in)::n
    real(dp),intent(out)::h0,h(n)
    logical,intent(in)::hetero
    integer::t
    if(.not.hetero)then;h0=mu;h=mu;return;end if
    h0=mu+sigma/sqrt(max(tiny_var,1.0_dp-phi*phi))*randn()
    h(1)=mu+phi*(h0-mu)+sigma*randn()
    do t=2,n;h(t)=mu+phi*(h(t-1)-mu)+sigma*randn();end do
  end subroutine simulate_logvol

  subroutine simulate_fsv(n,loadings,idio_params,factor_params,out,heteroskedastic,df)
    integer,intent(in)::n
    real(dp),intent(in)::loadings(:,:),idio_params(:,:),factor_params(:,:)
    type(fsv_sim_result),intent(out)::out
    logical,intent(in),optional::heteroskedastic(:)
    real(dp),intent(in),optional::df
    integer::m,r,i,t
    logical,allocatable::het(:)
    real(dp)::dof
    m=size(loadings,1);r=size(loadings,2);dof=huge(1.0_dp);if(present(df))dof=df
    allocate(het(m+r));het=.true.;if(present(heteroskedastic))het=heteroskedastic
    allocate(out%y(n,m),out%factors(n,r),out%loadings(m,r),out%h_idio(n,m),out%h_factor(n,r),out%h0_idio(m),out%h0_factor(r),out%idio_params(m,3),out%factor_params(r,2))
    out%loadings=loadings;out%idio_params=idio_params;out%factor_params=factor_params
    do i=1,m;call simulate_logvol(idio_params(i,1),idio_params(i,2),idio_params(i,3),n,out%h0_idio(i),out%h_idio(:,i),het(i));end do
    do i=1,r;call simulate_logvol(0.0_dp,factor_params(i,1),factor_params(i,2),n,out%h0_factor(i),out%h_factor(:,i),het(m+i));end do
    do t=1,n
      do i=1,r
        if(dof<huge(1.0_dp)/10.0_dp)then;out%factors(t,i)=exp(.5_dp*out%h_factor(t,i))*rand_student_t(dof);else;out%factors(t,i)=exp(.5_dp*out%h_factor(t,i))*randn();end if
      end do
      out%y(t,:)=matmul(loadings,out%factors(t,:))
      do i=1,m;out%y(t,i)=out%y(t,i)+exp(.5_dp*out%h_idio(t,i))*randn();end do
    end do
  end subroutine simulate_fsv

  subroutine fsv_covariance_path(loadings,h_idio,h_factor,cov)
    real(dp),intent(in)::loadings(:,:),h_idio(:,:),h_factor(:,:)
    real(dp),allocatable,intent(out)::cov(:,:,:)
    real(dp),allocatable::tmp(:,:)
    integer::t,i,m,r,n
    n=size(h_idio,1);m=size(h_idio,2);r=size(h_factor,2);allocate(cov(m,m,n),tmp(m,r))
    do t=1,n
      do i=1,r;tmp(:,i)=loadings(:,i)*exp(.5_dp*h_factor(t,i));end do
      cov(:,:,t)=matmul(tmp,transpose(tmp));do i=1,m;cov(i,i,t)=cov(i,i,t)+exp(h_idio(t,i));end do
    end do
  end subroutine fsv_covariance_path

  subroutine fsv_correlation_path(loadings,h_idio,h_factor,cor)
    real(dp),intent(in)::loadings(:,:),h_idio(:,:),h_factor(:,:)
    real(dp),allocatable,intent(out)::cor(:,:,:)
    real(dp),allocatable::cov(:,:,:)
    integer::t
    call fsv_covariance_path(loadings,h_idio,h_factor,cov);allocate(cor,source=cov)
    do t=1,size(cov,3);call correlation_from_covariance(cov(:,:,t),cor(:,:,t));end do
  end subroutine fsv_correlation_path

  subroutine covelement(loadings,h_idio,h_factor,i,j,x)
    real(dp),intent(in)::loadings(:,:),h_idio(:,:),h_factor(:,:)
    integer,intent(in)::i,j
    real(dp),intent(out)::x(:)
    integer::t,k
    do t=1,size(h_idio,1);x(t)=0.0_dp;do k=1,size(loadings,2);x(t)=x(t)+loadings(i,k)*loadings(j,k)*exp(h_factor(t,k));end do;if(i==j)x(t)=x(t)+exp(h_idio(t,i));end do
  end subroutine covelement
  subroutine corelement(loadings,h_idio,h_factor,i,j,x)
    real(dp),intent(in)::loadings(:,:),h_idio(:,:),h_factor(:,:)
    integer,intent(in)::i,j
    real(dp),intent(out)::x(:)
    real(dp),allocatable::cii(:),cjj(:)
    allocate(cii(size(x)),cjj(size(x)));call covelement(loadings,h_idio,h_factor,i,j,x);call covelement(loadings,h_idio,h_factor,i,i,cii);call covelement(loadings,h_idio,h_factor,j,j,cjj);x=x/sqrt(cii*cjj)
  end subroutine corelement

  subroutine expweightcov(dat,alpha,hist,cov)
    real(dp),intent(in)::dat(:,:),alpha
    integer,intent(in)::hist
    real(dp),intent(out)::cov(:,:)
    integer::n,k,t
    n=size(dat,1);k=min(hist,n);cov=0.0_dp
    cov=outer(dat(n-k+1,:),dat(n-k+1,:))
    do t=n-k+2,n;cov=(1.0_dp-alpha)*cov+alpha*outer(dat(t,:),dat(t,:));end do
  end subroutine expweightcov
  pure function outer(a,b) result(c)
    real(dp),intent(in)::a(:),b(:);real(dp)::c(size(a),size(b));integer::i
    do i=1,size(a);c(i,:)=a(i)*b;end do
  end function outer

  subroutine static_factor_initialize(dat,r,loadings,factors,uniq,info)
    real(dp),intent(in)::dat(:,:)
    integer,intent(in)::r
    real(dp),allocatable,intent(out)::loadings(:,:),factors(:,:),uniq(:)
    integer,intent(out)::info
    real(dp),allocatable::cov(:,:),eval(:),evec(:,:),xc(:,:),means(:),prec(:,:),rhs(:),coef(:)
    integer::m,n,i,j,inf
    n=size(dat,1);m=size(dat,2);allocate(cov(m,m));call covariance_matrix(dat,cov)
    if(r==0)then
      allocate(loadings(m,0),factors(n,0),uniq(m));uniq=diagv(cov);info=0;return
    end if
    call symmetric_eigen(cov,eval,evec,info);if(info/=0)return
    allocate(loadings(m,r),uniq(m),factors(n,r),means(m),xc(n,m));do j=1,m;means(j)=sum(dat(:,j))/real(n,dp);xc(:,j)=dat(:,j)-means(j);end do
    do j=1,r;loadings(:,j)=evec(:,m-j+1)*sqrt(max(eval(m-j+1),tiny_var));end do
    uniq=diagv(cov-matmul(loadings,transpose(loadings)));uniq=max(uniq,0.05_dp*diagv(cov))
    allocate(prec(r,r),rhs(r),coef(r));prec=matmul(transpose(loadings),loadings);do i=1,r;prec(i,i)=prec(i,i)+1.0e-6_dp;end do
    do i=1,n;rhs=matmul(transpose(loadings),xc(i,:));call solve_spd(prec,rhs,coef,inf);factors(i,:)=coef;end do
  end subroutine static_factor_initialize
  pure function diagv(a) result(d)
    real(dp),intent(in)::a(:,:);real(dp)::d(min(size(a,1),size(a,2)));integer::i
    do i=1,size(d);d(i)=a(i,i);end do
  end function diagv

  subroutine preorder(dat,r,ordering,info)
    real(dp),intent(in)::dat(:,:);integer,intent(in)::r;integer,intent(out)::ordering(:);integer,intent(out)::info
    real(dp),allocatable::l(:,:),f(:,:),u(:);logical,allocatable::used(:);integer::j,i,best,m
    m=size(dat,2);call static_factor_initialize(dat,r,l,f,u,info);if(info/=0)return;allocate(used(m));used=.false.
    do j=1,r
      best=0
      do i=1,m
        if(.not.used(i))then
          if(best==0)then
            best=i
          else if(abs(l(i,j))>abs(l(best,j)))then
            best=i
          end if
        end if
      end do
      ordering(j)=best;used(best)=.true.
    end do
    best=r;do i=1,m;if(.not.used(i))then;best=best+1;ordering(best)=i;end if;end do
  end subroutine preorder

  subroutine findrestrict(dat,r,restrict,info)
    real(dp),intent(in)::dat(:,:);integer,intent(in)::r;logical,intent(out)::restrict(:,:);integer,intent(out)::info
    integer,allocatable::ord(:);integer::j
    allocate(ord(size(dat,2)));call preorder(dat,r,ord,info);restrict=.false.;if(info/=0)return
    do j=1,r-1;restrict(ord(j),j+1:r)=.true.;end do
  end subroutine findrestrict

  subroutine fit_fsv(y,r,opt,out,restriction)
    real(dp),intent(in)::y(:,:)
    integer,intent(in)::r
    type(fsv_options),intent(in)::opt
    type(fsv_draws),intent(out)::out
    logical,intent(in),optional::restriction(:,:)
    integer::n,m,i,j,t,iter,total,keep,info,acc(7)
    real(dp),allocatable::b(:,:),f(:,:),uniq(:),hid(:,:),hf(:,:),h0i(:),h0f(:),tau(:),prec(:,:),rhs(:),mean(:),cov(:,:),draw(:),res(:),para(:,:),local_scale(:,:),global_shrink(:)
    logical,allocatable::active(:,:)
    type(sv_params)::sp
    type(sv_prior)::prior
    type(sv_mcmc_options)::sopt
    n=size(y,1);m=size(y,2);call static_factor_initialize(y,r,b,f,uniq,info);if(info/=0)error stop 'fit_fsv: initialization failed'
    allocate(active(m,r));active=.true.
    if(opt%lower_triangular)then;do j=1,r;do i=1,min(j-1,m);active(i,j)=.false.;end do;end do;end if
    if(present(restriction))then
      if(any(shape(restriction)/=[m,r]))error stop 'fit_fsv: restriction shape mismatch'
      active=restriction
    end if
    where(.not.active)b=0.0_dp
    if(opt%lower_triangular.and.r>0)then;do j=1,min(r,m);if(active(j,j).and.b(j,j)<0.0_dp)b(:,j)=-b(:,j);end do;end if
    allocate(hid(n,m),hf(n,r),h0i(m),h0f(r),para(3,m+r),tau(n),local_scale(m,r),global_shrink(r));tau=1.0_dp;local_scale=opt%loading_prior_sd**2;global_shrink=opt%global_shrinkage
    do i=1,m;hid(:,i)=log(max((y(:,i)-matmul(f,b(i,:)))**2,1.0e-4_dp));h0i(i)=sum(hid(:,i))/real(n,dp);para(:,i)=[h0i(i),.95_dp,.2_dp];end do
    do j=1,r;hf(:,j)=log(max(f(:,j)**2,1.0e-4_dp));h0f(j)=sum(hf(:,j))/real(n,dp);para(:,m+j)=[0.0_dp,.95_dp,.2_dp];end do
    total=opt%burnin+opt%draws*opt%thin;out%nobs=n;out%nseries=m;out%nfactors=r;out%ndraws=opt%draws
    allocate(out%loadings(m,r,opt%draws),out%para(3,m+r,opt%draws));if(opt%normal_gamma)allocate(out%local_scale(m,r,opt%draws),out%global_shrinkage(r,opt%draws));if(opt%store_factors)allocate(out%factors(n,r,opt%draws));if(opt%store_latent)allocate(out%latent(n,m+r,opt%draws))
    sopt%draws=1;sopt%burnin=0;sopt%thin=1;sopt%latent_sweeps=opt%sv_sweeps;sopt%use_mixture=.true.;sopt%sample_rho=.false.;sopt%sample_nu=.false.;sopt%store_latent=.false.
    keep=0;acc=0
    do iter=1,total
      if(r>0)then
      ! Draw factors conditionally, time by time.
      allocate(prec(r,r),rhs(r),mean(r),cov(r,r),draw(r))
      do t=1,n
        prec=0.0_dp;rhs=0.0_dp
        do j=1,r;prec(j,j)=exp(-hf(t,j));end do
        do i=1,m
          prec=prec+outer(b(i,:),b(i,:))*exp(-hid(t,i));rhs=rhs+b(i,:)*y(t,i)*exp(-hid(t,i))
        end do
        call inverse_spd(prec,cov,info);if(info/=0)error stop 'fit_fsv: factor precision failed';mean=matmul(cov,rhs);call mvn_draw(mean,cov,draw,info);f(t,:)=draw
      end do
      deallocate(prec,rhs,mean,cov,draw)
      ! Columnwise Normal-Gamma shrinkage update used by factorstochvol.
      if(opt%normal_gamma)then
        do j=1,r
          global_shrink(j)=rand_gamma(opt%ng_c+opt%ng_a*real(m,dp),opt%ng_d+.5_dp*opt%ng_a*sum(local_scale(:,j)))
          do i=1,m
            if(active(i,j))local_scale(i,j)=rand_gig_slice(opt%ng_a-.5_dp,b(i,j)*b(i,j),opt%ng_a*global_shrink(j),local_scale(i,j))
          end do
        end do
      end if
      ! Draw each loading row from weighted Gaussian regression.
      do i=1,m
        allocate(prec(r,r),rhs(r),mean(r),cov(r,r),draw(r));prec=0.0_dp;rhs=0.0_dp
        do j=1,r
          if(.not.active(i,j))then
            prec(j,j)=1.0e12_dp
          else if(opt%normal_gamma)then
            prec(j,j)=1.0_dp/max(tiny_var,local_scale(i,j))
          else
            prec(j,j)=opt%global_shrinkage/max(tiny_var,opt%loading_prior_sd**2)
          end if
        end do
        do t=1,n;prec=prec+outer(f(t,:),f(t,:))*exp(-hid(t,i));rhs=rhs+f(t,:)*y(t,i)*exp(-hid(t,i));end do
        call inverse_spd(prec,cov,info);if(info/=0)error stop 'fit_fsv: loading precision failed';mean=matmul(cov,rhs);call mvn_draw(mean,cov,draw,info);b(i,:)=draw
        do j=1,r;if(.not.active(i,j))b(i,j)=0.0_dp;end do
        deallocate(prec,rhs,mean,cov,draw)
      end do
      if(opt%lower_triangular)then;do j=1,min(r,m);if(active(j,j).and.b(j,j)<0.0_dp)then;b(:,j)=-b(:,j);f(:,j)=-f(:,j);end if;end do;end if
      end if
      ! Reuse univariate SV updates for idiosyncratic residuals and factors.
      do i=1,m
        allocate(res(n));res=y(:,i)-matmul(f,b(i,:));sp%mu=para(1,i);sp%phi=para(2,i);sp%sigma=para(3,i);sp%rho=0.0_dp;sp%nu=huge(1.0_dp)
        call sv_single_sweep(res,sp,hid(:,i),h0i(i),tau,prior,sopt,accept=acc);para(:,i)=[sp%mu,sp%phi,sp%sigma];deallocate(res)
      end do
      sopt%sample_mu=.false.
      do j=1,r
        sp%mu=0.0_dp;sp%phi=para(2,m+j);sp%sigma=para(3,m+j);sp%rho=0.0_dp;sp%nu=huge(1.0_dp)
        call sv_single_sweep(f(:,j),sp,hf(:,j),h0f(j),tau,prior,sopt,accept=acc);para(:,m+j)=[0.0_dp,sp%phi,sp%sigma]
      end do
      sopt%sample_mu=.true.
      if(iter>opt%burnin.and.mod(iter-opt%burnin,opt%thin)==0)then;keep=keep+1;out%loadings(:,:,keep)=b;out%para(:,:,keep)=para;if(opt%store_factors)out%factors(:,:,keep)=f;if(opt%store_latent)then;out%latent(:,1:m,keep)=hid;out%latent(:,m+1:m+r,keep)=hf;end if;if(opt%normal_gamma)then;out%local_scale(:,:,keep)=local_scale;out%global_shrinkage(:,keep)=global_shrink;end if;end if
    end do
  end subroutine fit_fsv

  subroutine predict_fsv(draws,steps,each,pred)
    type(fsv_draws),intent(in)::draws
    integer,intent(in)::steps,each
    type(fsv_prediction),intent(out)::pred
    integer::m,r,d,e,s,k,i,info
    real(dp),allocatable::hcur(:),hnew(:),cov(:,:),cor(:,:),inv(:,:)
    m=draws%nseries;r=draws%nfactors;allocate(pred%h(steps,m+r,draws%ndraws*each),pred%cov(m,m,steps,draws%ndraws*each),pred%cor(m,m,steps,draws%ndraws*each),pred%precision(m,m,steps,draws%ndraws*each),pred%logdet_precision(steps,draws%ndraws*each));allocate(hcur(m+r),hnew(m+r),cov(m,m),cor(m,m),inv(m,m))
    k=0;do d=1,draws%ndraws;do e=1,each;k=k+1;if(allocated(draws%latent))then;hcur=draws%latent(draws%nobs,:,d);else;hcur=draws%para(1,:,d);end if
      do s=1,steps
        do i=1,m+r;hnew(i)=draws%para(1,i,d)+draws%para(2,i,d)*(hcur(i)-draws%para(1,i,d))+draws%para(3,i,d)*randn();end do;hcur=hnew;pred%h(s,:,k)=hcur
        call one_cov(draws%loadings(:,:,d),hcur(1:m),hcur(m+1:m+r),cov);call correlation_from_covariance(cov,cor)
        call woodbury_precision(draws%loadings(:,:,d),hcur(1:m),hcur(m+1:m+r),inv,pred%logdet_precision(s,k),info)
        pred%cov(:,:,s,k)=cov;pred%cor(:,:,s,k)=cor;if(info==0)then;pred%precision(:,:,s,k)=inv;else;pred%precision(:,:,s,k)=0.0_dp;end if
      end do
    end do;end do
  end subroutine predict_fsv
  subroutine one_cov(b,hi,hf,cov)
    real(dp),intent(in)::b(:,:),hi(:),hf(:);real(dp),intent(out)::cov(:,:);real(dp),allocatable::tmp(:,:);integer::j,i
    allocate(tmp(size(b,1),size(b,2)));do j=1,size(b,2);tmp(:,j)=b(:,j)*exp(.5_dp*hf(j));end do;cov=matmul(tmp,transpose(tmp));do i=1,size(b,1);cov(i,i)=cov(i,i)+exp(hi(i));end do
  end subroutine one_cov

  subroutine predloglik_fsv(pred,y,loglik)
    type(fsv_prediction),intent(in)::pred
    real(dp),intent(in)::y(:,:)
    real(dp),intent(out)::loglik(:,:)
    real(dp),allocatable::zero(:)
    integer::s,k
    allocate(zero(size(y,2)));zero=0.0_dp
    do s=1,size(y,1);do k=1,size(pred%cov,4);loglik(s,k)=mvn_logpdf(y(s,:),zero,pred%cov(:,:,s,k));end do;end do
  end subroutine predloglik_fsv

  subroutine predloglik_fsv_woodbury(pred,y,loglik)
    type(fsv_prediction),intent(in)::pred
    real(dp),intent(in)::y(:,:)
    real(dp),intent(out)::loglik(:,:)
    integer::s,k,m
    m=size(y,2)
    do s=1,size(y,1);do k=1,size(pred%precision,4)
      loglik(s,k)=0.5_dp*(pred%logdet_precision(s,k)-real(m,dp)*log2pi-dot_product(y(s,:),matmul(pred%precision(:,:,s,k),y(s,:))))
    end do;end do
  end subroutine predloglik_fsv_woodbury

  subroutine aggregate_loglik_draws(logdraws,scores)
    real(dp),intent(in)::logdraws(:,:)
    real(dp),intent(out)::scores(:)
    real(dp)::mx
    integer::s
    do s=1,size(logdraws,1);mx=maxval(logdraws(s,:));scores(s)=mx+log(sum(exp(logdraws(s,:)-mx))/real(size(logdraws,2),dp));end do
  end subroutine aggregate_loglik_draws

  subroutine dmvnorm_columns(x,means,vars,loga,out)
    real(dp),intent(in)::x(:,:),means(:,:),vars(:,:,:);logical,intent(in)::loga;real(dp),intent(out)::out(:);integer::i
    do i=1,size(x,2);out(i)=mvn_logpdf(x(:,i),means(:,i),vars(:,:,i));if(.not.loga)out(i)=exp(out(i));end do
  end subroutine dmvnorm_columns

  subroutine sign_identify(loadings,factors)
    real(dp),intent(inout)::loadings(:,:,:)
    real(dp),intent(inout),optional::factors(:,:,:)
    integer::d,j,i
    do d=1,size(loadings,3);do j=1,size(loadings,2);i=maxloc(abs(loadings(:,j,d)),dim=1);if(loadings(i,j,d)<0.0_dp)then;loadings(:,j,d)=-loadings(:,j,d);if(present(factors))factors(:,j,d)=-factors(:,j,d);end if;end do;end do
  end subroutine sign_identify

  subroutine order_identify(loadings,factors)
    real(dp),intent(inout)::loadings(:,:,:)
    real(dp),intent(inout),optional::factors(:,:,:)
    real(dp),allocatable::score(:),tmpb(:,:),tmpf(:,:)
    integer,allocatable::ord(:);integer::d,j,r
    r=size(loadings,2);allocate(score(r),ord(r),tmpb(size(loadings,1),r));if(present(factors))allocate(tmpf(size(factors,1),r))
    do d=1,size(loadings,3);do j=1,r;score(j)=sum(loadings(:,j,d)**2);ord(j)=j;end do;call sort_desc(score,ord);tmpb=loadings(:,:,d);do j=1,r;loadings(:,j,d)=tmpb(:,ord(j));end do;if(present(factors))then;tmpf=factors(:,:,d);do j=1,r;factors(:,j,d)=tmpf(:,ord(j));end do;end if;end do
  end subroutine order_identify
  subroutine sort_desc(x,idx)
    real(dp),intent(inout)::x(:);integer,intent(inout)::idx(:);integer::i,j,it;real(dp)::v
    do i=2,size(x)
      v=x(i);it=idx(i);j=i-1
      do while(j>=1)
        if(x(j)>=v)exit
        x(j+1)=x(j);idx(j+1)=idx(j);j=j-1
      end do
      x(j+1)=v;idx(j+1)=it
    end do
  end subroutine sort_desc

  subroutine running_covariance(draws,time_index,covmean)
    type(fsv_draws),intent(in)::draws;integer,intent(in)::time_index;real(dp),intent(out)::covmean(:,:)
    real(dp),allocatable::cov(:,:);integer::d,m,r
    m=draws%nseries;r=draws%nfactors;allocate(cov(m,m));covmean=0.0_dp
    do d=1,draws%ndraws;call one_cov(draws%loadings(:,:,d),draws%latent(time_index,1:m,d),draws%latent(time_index,m+1:m+r,d),cov);covmean=covmean+cov;end do;covmean=covmean/real(draws%ndraws,dp)
  end subroutine running_covariance
  subroutine running_correlation(draws,time_index,cormean)
    type(fsv_draws),intent(in)::draws;integer,intent(in)::time_index;real(dp),intent(out)::cormean(:,:)
    real(dp),allocatable::cov(:,:),cor(:,:);integer::d,m,r
    m=draws%nseries;r=draws%nfactors;allocate(cov(m,m),cor(m,m));cormean=0.0_dp
    do d=1,draws%ndraws;call one_cov(draws%loadings(:,:,d),draws%latent(time_index,1:m,d),draws%latent(time_index,m+1:m+r,d),cov);call correlation_from_covariance(cov,cor);cormean=cormean+cor;end do;cormean=cormean/real(draws%ndraws,dp)
  end subroutine running_correlation

  subroutine woodbury_precision(loadings,h_idio,h_factor,precision,logdet_precision,info)
    real(dp),intent(in)::loadings(:,:),h_idio(:),h_factor(:)
    real(dp),intent(out)::precision(:,:),logdet_precision
    integer,intent(out)::info
    real(dp),allocatable::middle(:,:),middle_inv(:,:),dinvb(:,:)
    integer::m,r,i,j
    m=size(loadings,1);r=size(loadings,2)
    if(r==0)then
      precision=0.0_dp;do i=1,m;precision(i,i)=exp(-h_idio(i));end do;logdet_precision=-sum(h_idio);info=0;return
    end if
    allocate(middle(r,r),dinvb(m,r))
    do j=1,r;do i=1,m;dinvb(i,j)=loadings(i,j)*exp(-h_idio(i));end do;end do
    middle=matmul(transpose(loadings),dinvb)
    do j=1,r;middle(j,j)=middle(j,j)+exp(-h_factor(j));end do
    call inverse_spd(middle,middle_inv,info);if(info/=0)return
    precision=-matmul(dinvb,matmul(middle_inv,transpose(dinvb)))
    do i=1,m;precision(i,i)=precision(i,i)+exp(-h_idio(i));end do
    logdet_precision=-sum(h_idio)-sum(h_factor)-logdet_small(middle,info)
  end subroutine woodbury_precision

  real(dp) function logdet_small(a,info) result(v)
    real(dp),intent(in)::a(:,:);integer,intent(out)::info
    real(dp),allocatable::ev(:),vec(:,:);integer::i
    call symmetric_eigen(a,ev,vec,info);v=0.0_dp;if(info/=0)return
    if(any(ev<=0.0_dp))then;info=1;return;end if
    do i=1,size(ev);v=v+log(ev(i));end do
  end function logdet_small

  subroutine predict_conditional_fsv(draws,steps,each,means,vols)
    type(fsv_draws),intent(in)::draws
    integer,intent(in)::steps,each
    real(dp),allocatable,intent(out)::means(:,:,:),vols(:,:,:)
    real(dp),allocatable::hcur(:),factor(:)
    integer::m,r,d,e,s,k,i
    m=draws%nseries;r=draws%nfactors;allocate(means(m,steps,draws%ndraws*each),vols(m,steps,draws%ndraws*each),hcur(m+r),factor(r))
    k=0
    do d=1,draws%ndraws;do e=1,each
      k=k+1;hcur=draws%latent(draws%nobs,:,d)
      do s=1,steps
        do i=1,m+r;hcur(i)=draws%para(1,i,d)+draws%para(2,i,d)*(hcur(i)-draws%para(1,i,d))+draws%para(3,i,d)*randn();end do
        do i=1,r;factor(i)=exp(.5_dp*hcur(m+i))*randn();end do
        means(:,s,k)=matmul(draws%loadings(:,:,d),factor);vols(:,s,k)=exp(.5_dp*hcur(1:m))
      end do
    end do;end do
  end subroutine predict_conditional_fsv

  subroutine eigen_loading_diagnostics(draws,eigenvalues,info)
    type(fsv_draws),intent(in)::draws
    real(dp),allocatable,intent(out)::eigenvalues(:,:)
    integer,intent(out)::info
    real(dp),allocatable::a(:,:),ev(:),vec(:,:)
    integer::d,r
    r=draws%nfactors;allocate(eigenvalues(r,draws%ndraws),a(r,r));info=0
    do d=1,draws%ndraws
      a=matmul(transpose(draws%loadings(:,:,d)),draws%loadings(:,:,d));call symmetric_eigen(a,ev,vec,info);if(info/=0)return;eigenvalues(:,d)=ev
    end do
  end subroutine eigen_loading_diagnostics
end module fsv_core
