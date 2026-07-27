! SPDX-License-Identifier: GPL-2.0-or-later
! Numerical translation derived from tsDyn, copyright its original authors.
! This file may be redistributed and/or modified under GPL version 2 or later.
module tsdyn_girf
  use tsdyn_kinds, only: dp, n_deterministic
  use tsdyn_utils, only: seed_random, random_normal_vector
  use tsdyn_setar, only: setar_model, simulate_setar
  use tsdyn_tvar, only: tvar_model, simulate_tvar
  use tsdyn_tvecm, only: tvecm_model, simulate_tvecm
  implicit none
  private
  public :: regime_irf_setar, regime_irf_tvar
  public :: girf_setar, girf_tvar, girf_tvecm
contains
  subroutine regime_irf_setar(model,h,regime,irf,info,cumulative)
    type(setar_model),intent(in)::model
    integer,intent(in)::h,regime
    real(dp),allocatable,intent(out)::irf(:)
    integer,intent(out)::info
    logical,intent(in),optional::cumulative
    real(dp),allocatable::work(:)
    integer::j,t,nd
    logical::cum
    if(regime<1.or.regime>model%nregime.or.h<0)then;info=-1;allocate(irf(0));return;end if
    nd=n_deterministic(model%include);allocate(irf(0:h),work(model%pmax+h+1));work=0.0_dp;work(model%pmax+1)=1.0_dp;irf(0)=1.0_dp
    do t=1,h
      do j=1,model%orders(regime);work(model%pmax+1+t)=work(model%pmax+1+t)+model%coefficients(nd+j,regime)*work(model%pmax+1+t-j);end do
      irf(t)=work(model%pmax+1+t)
    end do
    cum=.false.;if(present(cumulative))cum=cumulative;if(cum)then;do t=1,h;irf(t)=irf(t)+irf(t-1);end do;end if;info=0
  end subroutine regime_irf_setar

  subroutine regime_irf_tvar(model,h,regime,irf,info,cumulative,orthogonal)
    use tsdyn_linalg, only: cholesky_lower
    type(tvar_model),intent(in)::model
    integer,intent(in)::h,regime
    real(dp),allocatable,intent(out)::irf(:,:,:)
    integer,intent(out)::info
    logical,intent(in),optional::cumulative,orthogonal
    real(dp),allocatable::a(:,:,:),psi(:,:,:),impact(:,:),l(:,:)
    integer::k,p,nd,t,j,istat
    logical::cum,ortho
    if(regime<1.or.regime>model%nregime.or.h<0)then;info=-1;allocate(irf(0,0,0));return;end if
    k=model%nvar;p=model%order;nd=n_deterministic(model%include);allocate(a(k,k,p))
    do j=1,p;a(:,:,j)=transpose(model%coefficients(nd+(j-1)*k+1:nd+j*k,regime,:));end do
    allocate(psi(0:h,k,k),impact(k,k));psi=0.0_dp;impact=identity(k);ortho=.false.;if(present(orthogonal))ortho=orthogonal
    if(ortho)then;call cholesky_lower(model%sigma,l,istat);if(istat/=0)then;info=istat;return;end if;impact=l;end if
    psi(0,:,:)=identity(k)
    do t=1,h
      do j=1,min(p,t);psi(t,:,:)=psi(t,:,:)+matmul(a(:,:,j),psi(t-j,:,:));end do
    end do
    allocate(irf(0:h,k,k));do t=0,h;irf(t,:,:)=matmul(psi(t,:,:),impact);end do
    cum=.false.;if(present(cumulative))cum=cumulative;if(cum)then;do t=1,h;irf(t,:,:)=irf(t,:,:)+irf(t-1,:,:);end do;end if;info=0
  end subroutine regime_irf_tvar

  pure function identity(n) result(a)
    integer,intent(in)::n
    real(dp)::a(n,n)
    integer::i
    a=0.0_dp;do i=1,n;a(i,i)=1.0_dp;end do
  end function identity

  subroutine girf_setar(model,history,shock,h,nrep,response,info,seed)
    type(setar_model),intent(in)::model
    real(dp),intent(in)::history(:),shock
    integer,intent(in)::h,nrep
    real(dp),allocatable,intent(out)::response(:)
    integer,intent(out)::info
    integer,intent(in),optional::seed
    real(dp),allocatable::innov(:),innov_shock(:),base(:),hit(:)
    integer::r,istat,i
    if(size(history)<model%pmax.or.h<1.or.nrep<1)then;info=-1;allocate(response(0));return;end if
    if(present(seed))call seed_random(seed);allocate(response(h),innov(h),innov_shock(h));response=0.0_dp
    do r=1,nrep
      do i=1,h;innov(i)=sqrt(max(model%sigma2,tiny(1.0_dp)))*random_scalar_normal();end do
      innov_shock=innov;innov_shock(1)=innov_shock(1)+shock
      call simulate_setar(model,h,base,innov=innov,start=history,info=istat);if(istat/=0)then;info=istat;return;end if
      call simulate_setar(model,h,hit,innov=innov_shock,start=history,info=istat);if(istat/=0)then;info=istat;return;end if
      response=response+(hit-base)
    end do
    response=response/real(nrep,dp);info=0
  end subroutine girf_setar

  subroutine girf_tvar(model,history,shock,h,nrep,response,info,seed)
    type(tvar_model),intent(in)::model
    real(dp),intent(in)::history(:,:),shock(:)
    integer,intent(in)::h,nrep
    real(dp),allocatable,intent(out)::response(:,:)
    integer,intent(out)::info
    integer,intent(in),optional::seed
    real(dp),allocatable::innov(:,:),innov_shock(:,:),base(:,:),hit(:,:),e(:)
    integer::r,istat,i
    if(size(history,1)<model%order.or.size(history,2)/=model%nvar.or.size(shock)/=model%nvar.or.h<1.or.nrep<1)then;info=-1;allocate(response(0,0));return;end if
    if(present(seed))call seed_random(seed);allocate(response(h,model%nvar),innov(h,model%nvar),innov_shock(h,model%nvar),e(model%nvar));response=0.0_dp
    do r=1,nrep
      do i=1,h;call random_normal_vector(e);innov(i,:)=matmul_chol(model%sigma,e,istat);if(istat/=0)then;info=istat;return;end if;end do
      innov_shock=innov;innov_shock(1,:)=innov_shock(1,:)+shock
      call simulate_tvar(model,h,base,innov=innov,start=history,info=istat);if(istat/=0)then;info=istat;return;end if
      call simulate_tvar(model,h,hit,innov=innov_shock,start=history,info=istat);if(istat/=0)then;info=istat;return;end if
      response=response+(hit-base)
    end do
    response=response/real(nrep,dp);info=0
  end subroutine girf_tvar

  subroutine girf_tvecm(model,history,shock,h,nrep,response,info,seed)
    type(tvecm_model),intent(in)::model
    real(dp),intent(in)::history(:,:),shock(:)
    integer,intent(in)::h,nrep
    real(dp),allocatable,intent(out)::response(:,:)
    integer,intent(out)::info
    integer,intent(in),optional::seed
    real(dp),allocatable::innov(:,:),innov_shock(:,:),base(:,:),hit(:,:),e(:)
    integer::r,istat,i
    if(size(history,1)<model%lag_diff+1.or.size(history,2)/=model%nvar.or.size(shock)/=model%nvar.or.h<1.or.nrep<1)then;info=-1;allocate(response(0,0));return;end if
    if(present(seed))call seed_random(seed);allocate(response(h,model%nvar),innov(h,model%nvar),innov_shock(h,model%nvar),e(model%nvar));response=0.0_dp
    do r=1,nrep
      do i=1,h;call random_normal_vector(e);innov(i,:)=matmul_chol(model%sigma,e,istat);if(istat/=0)then;info=istat;return;end if;end do
      innov_shock=innov;innov_shock(1,:)=innov_shock(1,:)+shock
      call simulate_tvecm(model,h,base,innov=innov,start=history,info=istat);if(istat/=0)then;info=istat;return;end if
      call simulate_tvecm(model,h,hit,innov=innov_shock,start=history,info=istat);if(istat/=0)then;info=istat;return;end if
      response=response+(hit-base)
    end do
    response=response/real(nrep,dp);info=0
  end subroutine girf_tvecm

  real(dp) function random_scalar_normal() result(z)
    real(dp)::v(1)
    call random_normal_vector(v);z=v(1)
  end function random_scalar_normal

  function matmul_chol(sigma,e,info) result(v)
    use tsdyn_linalg, only: cholesky_lower
    real(dp),intent(in)::sigma(:,:),e(:)
    integer,intent(out)::info
    real(dp)::v(size(e))
    real(dp),allocatable::l(:,:)
    call cholesky_lower(sigma,l,info);if(info==0)then;v=matmul(l,e);else;v=0.0_dp;end if
  end function matmul_chol
end module tsdyn_girf
