! Modern Fortran translation of the computational core of GPareto 1.1.9.
! GPareto is GPL-3.0-only; see LICENSE and UPSTREAM.md.
module gpareto_sur
  use gpareto_kinds, only : dp, i8
  use gpareto_math, only : rng_state
  use gpareto_models, only : gp_model_set, make_trend
  use gpareto_probability, only : probability_nondomination
  use gpareto_criteria, only : check_predict
  use dk_model, only : km_prediction, km_predict
  implicit none
  private
  public :: crit_sur
contains
  subroutine crit_sur(x,models,front,integration,value,weights,nsamp,seed,kind,threshold)
    real(dp),intent(in)::x(:,:),front(:,:),integration(:,:)
    type(gp_model_set),intent(in)::models
    real(dp),allocatable,intent(out)::value(:)
    real(dp),intent(in),optional::weights(:),threshold
    integer,intent(in),optional::nsamp
    integer(i8),intent(in),optional::seed
    character(len=*),intent(in),optional::kind
    real(dp),allocatable::base_mean(:,:),base_sd(:,:),pn(:),w(:),joint(:,:),f(:,:),c(:)
    real(dp),allocatable::mnew(:,:),snew(:,:)
    type(km_prediction)::pr
    type(rng_state)::rng
    character(len=2)::kt
    real(dp)::u0,ue,varx,z
    integer::nq,m,nc,ns,i,j,k
    nq=size(integration,1)
    m=models%nobj()
    nc=size(x,1)
    ns=20
    if(present(nsamp))ns=nsamp
    kt='UK'
    if(present(kind))kt=kind
    allocate(w(nq))
    if(present(weights))then
    if(size(weights)/=nq)error stop 'crit_sur: weights'
    w=weights
    else
    w=1.0_dp/real(nq,dp)
    end if
    allocate(base_mean(nq,m),base_sd(nq,m))
    do j=1,m
      call make_trend(models%model(j)%trend_kind,integration,f)
      call km_predict(models%model(j)%km,integration,f,kt,pr,se_compute=.true.,cov_compute=.false.)
      base_mean(:,j)=pr%mean
      base_sd(:,j)=pr%sd
    end do
    call probability_nondomination(front,base_mean,base_sd,pn)
    u0=sum(w*pn*(1.0_dp-pn))
    allocate(value(nc))
    call rng%seed(42_i8)
    if(present(seed))call rng%seed(seed)
    do i=1,nc
      if(check_predict(x(i,:),models,threshold))then
      value(i)=-1.0_dp
      cycle
      end if
      ue=0.0_dp
      do k=1,ns
        allocate(mnew(nq,m),snew(nq,m))
        do j=1,m
          allocate(joint(nq+1,size(x,2)))
          joint(1:nq,:)=integration
          joint(nq+1,:)=x(i,:)
          call make_trend(models%model(j)%trend_kind,joint,f)
          call km_predict(models%model(j)%km,joint,f,kt,pr,se_compute=.true.,cov_compute=.true.)
          varx=max(pr%cov(nq+1,nq+1),tiny(1.0_dp))
          allocate(c(nq))
          c=pr%cov(1:nq,nq+1)
          z=rng%normal()
          mnew(:,j)=base_mean(:,j)+c/sqrt(varx)*z
          snew(:,j)=sqrt(max(base_sd(:,j)**2-c*c/varx,0.0_dp))
          deallocate(joint,c)
        end do
        call probability_nondomination(front,mnew,snew,pn)
        ue=ue+sum(w*pn*(1.0_dp-pn))
        deallocate(mnew,snew)
      end do
      value=u0-ue/real(ns,dp)
    end do
  end subroutine crit_sur
end module gpareto_sur
