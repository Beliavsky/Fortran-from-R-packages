! SPDX-License-Identifier: GPL-2.0-or-later
module flexsurv_shared_multistate
  use flexsurv_kinds, only : dp
  use flexsurv_fit, only : flexsurv_spec, flexsurv_result, predict_cumhaz, parameter_count
  use flexsurv_multistate, only : flexsurv_transition, msm_parametric, pmatrix_flexsurv
  use flexsurv_math, only : near_positive_definite, rng_normal
  implicit none
  private

  type, public :: flexsurv_shared_msm
    integer :: nstate = 0
    type(flexsurv_spec) :: spec
    real(dp), allocatable :: theta(:)
    real(dp), allocatable :: covariance(:,:)
    integer, allocatable :: from(:), to(:), row(:)
  end type flexsurv_shared_msm

  public :: make_shared_msm, shared_transitions, draw_shared_transitions
  public :: shared_cumhaz_covariance, shared_cumhaz_bootstrap_covariance
  public :: pmatrix_shared_ci

contains

  subroutine make_shared_msm(model,nstate,spec,fit,from,to,row,status)
    type(flexsurv_shared_msm),intent(out)::model
    integer,intent(in)::nstate,from(:),to(:),row(:)
    type(flexsurv_spec),intent(in)::spec
    type(flexsurv_result),intent(in)::fit
    integer,intent(out),optional::status
    integer::st
    st=0
    if(size(from)/=size(to).or.size(from)/=size(row))st=1
    if(.not.allocated(fit%theta))st=2
    if(st==0.and.size(fit%theta)/=parameter_count(spec))st=3
    if(st==0)then
      model%nstate=nstate;model%spec=spec;model%theta=fit%theta
      if(allocated(fit%covariance))model%covariance=fit%covariance
      model%from=from;model%to=to;model%row=row
    end if
    if(present(status))status=st
  end subroutine make_shared_msm

  subroutine shared_transitions(model,trans,theta)
    type(flexsurv_shared_msm),intent(in)::model
    type(flexsurv_transition),allocatable,intent(out)::trans(:)
    real(dp),intent(in),optional::theta(:)
    integer::j,n
    n=size(model%from);allocate(trans(n))
    do j=1,n
      trans(j)%from=model%from(j);trans(j)%to=model%to(j);trans(j)%row=model%row(j)
      trans(j)%model_kind=msm_parametric;trans(j)%spec=model%spec
      if(present(theta))then;trans(j)%theta=theta;else;trans(j)%theta=model%theta;end if
      if(allocated(model%covariance))trans(j)%covariance=model%covariance
    end do
  end subroutine shared_transitions

  subroutine draw_shared_transitions(model,trans)
    type(flexsurv_shared_msm),intent(in)::model
    type(flexsurv_transition),allocatable,intent(out)::trans(:)
    real(dp),allocatable::draw(:)
    if(allocated(model%covariance))then
      call mvn_draw_shared(model%theta,model%covariance,draw)
      call shared_transitions(model,trans,draw)
    else
      call shared_transitions(model,trans)
    end if
  end subroutine draw_shared_transitions

  subroutine shared_cumhaz_covariance(model,times,covhaz,status)
    type(flexsurv_shared_msm),intent(in)::model
    real(dp),intent(in)::times(:)
    real(dp),allocatable,intent(out)::covhaz(:,:,:)
    integer,intent(out),optional::status
    real(dp),allocatable::grad(:,:,:),xp(:),xm(:)
    real(dp)::h,hp,hm
    integer::nt,ntr,np,a,i,r,s,st
    nt=size(times);ntr=size(model%from);np=size(model%theta);st=0
    allocate(covhaz(ntr,ntr,nt));covhaz=0.0_dp
    if(.not.allocated(model%covariance))then;st=1;if(present(status))status=st;return;end if
    if(any(shape(model%covariance)/=[np,np]))then;st=2;if(present(status))status=st;return;end if
    allocate(grad(ntr,np,nt),xp(np),xm(np));grad=0.0_dp
    do a=1,np
      h=1.0e-5_dp*max(1.0_dp,abs(model%theta(a)));xp=model%theta;xm=model%theta
      xp(a)=xp(a)+h;xm(a)=xm(a)-h
      do r=1,ntr
        do i=1,nt
          hp=predict_cumhaz(model%spec,xp,model%row(r),times(i))
          hm=predict_cumhaz(model%spec,xm,model%row(r),times(i))
          grad(r,a,i)=(hp-hm)/(2.0_dp*h)
        end do
      end do
    end do
    do i=1,nt
      do r=1,ntr
        do s=r,ntr
          covhaz(r,s,i)=dot_product(grad(r,:,i),matmul(model%covariance,grad(s,:,i)))
          covhaz(s,r,i)=covhaz(r,s,i)
        end do
      end do
    end do
    if(present(status))status=st
  end subroutine shared_cumhaz_covariance

  subroutine shared_cumhaz_bootstrap_covariance(model,times,b,covhaz,seed,status)
    type(flexsurv_shared_msm),intent(in)::model
    real(dp),intent(in)::times(:)
    integer,intent(in)::b
    real(dp),allocatable,intent(out)::covhaz(:,:,:)
    integer,intent(in),optional::seed
    integer,intent(out),optional::status
    real(dp),allocatable::vals(:,:,:),draw(:),means(:,:)
    integer::bb,r,s,i,ntr,nt,st
    ntr=size(model%from);nt=size(times);st=0
    allocate(covhaz(ntr,ntr,nt));covhaz=0.0_dp
    if(.not.allocated(model%covariance).or.b<2)then;st=1;if(present(status))status=st;return;end if
    if(present(seed))call set_seed_shared(seed)
    allocate(vals(b,ntr,nt),means(ntr,nt))
    do bb=1,b
      call mvn_draw_shared(model%theta,model%covariance,draw)
      do r=1,ntr
        do i=1,nt;vals(bb,r,i)=predict_cumhaz(model%spec,draw,model%row(r),times(i));end do
      end do
    end do
    means=sum(vals,dim=1)/real(b,dp)
    do i=1,nt
      do r=1,ntr
        do s=r,ntr
          covhaz(r,s,i)=sum((vals(:,r,i)-means(r,i))*(vals(:,s,i)-means(s,i)))/real(b-1,dp)
          covhaz(s,r,i)=covhaz(r,s,i)
        end do
      end do
    end do
    if(present(status))status=st
  end subroutine shared_cumhaz_bootstrap_covariance

  subroutine pmatrix_shared_ci(model,times,b,lower,upper,meanp,cl,seed)
    type(flexsurv_shared_msm),intent(in)::model
    real(dp),intent(in)::times(:)
    integer,intent(in)::b
    real(dp),allocatable,intent(out)::lower(:,:,:),upper(:,:,:),meanp(:,:,:)
    real(dp),intent(in),optional::cl
    integer,intent(in),optional::seed
    real(dp),allocatable::boot(:,:,:,:),p(:,:,:)
    type(flexsurv_transition),allocatable::tr(:)
    real(dp)::cc,plo,phi
    integer::bb,i,j,k,n
    cc=0.95_dp;if(present(cl))cc=cl;plo=0.5_dp*(1.0_dp-cc);phi=1.0_dp-plo
    if(present(seed))call set_seed_shared(seed)
    n=model%nstate;allocate(boot(b,n,n,size(times)))
    do bb=1,b
      call draw_shared_transitions(model,tr);call pmatrix_flexsurv(tr,n,times,p);boot(bb,:,:,:)=p
    end do
    allocate(lower(n,n,size(times)),upper(n,n,size(times)),meanp(n,n,size(times)))
    do i=1,n;do j=1,n;do k=1,size(times)
      meanp(i,j,k)=sum(boot(:,i,j,k))/real(b,dp)
      lower(i,j,k)=sample_quantile(boot(:,i,j,k),plo);upper(i,j,k)=sample_quantile(boot(:,i,j,k),phi)
    end do;end do;end do
  end subroutine pmatrix_shared_ci

  subroutine mvn_draw_shared(mean,cov,draw)
    real(dp),intent(in)::mean(:),cov(:,:)
    real(dp),allocatable,intent(out)::draw(:)
    real(dp),allocatable::pd(:,:),l(:,:),z(:)
    integer::n,i,j,k
    real(dp)::v
    n=size(mean);allocate(pd(n,n),l(n,n),z(n),draw(n));call near_positive_definite(cov,pd,1.0e-10_dp)
    l=0.0_dp
    do i=1,n
      do j=1,i
        v=pd(i,j);do k=1,j-1;v=v-l(i,k)*l(j,k);end do
        if(i==j)then;l(i,j)=sqrt(max(v,1.0e-14_dp));else;l(i,j)=v/l(j,j);end if
      end do
    end do
    do i=1,n;z(i)=rng_normal();end do;draw=mean+matmul(l,z)
  end subroutine mvn_draw_shared

  real(dp) function sample_quantile(x,p) result(q)
    real(dp),intent(in)::x(:),p
    real(dp),allocatable::z(:)
    real(dp)::h,a,key
    integer::i,j,lo,hi,n
    z=x;n=size(z)
    do i=2,n
      key=z(i);j=i-1
      do while(j>=1)
        if(z(j)<=key)exit
        z(j+1)=z(j);j=j-1
      end do
      z(j+1)=key
    end do
    h=1.0_dp+real(n-1,dp)*p;lo=max(1,floor(h));hi=min(n,ceiling(h));a=h-real(lo,dp)
    q=(1.0_dp-a)*z(lo)+a*z(hi)
  end function sample_quantile

  subroutine set_seed_shared(seed)
    integer,intent(in)::seed
    integer::n,i
    integer,allocatable::put(:)
    call random_seed(size=n);allocate(put(n))
    do i=1,n;put(i)=mod(abs(seed)+19001*i,2147483646)+1;end do
    call random_seed(put=put)
  end subroutine set_seed_shared

end module flexsurv_shared_multistate
