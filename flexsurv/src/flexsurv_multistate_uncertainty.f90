! SPDX-License-Identifier: GPL-2.0-or-later
module flexsurv_multistate_uncertainty
  use flexsurv_kinds, only : dp
  use flexsurv_fit, only : parameter_row, parameter_row_tcov
  use flexsurv_distributions, only : dist_random
  use flexsurv_spline, only : survspline_random, survspline_model
  use flexsurv_multistate, only : flexsurv_transition, msm_path, msm_parametric, msm_spline, &
    pmatrix_flexsurv, expected_length_of_stay
  use flexsurv_math, only : near_positive_definite, rng_normal
  implicit none
  private

  public :: simulate_msm_path_reset, pmatrix_simfs, totlos_simfs
  public :: pmatrix_simfs_ci, pmatrix_flexsurv_ci, totlos_simfs_ci
  public :: draw_transition_parameters, transition_parameters_at_time
  public :: transition_spline_model_at_time

contains

  function simulate_msm_path_reset(trans,nstate,start_state,tmax,max_events) result(path)
    type(flexsurv_transition),intent(in)::trans(:)
    integer,intent(in)::nstate,start_state
    real(dp),intent(in)::tmax
    integer,intent(in),optional::max_events
    type(msm_path)::path
    real(dp),allocatable::tt(:),wait(:)
    integer,allocatable::ss(:),idx(:)
    integer::m,nev,state,j,nout,k,best
    real(dp)::cur,wmin
    m=1000;if(present(max_events))m=max_events
    allocate(tt(m+1),ss(m+1));nev=1;tt(1)=0.0_dp;ss(1)=start_state
    cur=0.0_dp;state=start_state
    do while(cur<tmax.and.nev<=m)
      nout=count([(trans(j)%from==state,j=1,size(trans))])
      if(nout==0)exit
      allocate(idx(nout),wait(nout));k=0
      do j=1,size(trans)
        if(trans(j)%from==state)then
          k=k+1;idx(k)=j;wait(k)=transition_random(trans(j),cur)
          if(wait(k)<0.0_dp.or.wait(k)/=wait(k))wait(k)=huge(1.0_dp)
        end if
      end do
      best=1;wmin=wait(1)
      do k=2,nout
        if(wait(k)<wmin)then;wmin=wait(k);best=k;end if
      end do
      if(.not.(wmin<huge(1.0_dp)).or.cur+wmin>tmax)then
        deallocate(idx,wait);exit
      end if
      cur=cur+wmin;state=trans(idx(best))%to
      nev=nev+1;tt(nev)=cur;ss(nev)=state
      deallocate(idx,wait)
    end do
    path%n=nev;allocate(path%time(nev),path%state(nev))
    path%time=tt(1:nev);path%state=ss(1:nev)
  end function simulate_msm_path_reset

  subroutine pmatrix_simfs(trans,nstate,times,m,p,seed)
    type(flexsurv_transition),intent(in)::trans(:)
    integer,intent(in)::nstate,m
    real(dp),intent(in)::times(:)
    real(dp),allocatable,intent(out)::p(:,:,:)
    integer,intent(in),optional::seed
    type(msm_path)::path
    integer::i,j,k,r,st
    real(dp)::tmax
    if(present(seed))call set_seed(seed)
    allocate(p(nstate,nstate,size(times)));p=0.0_dp
    if(size(times)==0.or.m<=0)return
    tmax=maxval(times)
    do i=1,nstate
      do r=1,m
        path=simulate_msm_path_reset(trans,nstate,i,tmax)
        do k=1,size(times)
          st=path%state(1)
          do j=2,path%n
            if(path%time(j)<=times(k))then;st=path%state(j);else;exit;end if
          end do
          if(st>=1.and.st<=nstate)p(i,st,k)=p(i,st,k)+1.0_dp
        end do
      end do
    end do
    p=p/real(m,dp)
  end subroutine pmatrix_simfs

  subroutine totlos_simfs(trans,nstate,start_state,tmax,m,elos,seed)
    type(flexsurv_transition),intent(in)::trans(:)
    integer,intent(in)::nstate,start_state,m
    real(dp),intent(in)::tmax
    real(dp),intent(out)::elos(nstate)
    integer,intent(in),optional::seed
    type(msm_path)::path
    integer::r,j,st
    real(dp)::a,b
    if(present(seed))call set_seed(seed)
    elos=0.0_dp
    do r=1,m
      path=simulate_msm_path_reset(trans,nstate,start_state,tmax)
      do j=1,path%n
        st=path%state(j);a=path%time(j)
        if(j<path%n)then;b=min(tmax,path%time(j+1));else;b=tmax;end if
        if(st>=1.and.st<=nstate.and.b>a)elos(st)=elos(st)+(b-a)
      end do
    end do
    if(m>0)elos=elos/real(m,dp)
  end subroutine totlos_simfs

  subroutine pmatrix_simfs_ci(trans,nstate,times,m,b,cl,estimate,lower,upper,seed)
    type(flexsurv_transition),intent(in)::trans(:)
    integer,intent(in)::nstate,m,b
    real(dp),intent(in)::times(:),cl
    real(dp),allocatable,intent(out)::estimate(:,:,:),lower(:,:,:),upper(:,:,:)
    integer,intent(in),optional::seed
    type(flexsurv_transition),allocatable::tb(:)
    real(dp),allocatable::pb(:,:,:),rep(:,:,:,:)
    integer::ib,s0
    s0=13579;if(present(seed))s0=seed;call set_seed(s0)
    call pmatrix_simfs(trans,nstate,times,m,estimate)
    allocate(rep(b,nstate,nstate,size(times)))
    do ib=1,b
      call draw_transition_parameters(trans,tb)
      call pmatrix_simfs(tb,nstate,times,m,pb)
      rep(ib,:,:,:)=pb
    end do
    call quantile_array4(rep,0.5_dp*(1.0_dp-cl),lower)
    call quantile_array4(rep,1.0_dp-0.5_dp*(1.0_dp-cl),upper)
  end subroutine pmatrix_simfs_ci

  subroutine pmatrix_flexsurv_ci(trans,nstate,times,b,cl,estimate,lower,upper,seed)
    type(flexsurv_transition),intent(in)::trans(:)
    integer,intent(in)::nstate,b
    real(dp),intent(in)::times(:),cl
    real(dp),allocatable,intent(out)::estimate(:,:,:),lower(:,:,:),upper(:,:,:)
    integer,intent(in),optional::seed
    type(flexsurv_transition),allocatable::tb(:)
    real(dp),allocatable::pb(:,:,:),rep(:,:,:,:)
    integer::ib,s0
    s0=24681;if(present(seed))s0=seed;call set_seed(s0)
    call pmatrix_flexsurv(trans,nstate,times,estimate)
    allocate(rep(b,nstate,nstate,size(times)))
    do ib=1,b
      call draw_transition_parameters(trans,tb)
      call pmatrix_flexsurv(tb,nstate,times,pb)
      rep(ib,:,:,:)=pb
    end do
    call quantile_array4(rep,0.5_dp*(1.0_dp-cl),lower)
    call quantile_array4(rep,1.0_dp-0.5_dp*(1.0_dp-cl),upper)
  end subroutine pmatrix_flexsurv_ci

  subroutine totlos_simfs_ci(trans,nstate,start_state,tmax,m,b,cl,estimate,lower,upper,seed)
    type(flexsurv_transition),intent(in)::trans(:)
    integer,intent(in)::nstate,start_state,m,b
    real(dp),intent(in)::tmax,cl
    real(dp),intent(out)::estimate(nstate),lower(nstate),upper(nstate)
    integer,intent(in),optional::seed
    type(flexsurv_transition),allocatable::tb(:)
    real(dp),allocatable::rep(:,:),eb(:)
    integer::ib,j,s0
    s0=97531;if(present(seed))s0=seed;call set_seed(s0)
    call totlos_simfs(trans,nstate,start_state,tmax,m,estimate)
    allocate(rep(b,nstate),eb(nstate))
    do ib=1,b
      call draw_transition_parameters(trans,tb)
      call totlos_simfs(tb,nstate,start_state,tmax,m,eb)
      rep(ib,:)=eb
    end do
    do j=1,nstate
      lower(j)=sample_quantile(rep(:,j),0.5_dp*(1.0_dp-cl))
      upper(j)=sample_quantile(rep(:,j),1.0_dp-0.5_dp*(1.0_dp-cl))
    end do
  end subroutine totlos_simfs_ci

  subroutine draw_transition_parameters(trans,out)
    type(flexsurv_transition),intent(in)::trans(:)
    type(flexsurv_transition),allocatable,intent(out)::out(:)
    real(dp),allocatable::draw(:)
    integer::j,ng
    allocate(out(size(trans)));out=trans
    do j=1,size(trans)
      select case(trans(j)%model_kind)
      case(msm_parametric)
        if(allocated(trans(j)%theta).and.allocated(trans(j)%covariance))then
          call mvn_draw(trans(j)%theta,trans(j)%covariance,draw)
          out(j)%theta=draw
        end if
      case(msm_spline)
        ng=size(trans(j)%spline%model%gamma)
        if(allocated(trans(j)%spline%covariance))then
          if(size(trans(j)%spline%covariance,1)>=ng)then
            call mvn_draw(trans(j)%spline%model%gamma, &
              trans(j)%spline%covariance(1:ng,1:ng),draw)
            out(j)%spline%model%gamma=draw
          end if
        end if
      end select
    end do
  end subroutine draw_transition_parameters

  real(dp) function transition_random(tr,current_time) result(v)
    type(flexsurv_transition),intent(in)::tr
    real(dp),intent(in),optional::current_time
    real(dp),allocatable::par(:)
    select case(tr%model_kind)
    case(msm_parametric)
      if(present(current_time).and.allocated(tr%tcov_parameter).and. &
          allocated(tr%tcov_column).and.allocated(tr%tcov_rate))then
        par=parameter_row_tcov(tr%spec,tr%theta,tr%row,current_time, &
          tr%tcov_parameter,tr%tcov_column,tr%tcov_rate)
      else
        par=parameter_row(tr%spec,tr%theta,tr%row)
      end if
      v=dist_random(tr%spec%dist,par)
    case(msm_spline)
      if(present(current_time))then
        v=survspline_random(transition_spline_model_at_time(tr,current_time))
      else
        v=survspline_random(tr%spline%model)
      end if
    case default
      v=huge(1.0_dp)
    end select
  end function transition_random

  function transition_parameters_at_time(tr,current_time) result(par)
    type(flexsurv_transition),intent(in)::tr
    real(dp),intent(in)::current_time
    real(dp),allocatable::par(:)
    if(tr%model_kind/=msm_parametric)then
      allocate(par(0));return
    end if
    if(allocated(tr%tcov_parameter).and.allocated(tr%tcov_column).and. &
        allocated(tr%tcov_rate))then
      par=parameter_row_tcov(tr%spec,tr%theta,tr%row,current_time, &
        tr%tcov_parameter,tr%tcov_column,tr%tcov_rate)
    else
      par=parameter_row(tr%spec,tr%theta,tr%row)
    end if
  end function transition_parameters_at_time


  function transition_spline_model_at_time(tr,current_time) result(model)
    type(flexsurv_transition),intent(in)::tr
    real(dp),intent(in)::current_time
    type(survspline_model)::model
    model=tr%spline%model
    if(allocated(tr%spline_gamma_tcov_rate))then
      if(size(tr%spline_gamma_tcov_rate)==size(model%gamma)) &
        model%gamma=model%gamma+current_time*tr%spline_gamma_tcov_rate
    end if
  end function transition_spline_model_at_time

  subroutine mvn_draw(mean,cov,draw)
    real(dp),intent(in)::mean(:),cov(:,:)
    real(dp),allocatable,intent(out)::draw(:)
    real(dp),allocatable::pd(:,:),l(:,:),z(:)
    integer::n,i,st
    n=size(mean);allocate(pd(n,n),l(n,n),z(n),draw(n))
    call near_positive_definite(cov,pd,1.0e-10_dp);call chol_lower(pd,l,st)
    if(st/=0)then;draw=mean;return;end if
    do i=1,n;z(i)=rng_normal();end do
    draw=mean+matmul(l,z)
  end subroutine mvn_draw

  subroutine chol_lower(a,l,status)
    real(dp),intent(in)::a(:,:)
    real(dp),intent(out)::l(size(a,1),size(a,2))
    integer,intent(out)::status
    real(dp)::s
    integer::i,j,k,n
    n=size(a,1);l=0.0_dp;status=0
    do i=1,n
      do j=1,i
        s=a(i,j)
        do k=1,j-1;s=s-l(i,k)*l(j,k);end do
        if(i==j)then
          if(s<=0.0_dp)then;status=1;return;end if
          l(i,j)=sqrt(s)
        else
          l(i,j)=s/l(j,j)
        end if
      end do
    end do
  end subroutine chol_lower

  subroutine quantile_array4(rep,p,out)
    real(dp),intent(in)::rep(:,:,:,:),p
    real(dp),allocatable,intent(out)::out(:,:,:)
    integer::i,j,k
    allocate(out(size(rep,2),size(rep,3),size(rep,4)))
    do k=1,size(rep,4);do j=1,size(rep,3);do i=1,size(rep,2)
      out(i,j,k)=sample_quantile(rep(:,i,j,k),p)
    end do;end do;end do
  end subroutine quantile_array4

  real(dp) function sample_quantile(x,p) result(q)
    real(dp),intent(in)::x(:),p
    real(dp),allocatable::a(:)
    real(dp)::key,pos,frac
    integer::i,j,k,n
    n=size(x);if(n==0)then;q=0.0_dp;return;end if
    allocate(a(n));a=x
    do i=2,n
      key=a(i);j=i-1
      do while(j>=1)
        if(a(j)<=key)exit
        a(j+1)=a(j);j=j-1
      end do
      a(j+1)=key
    end do
    if(n==1)then;q=a(1);return;end if
    pos=1.0_dp+max(0.0_dp,min(1.0_dp,p))*real(n-1,dp)
    k=min(n-1,max(1,int(floor(pos))));frac=pos-real(k,dp)
    q=(1.0_dp-frac)*a(k)+frac*a(k+1)
  end function sample_quantile

  subroutine set_seed(seed)
    integer,intent(in)::seed
    integer::n,i
    integer,allocatable::put(:)
    call random_seed(size=n);allocate(put(n))
    do i=1,n;put(i)=mod(abs(seed)+104729*i,2147483646)+1;end do
    call random_seed(put=put)
  end subroutine set_seed

end module flexsurv_multistate_uncertainty
