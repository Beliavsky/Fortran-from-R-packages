! SPDX-License-Identifier: GPL-2.0-or-later
module flexsurv_final_states
  use flexsurv_kinds, only : dp
  use flexsurv_multistate, only : flexsurv_transition, msm_path, pmatrix_flexsurv
  use flexsurv_multistate_uncertainty, only : simulate_msm_path_reset, draw_transition_parameters
  implicit none
  private

  type, public :: final_state_summary
    integer, allocatable :: state(:)
    real(dp), allocatable :: probability(:)
    real(dp), allocatable :: mean_time(:)
    real(dp), allocatable :: quantile(:,:)
  end type final_state_summary

  public :: simfinal_fmsm, simfinal_fmsm_ci
  public :: pfinal_fmsm, pfinal_fmsm_ci

contains

  subroutine simfinal_fmsm(trans,nstate,start_state,tmax,m,probs,res,seed)
    type(flexsurv_transition),intent(in)::trans(:)
    integer,intent(in)::nstate,start_state,m
    real(dp),intent(in)::tmax,probs(:)
    type(final_state_summary),intent(out)::res
    integer,intent(in),optional::seed
    logical,allocatable::isabs(:)
    integer,allocatable::ast(:),cnt(:)
    real(dp),allocatable::vals(:,:)
    type(msm_path)::path
    integer::i,j,na,st,tot
    if(present(seed))call set_seed_final(seed)
    call absorbing_states(trans,nstate,isabs,ast)
    na=size(ast)
    allocate(res%state(na),res%probability(na),res%mean_time(na))
    allocate(res%quantile(na,size(probs)),cnt(na),vals(max(1,m),na))
    res%state=ast;res%probability=0.0_dp;res%mean_time=0.0_dp
    res%quantile=0.0_dp;cnt=0;vals=0.0_dp;tot=0
    do i=1,max(0,m)
      path=simulate_msm_path_reset(trans,nstate,start_state,tmax)
      st=path%state(path%n)
      if(st>=1.and.st<=nstate.and.isabs(st))then
        j=find_state(ast,st)
        if(j>0)then
          cnt(j)=cnt(j)+1;tot=tot+1
          vals(cnt(j),j)=path%time(path%n)
        end if
      end if
    end do
    if(tot>0)res%probability=real(cnt,dp)/real(tot,dp)
    do j=1,na
      if(cnt(j)>0)then
        res%mean_time(j)=sum(vals(1:cnt(j),j))/real(cnt(j),dp)
        call quantiles(vals(1:cnt(j),j),probs,res%quantile(j,:))
      end if
    end do
  end subroutine simfinal_fmsm

  subroutine simfinal_fmsm_ci(trans,nstate,start_state,tmax,m,probs,b,cl, &
      estimate,prob_lower,prob_upper,mean_lower,mean_upper,q_lower,q_upper,seed)
    type(flexsurv_transition),intent(in)::trans(:)
    integer,intent(in)::nstate,start_state,m,b
    real(dp),intent(in)::tmax,probs(:),cl
    type(final_state_summary),intent(out)::estimate
    real(dp),allocatable,intent(out)::prob_lower(:),prob_upper(:)
    real(dp),allocatable,intent(out)::mean_lower(:),mean_upper(:)
    real(dp),allocatable,intent(out)::q_lower(:,:),q_upper(:,:)
    integer,intent(in),optional::seed
    type(flexsurv_transition),allocatable::tb(:)
    type(final_state_summary)::rb
    real(dp),allocatable::pr(:,:),mn(:,:),qr(:,:,:)
    integer::ib,i,j,nb
    if(present(seed))call set_seed_final(seed)
    call simfinal_fmsm(trans,nstate,start_state,tmax,m,probs,estimate)
    nb=max(1,b)
    allocate(pr(nb,size(estimate%state)),mn(nb,size(estimate%state)))
    allocate(qr(nb,size(estimate%state),size(probs)))
    do ib=1,nb
      call draw_transition_parameters(trans,tb)
      call simfinal_fmsm(tb,nstate,start_state,tmax,m,probs,rb)
      pr(ib,:)=rb%probability;mn(ib,:)=rb%mean_time;qr(ib,:,:)=rb%quantile
    end do
    allocate(prob_lower(size(estimate%state)),prob_upper(size(estimate%state)))
    allocate(mean_lower(size(estimate%state)),mean_upper(size(estimate%state)))
    allocate(q_lower(size(estimate%state),size(probs)),q_upper(size(estimate%state),size(probs)))
    do i=1,size(estimate%state)
      prob_lower(i)=sample_quantile(pr(:,i),0.5_dp*(1.0_dp-cl))
      prob_upper(i)=sample_quantile(pr(:,i),1.0_dp-0.5_dp*(1.0_dp-cl))
      mean_lower(i)=sample_quantile(mn(:,i),0.5_dp*(1.0_dp-cl))
      mean_upper(i)=sample_quantile(mn(:,i),1.0_dp-0.5_dp*(1.0_dp-cl))
      do j=1,size(probs)
        q_lower(i,j)=sample_quantile(qr(:,i,j),0.5_dp*(1.0_dp-cl))
        q_upper(i,j)=sample_quantile(qr(:,i,j),1.0_dp-0.5_dp*(1.0_dp-cl))
      end do
    end do
  end subroutine simfinal_fmsm_ci

  subroutine pfinal_fmsm(trans,nstate,fromstate,maxt,prob,status)
    type(flexsurv_transition),intent(in)::trans(:)
    integer,intent(in)::nstate,fromstate
    real(dp),intent(in)::maxt
    real(dp),intent(out)::prob(nstate)
    integer,intent(out),optional::status
    real(dp),allocatable::p(:,:,:)
    integer::st
    call pmatrix_flexsurv(trans,nstate,[0.0_dp,maxt],p,st)
    prob=p(fromstate,:,2)
    if(present(status))status=st
  end subroutine pfinal_fmsm

  subroutine pfinal_fmsm_ci(trans,nstate,fromstate,maxt,b,cl,estimate,lower,upper,seed)
    type(flexsurv_transition),intent(in)::trans(:)
    integer,intent(in)::nstate,fromstate,b
    real(dp),intent(in)::maxt,cl
    real(dp),intent(out)::estimate(nstate),lower(nstate),upper(nstate)
    integer,intent(in),optional::seed
    type(flexsurv_transition),allocatable::tb(:)
    real(dp),allocatable::rep(:,:),v(:)
    integer::ib,j
    if(present(seed))call set_seed_final(seed)
    call pfinal_fmsm(trans,nstate,fromstate,maxt,estimate)
    allocate(rep(max(1,b),nstate),v(nstate))
    do ib=1,max(1,b)
      call draw_transition_parameters(trans,tb)
      call pfinal_fmsm(tb,nstate,fromstate,maxt,v)
      rep(ib,:)=v
    end do
    do j=1,nstate
      lower(j)=sample_quantile(rep(:,j),0.5_dp*(1.0_dp-cl))
      upper(j)=sample_quantile(rep(:,j),1.0_dp-0.5_dp*(1.0_dp-cl))
    end do
  end subroutine pfinal_fmsm_ci

  subroutine absorbing_states(trans,nstate,isabs,states)
    type(flexsurv_transition),intent(in)::trans(:)
    integer,intent(in)::nstate
    logical,allocatable,intent(out)::isabs(:)
    integer,allocatable,intent(out)::states(:)
    integer::i,n
    allocate(isabs(nstate));isabs=.true.
    do i=1,size(trans)
      if(trans(i)%from>=1.and.trans(i)%from<=nstate)isabs(trans(i)%from)=.false.
    end do
    n=count(isabs);allocate(states(n));n=0
    do i=1,nstate
      if(isabs(i))then;n=n+1;states(n)=i;end if
    end do
  end subroutine absorbing_states

  integer function find_state(states,state) result(idx)
    integer,intent(in)::states(:),state
    integer::i
    idx=0
    do i=1,size(states)
      if(states(i)==state)then;idx=i;return;end if
    end do
  end function find_state

  subroutine quantiles(x,p,q)
    real(dp),intent(in)::x(:),p(:)
    real(dp),intent(out)::q(size(p))
    integer::i
    do i=1,size(p);q(i)=sample_quantile(x,p(i));end do
  end subroutine quantiles

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
    k=max(1,min(n-1,int(floor(pos))));frac=pos-real(k,dp)
    q=(1.0_dp-frac)*a(k)+frac*a(k+1)
  end function sample_quantile

  subroutine set_seed_final(seed)
    integer,intent(in)::seed
    integer::n,i
    integer,allocatable::put(:)
    call random_seed(size=n);allocate(put(n))
    do i=1,n;put(i)=mod(abs(seed)+31337*i,2147483646)+1;end do
    call random_seed(put=put)
  end subroutine set_seed_final

end module flexsurv_final_states
