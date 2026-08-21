! SPDX-License-Identifier: GPL-2.0-or-later
module flexsurv_fmixmsm
  use flexsurv_kinds, only : dp
  use flexsurv_mixture_full, only : flexsurvmix_full_result, mix_probability_at, &
    mix_component_mean, mix_component_random, resample_flexsurvmix_full
  implicit none
  private

  type, public :: fmix_transition_model
    integer :: from = 0
    integer, allocatable :: to(:)
    type(flexsurvmix_full_result) :: mix
  end type fmix_transition_model

  type, public :: fmix_pathways
    integer, allocatable :: state(:,:)
    integer, allocatable :: length(:)
    integer :: npath = 0
    logical :: has_cycle = .false.
  end type fmix_pathways

  type, public :: fmixmsm_model
    integer :: nstate = 0
    integer :: start_state = 1
    type(fmix_transition_model), allocatable :: node(:)
    type(fmix_pathways) :: pathways
  end type fmixmsm_model

  public :: build_fmixmsm, enumerate_fmix_pathways
  public :: fmix_path_probabilities, fmix_final_probabilities
  public :: fmix_path_mean_times, fmix_final_mean_times
  public :: fmix_path_time_quantiles, fmix_final_time_quantiles
  public :: fmix_path_probabilities_ci, fmix_final_probabilities_ci
  public :: fmix_path_mean_times_ci, fmix_final_mean_times_ci
  public :: fmix_path_time_quantiles_ci, fmix_final_time_quantiles_ci
  public :: resample_fmixmsm

contains

  function build_fmixmsm(node,nstate,start_state) result(model)
    type(fmix_transition_model),intent(in)::node(:)
    integer,intent(in)::nstate
    integer,intent(in),optional::start_state
    type(fmixmsm_model)::model
    model%nstate=nstate;model%start_state=1;if(present(start_state))model%start_state=start_state
    allocate(model%node(size(node)));model%node=node
    call enumerate_fmix_pathways(model,model%pathways)
  end function build_fmixmsm

  subroutine enumerate_fmix_pathways(model,paths)
    type(fmixmsm_model),intent(in)::model
    type(fmix_pathways),intent(out)::paths
    integer,allocatable::tmp(:,:),lens(:),cur(:)
    integer::cap,nfound
    cap=max(16,model%nstate*model%nstate);allocate(tmp(cap,model%nstate+1),lens(cap),cur(model%nstate+1))
    tmp=0;lens=0;cur=0;nfound=0;paths%has_cycle=.false.;cur(1)=model%start_state
    call dfs(model%start_state,1)
    paths%npath=nfound
    allocate(paths%state(max(0,nfound),model%nstate+1),paths%length(max(0,nfound)))
    if(nfound>0)then;paths%state=tmp(1:nfound,:);paths%length=lens(1:nfound);end if
  contains
    recursive subroutine dfs(state,depth)
      integer,intent(in)::state,depth
      integer::idx,j,next
      idx=node_index(model,state)
      if(idx==0)then
        call add_path(depth);return
      end if
      do j=1,size(model%node(idx)%to)
        next=model%node(idx)%to(j)
        if(any(cur(1:depth)==next))then
          paths%has_cycle=.true.;return
        end if
        if(depth>=size(cur))then;paths%has_cycle=.true.;return;end if
        cur(depth+1)=next;call dfs(next,depth+1);cur(depth+1)=0
        if(paths%has_cycle)return
      end do
    end subroutine dfs
    subroutine add_path(depth)
      integer,intent(in)::depth
      integer,allocatable::nt(:,:),nl(:)
      if(nfound==cap)then
        allocate(nt(2*cap,model%nstate+1),nl(2*cap));nt=0;nl=0
        nt(1:cap,:)=tmp;nl(1:cap)=lens;call move_alloc(nt,tmp);call move_alloc(nl,lens);cap=2*cap
      end if
      nfound=nfound+1;tmp(nfound,1:depth)=cur(1:depth);lens(nfound)=depth
    end subroutine add_path
  end subroutine enumerate_fmix_pathways

  subroutine fmix_path_probabilities(model,row,prob)
    type(fmixmsm_model),intent(in)::model
    integer,intent(in)::row
    real(dp),intent(out)::prob(model%pathways%npath)
    real(dp),allocatable::pr(:)
    integer::ip,j,idx,k,next
    prob=1.0_dp
    do ip=1,model%pathways%npath
      do j=1,model%pathways%length(ip)-1
        idx=node_index(model,model%pathways%state(ip,j));next=model%pathways%state(ip,j+1)
        if(idx==0)then;prob(ip)=0.0_dp;exit;end if
        allocate(pr(model%node(idx)%mix%k));call row_prob(model%node(idx)%mix,row,pr)
        k=find_destination(model%node(idx),next)
        if(k==0)then;prob(ip)=0.0_dp;deallocate(pr);exit;end if
        prob(ip)=prob(ip)*pr(k);deallocate(pr)
      end do
    end do
  end subroutine fmix_path_probabilities

  subroutine fmix_final_probabilities(model,row,prob)
    type(fmixmsm_model),intent(in)::model
    integer,intent(in)::row
    real(dp),intent(out)::prob(model%nstate)
    real(dp),allocatable::pp(:)
    integer::i,st
    allocate(pp(model%pathways%npath));call fmix_path_probabilities(model,row,pp);prob=0.0_dp
    do i=1,model%pathways%npath
      st=model%pathways%state(i,model%pathways%length(i));prob(st)=prob(st)+pp(i)
    end do
  end subroutine fmix_final_probabilities

  subroutine fmix_path_mean_times(model,row,means)
    type(fmixmsm_model),intent(in)::model
    integer,intent(in)::row
    real(dp),intent(out)::means(model%pathways%npath)
    integer::ip,j,idx,k,next
    means=0.0_dp
    do ip=1,model%pathways%npath
      do j=1,model%pathways%length(ip)-1
        idx=node_index(model,model%pathways%state(ip,j));next=model%pathways%state(ip,j+1)
        k=find_destination(model%node(idx),next)
        means(ip)=means(ip)+mix_component_mean(model%node(idx)%mix,k,row)
      end do
    end do
  end subroutine fmix_path_mean_times

  subroutine fmix_final_mean_times(model,row,means)
    type(fmixmsm_model),intent(in)::model
    integer,intent(in)::row
    real(dp),intent(out)::means(model%nstate)
    real(dp),allocatable::pp(:),mp(:),den(:)
    integer::i,st
    allocate(pp(model%pathways%npath),mp(model%pathways%npath),den(model%nstate))
    call fmix_path_probabilities(model,row,pp);call fmix_path_mean_times(model,row,mp)
    means=0.0_dp;den=0.0_dp
    do i=1,model%pathways%npath
      st=model%pathways%state(i,model%pathways%length(i))
      means(st)=means(st)+pp(i)*mp(i);den(st)=den(st)+pp(i)
    end do
    do st=1,model%nstate;if(den(st)>0.0_dp)means(st)=means(st)/den(st);end do
  end subroutine fmix_final_mean_times

  subroutine fmix_path_time_quantiles(model,row,nsim,probs,q,seed)
    type(fmixmsm_model),intent(in)::model
    integer,intent(in)::row,nsim
    real(dp),intent(in)::probs(:)
    real(dp),intent(out)::q(model%pathways%npath,size(probs))
    integer,intent(in),optional::seed
    real(dp),allocatable::s(:)
    integer::ip,m,j,idx,k,next
    if(present(seed))call set_seed_local(seed)
    allocate(s(nsim))
    do ip=1,model%pathways%npath
      do m=1,nsim
        s(m)=0.0_dp
        do j=1,model%pathways%length(ip)-1
          idx=node_index(model,model%pathways%state(ip,j));next=model%pathways%state(ip,j+1)
          k=find_destination(model%node(idx),next)
          s(m)=s(m)+mix_component_random(model%node(idx)%mix,k,row)
        end do
      end do
      call sample_quantiles(s,probs,q(ip,:))
    end do
  end subroutine fmix_path_time_quantiles

  subroutine fmix_final_time_quantiles(model,row,nsim,probs,q,seed)
    type(fmixmsm_model),intent(in)::model
    integer,intent(in)::row,nsim
    real(dp),intent(in)::probs(:)
    real(dp),intent(out)::q(model%nstate,size(probs))
    integer,intent(in),optional::seed
    real(dp),allocatable::pp(:),work(:,:),u(:)
    integer,allocatable::nby(:)
    integer::st,ip,m,j,idx,k,next,np,sel
    real(dp)::den,r,cum
    if(present(seed))call set_seed_local(seed)
    allocate(pp(model%pathways%npath),work(nsim,model%nstate),nby(model%nstate));work=0.0_dp;nby=0
    call fmix_path_probabilities(model,row,pp)
    do st=1,model%nstate
      den=0.0_dp
      do ip=1,model%pathways%npath
        if(model%pathways%state(ip,model%pathways%length(ip))==st)den=den+pp(ip)
      end do
      if(den<=0.0_dp)cycle
      do m=1,nsim
        call random_number(r);cum=0.0_dp;sel=0
        do ip=1,model%pathways%npath
          if(model%pathways%state(ip,model%pathways%length(ip))/=st)cycle
          cum=cum+pp(ip)/den
          if(r<=cum)then;sel=ip;exit;end if
        end do
        if(sel==0)cycle
        nby(st)=nby(st)+1
        do j=1,model%pathways%length(sel)-1
          idx=node_index(model,model%pathways%state(sel,j));next=model%pathways%state(sel,j+1)
          k=find_destination(model%node(idx),next)
          work(nby(st),st)=work(nby(st),st)+mix_component_random(model%node(idx)%mix,k,row)
        end do
      end do
      if(nby(st)>0)call sample_quantiles(work(1:nby(st),st),probs,q(st,:))
    end do
  end subroutine fmix_final_time_quantiles


  subroutine resample_fmixmsm(model,out,seed)
    type(fmixmsm_model),intent(in)::model
    type(fmixmsm_model),intent(out)::out
    integer,intent(in),optional::seed
    integer::j
    out=model
    if(present(seed))call set_seed_local(seed)
    do j=1,size(model%node)
      call resample_flexsurvmix_full(model%node(j)%mix,out%node(j)%mix)
    end do
  end subroutine resample_fmixmsm

  subroutine fmix_path_probabilities_ci(model,row,b,cl,estimate,lower,upper,seed)
    type(fmixmsm_model),intent(in)::model
    integer,intent(in)::row,b
    real(dp),intent(in)::cl
    real(dp),intent(out)::estimate(model%pathways%npath)
    real(dp),intent(out)::lower(model%pathways%npath),upper(model%pathways%npath)
    integer,intent(in),optional::seed
    type(fmixmsm_model)::mb
    real(dp),allocatable::rep(:,:),v(:)
    integer::ib,j
    if(present(seed))call set_seed_local(seed)
    call fmix_path_probabilities(model,row,estimate)
    allocate(rep(max(1,b),size(estimate)),v(size(estimate)))
    do ib=1,max(1,b)
      call resample_fmixmsm(model,mb)
      call fmix_path_probabilities(mb,row,v)
      rep(ib,:)=v
    end do
    do j=1,size(estimate)
      lower(j)=sample_quantile1(rep(:,j),0.5_dp*(1.0_dp-cl))
      upper(j)=sample_quantile1(rep(:,j),1.0_dp-0.5_dp*(1.0_dp-cl))
    end do
  end subroutine fmix_path_probabilities_ci

  subroutine fmix_final_probabilities_ci(model,row,b,cl,estimate,lower,upper,seed)
    type(fmixmsm_model),intent(in)::model
    integer,intent(in)::row,b
    real(dp),intent(in)::cl
    real(dp),intent(out)::estimate(model%nstate),lower(model%nstate),upper(model%nstate)
    integer,intent(in),optional::seed
    type(fmixmsm_model)::mb
    real(dp),allocatable::rep(:,:),v(:)
    integer::ib,j
    if(present(seed))call set_seed_local(seed)
    call fmix_final_probabilities(model,row,estimate)
    allocate(rep(max(1,b),model%nstate),v(model%nstate))
    do ib=1,max(1,b)
      call resample_fmixmsm(model,mb)
      call fmix_final_probabilities(mb,row,v)
      rep(ib,:)=v
    end do
    do j=1,model%nstate
      lower(j)=sample_quantile1(rep(:,j),0.5_dp*(1.0_dp-cl))
      upper(j)=sample_quantile1(rep(:,j),1.0_dp-0.5_dp*(1.0_dp-cl))
    end do
  end subroutine fmix_final_probabilities_ci

  subroutine fmix_path_mean_times_ci(model,row,b,cl,estimate,lower,upper,seed)
    type(fmixmsm_model),intent(in)::model
    integer,intent(in)::row,b
    real(dp),intent(in)::cl
    real(dp),intent(out)::estimate(model%pathways%npath)
    real(dp),intent(out)::lower(model%pathways%npath),upper(model%pathways%npath)
    integer,intent(in),optional::seed
    type(fmixmsm_model)::mb
    real(dp),allocatable::rep(:,:),v(:)
    integer::ib,j
    if(present(seed))call set_seed_local(seed)
    call fmix_path_mean_times(model,row,estimate)
    allocate(rep(max(1,b),size(estimate)),v(size(estimate)))
    do ib=1,max(1,b)
      call resample_fmixmsm(model,mb)
      call fmix_path_mean_times(mb,row,v)
      rep(ib,:)=v
    end do
    do j=1,size(estimate)
      lower(j)=sample_quantile1(rep(:,j),0.5_dp*(1.0_dp-cl))
      upper(j)=sample_quantile1(rep(:,j),1.0_dp-0.5_dp*(1.0_dp-cl))
    end do
  end subroutine fmix_path_mean_times_ci

  subroutine fmix_final_mean_times_ci(model,row,b,cl,estimate,lower,upper,seed)
    type(fmixmsm_model),intent(in)::model
    integer,intent(in)::row,b
    real(dp),intent(in)::cl
    real(dp),intent(out)::estimate(model%nstate),lower(model%nstate),upper(model%nstate)
    integer,intent(in),optional::seed
    type(fmixmsm_model)::mb
    real(dp),allocatable::rep(:,:),v(:)
    integer::ib,j
    if(present(seed))call set_seed_local(seed)
    call fmix_final_mean_times(model,row,estimate)
    allocate(rep(max(1,b),model%nstate),v(model%nstate))
    do ib=1,max(1,b)
      call resample_fmixmsm(model,mb)
      call fmix_final_mean_times(mb,row,v)
      rep(ib,:)=v
    end do
    do j=1,model%nstate
      lower(j)=sample_quantile1(rep(:,j),0.5_dp*(1.0_dp-cl))
      upper(j)=sample_quantile1(rep(:,j),1.0_dp-0.5_dp*(1.0_dp-cl))
    end do
  end subroutine fmix_final_mean_times_ci

  subroutine fmix_path_time_quantiles_ci(model,row,nsim,probs,b,cl,estimate,lower,upper,seed)
    type(fmixmsm_model),intent(in)::model
    integer,intent(in)::row,nsim,b
    real(dp),intent(in)::probs(:),cl
    real(dp),intent(out)::estimate(model%pathways%npath,size(probs))
    real(dp),intent(out)::lower(model%pathways%npath,size(probs))
    real(dp),intent(out)::upper(model%pathways%npath,size(probs))
    integer,intent(in),optional::seed
    type(fmixmsm_model)::mb
    real(dp),allocatable::rep(:,:,:),v(:,:)
    integer::ib,i,j
    if(present(seed))call set_seed_local(seed)
    call fmix_path_time_quantiles(model,row,nsim,probs,estimate)
    allocate(rep(max(1,b),size(estimate,1),size(estimate,2)))
    allocate(v(size(estimate,1),size(estimate,2)))
    do ib=1,max(1,b)
      call resample_fmixmsm(model,mb)
      call fmix_path_time_quantiles(mb,row,nsim,probs,v)
      rep(ib,:,:)=v
    end do
    do j=1,size(estimate,2)
      do i=1,size(estimate,1)
        lower(i,j)=sample_quantile1(rep(:,i,j),0.5_dp*(1.0_dp-cl))
        upper(i,j)=sample_quantile1(rep(:,i,j),1.0_dp-0.5_dp*(1.0_dp-cl))
      end do
    end do
  end subroutine fmix_path_time_quantiles_ci

  subroutine fmix_final_time_quantiles_ci(model,row,nsim,probs,b,cl,estimate,lower,upper,seed)
    type(fmixmsm_model),intent(in)::model
    integer,intent(in)::row,nsim,b
    real(dp),intent(in)::probs(:),cl
    real(dp),intent(out)::estimate(model%nstate,size(probs))
    real(dp),intent(out)::lower(model%nstate,size(probs)),upper(model%nstate,size(probs))
    integer,intent(in),optional::seed
    type(fmixmsm_model)::mb
    real(dp),allocatable::rep(:,:,:),v(:,:)
    integer::ib,i,j
    if(present(seed))call set_seed_local(seed)
    call fmix_final_time_quantiles(model,row,nsim,probs,estimate)
    allocate(rep(max(1,b),model%nstate,size(probs)),v(model%nstate,size(probs)))
    do ib=1,max(1,b)
      call resample_fmixmsm(model,mb)
      call fmix_final_time_quantiles(mb,row,nsim,probs,v)
      rep(ib,:,:)=v
    end do
    do j=1,size(probs)
      do i=1,model%nstate
        lower(i,j)=sample_quantile1(rep(:,i,j),0.5_dp*(1.0_dp-cl))
        upper(i,j)=sample_quantile1(rep(:,i,j),1.0_dp-0.5_dp*(1.0_dp-cl))
      end do
    end do
  end subroutine fmix_final_time_quantiles_ci

  real(dp) function sample_quantile1(x,p) result(q)
    real(dp),intent(in)::x(:),p
    real(dp),allocatable::y(:)
    real(dp)::pos,frac
    integer::i,k,n
    n=size(x);if(n==0)then;q=0.0_dp;return;end if
    allocate(y(n));y=x;call sort_real(y)
    if(n==1)then;q=y(1);return;end if
    pos=1.0_dp+max(0.0_dp,min(1.0_dp,p))*real(n-1,dp)
    i=max(1,min(n-1,int(floor(pos))));k=i+1;frac=pos-real(i,dp)
    q=(1.0_dp-frac)*y(i)+frac*y(k)
  end function sample_quantile1

  integer function node_index(model,state) result(idx)
    type(fmixmsm_model),intent(in)::model
    integer,intent(in)::state
    integer::i
    idx=0;do i=1,size(model%node);if(model%node(i)%from==state)then;idx=i;return;end if;end do
  end function node_index

  integer function find_destination(node,state) result(idx)
    type(fmix_transition_model),intent(in)::node
    integer,intent(in)::state
    integer::i
    idx=0;do i=1,size(node%to);if(node%to(i)==state)then;idx=i;return;end if;end do
  end function find_destination

  subroutine row_prob(mix,row,pr)
    type(flexsurvmix_full_result),intent(in)::mix
    integer,intent(in)::row
    real(dp),intent(out)::pr(mix%k)
    real(dp),allocatable::x(:)
    if(allocated(mix%prob_x).and.row>=1.and.row<=size(mix%prob_x,1))then
      allocate(x(size(mix%prob_x,2)));x=mix%prob_x(row,:);call mix_probability_at(mix,x,pr)
    else
      allocate(x(mix%nprob_cov));x=0.0_dp;call mix_probability_at(mix,x,pr)
    end if
  end subroutine row_prob

  subroutine sample_quantiles(x,p,q)
    real(dp),intent(in)::x(:),p(:)
    real(dp),intent(out)::q(size(p))
    real(dp),allocatable::y(:)
    real(dp)::h,frac
    integer::i,j,k,n
    n=size(x);allocate(y(n));y=x;call sort_real(y)
    do j=1,size(p)
      if(p(j)<=0.0_dp)then;q(j)=y(1)
      else if(p(j)>=1.0_dp)then;q(j)=y(n)
      else
        h=1.0_dp+real(n-1,dp)*p(j);i=floor(h);k=min(n,i+1);frac=h-real(i,dp)
        q(j)=(1.0_dp-frac)*y(i)+frac*y(k)
      end if
    end do
  end subroutine sample_quantiles

  subroutine sort_real(x)
    real(dp),intent(inout)::x(:)
    integer::i,j
    real(dp)::key
    do i=2,size(x)
      key=x(i);j=i-1
      do while(j>=1)
        if(x(j)<=key)exit
        x(j+1)=x(j);j=j-1
      end do
      x(j+1)=key
    end do
  end subroutine sort_real

  subroutine set_seed_local(seed)
    integer,intent(in)::seed
    integer::n,i
    integer,allocatable::put(:)
    call random_seed(size=n);allocate(put(n));do i=1,n;put(i)=mod(abs(seed)+7919*i,2147483646)+1;end do
    call random_seed(put=put)
  end subroutine set_seed_local

end module flexsurv_fmixmsm
