module flsss_gap
  use flsss_kinds, only : dp, i8
  use flsss_types, only : gap_result
  use flsss_util, only : timer_type
  implicit none
  private
  public :: gap_solve, aux_gap_bb, aux_gap_bbdp, aux_gap_ga

contains

  function gap_solve(cost, profit, budget, tlimit, heuristic) result(r)
    real(dp), intent(in) :: cost(:,:), profit(:,:), budget(:)
    real(dp), intent(in), optional :: tlimit
    logical, intent(in), optional :: heuristic
    type(gap_result) :: r
    real(dp) :: lim
    if (present(heuristic)) then
      if (heuristic .and. size(cost,2) < 0) error stop "gap_solve: invalid problem"
    end if
    lim=huge(1.0_dp); if(present(tlimit)) lim=tlimit
    call gap_exact(cost,profit,budget,.true.,lim,r)
    r%unconstrained_max_profit=sum(maxval(profit,dim=1))
  end function gap_solve

  function aux_gap_bb(cost, profit_or_loss, budget, optim, tlimit, greedy_branching) result(r)
    real(dp),intent(in)::cost(:,:),profit_or_loss(:,:),budget(:)
    character(len=*),intent(in),optional::optim
    real(dp),intent(in),optional::tlimit
    logical,intent(in),optional::greedy_branching
    type(gap_result)::r
    logical::maximize
    real(dp)::lim
    maximize=.true.; if(present(optim)) maximize=trim(adjustl(optim))/='min'
    if (present(greedy_branching)) then
      if (greedy_branching .and. size(cost,2) < 0) error stop "aux_gap_bb: invalid problem"
    end if
    lim=huge(1.0_dp); if(present(tlimit)) lim=tlimit
    call gap_exact(cost,profit_or_loss,budget,maximize,lim,r)
  end function aux_gap_bb

  function aux_gap_bbdp(cost, profit_or_loss, budget, optim, tlimit, greedy_branching) result(r)
    real(dp),intent(in)::cost(:,:),profit_or_loss(:,:),budget(:)
    character(len=*),intent(in),optional::optim
    real(dp),intent(in),optional::tlimit
    logical,intent(in),optional::greedy_branching
    type(gap_result)::r
    r=aux_gap_bb(cost,profit_or_loss,budget,optim,tlimit,greedy_branching)
  end function aux_gap_bbdp

  function aux_gap_ga(cost, profit_or_loss, budget, trials, population_size, generations, seed, optim) result(r)
    real(dp),intent(in)::cost(:,:),profit_or_loss(:,:),budget(:)
    integer,intent(in)::trials,population_size,generations
    integer(i8),intent(in),optional::seed
    character(len=*),intent(in),optional::optim
    type(gap_result)::r
    integer::a,t,nagent,ntask,p,g,i,j,bestp
    integer,allocatable::pop(:,:),child(:),best(:)
    real(dp),allocatable::fit(:),acost(:)
    real(dp)::bestfit,u
    logical::maximize
    integer(i8)::state

    nagent=size(cost,1); ntask=size(cost,2)
    if(size(profit_or_loss,1)/=nagent .or. size(profit_or_loss,2)/=ntask .or. size(budget)/=nagent) &
      error stop "aux_gap_ga: dimension mismatch"
    maximize=.true.; if(present(optim)) maximize=trim(adjustl(optim))/='min'
    state=88172645463325252_i8; if(present(seed)) state=seed
    bestfit=merge(-huge(1.0_dp),huge(1.0_dp),maximize)
    allocate(best(ntask),acost(nagent)); best=1

    do t=1,max(1,trials)
      allocate(pop(max(2,population_size),ntask),fit(max(2,population_size)),child(ntask))
      do p=1,size(pop,1)
        do j=1,ntask
          pop(p,j)=1+int(rand_u(state)*real(nagent,dp)); pop(p,j)=min(nagent,pop(p,j))
        end do
        call repair_assignment(pop(p,:),cost,budget,profit_or_loss,maximize,state)
        fit(p)=assignment_value(pop(p,:),profit_or_loss)
      end do
      do g=1,max(1,generations)
        do p=1,size(pop,1)
          i=tournament(fit,maximize,state); j=tournament(fit,maximize,state)
          child=pop(i,:)
          do a=1,ntask
            if(rand_u(state)<0.5_dp) child(a)=pop(j,a)
            if(rand_u(state)<1.0_dp/real(max(1,ntask),dp)) then
              child(a)=1+int(rand_u(state)*real(nagent,dp)); child(a)=min(nagent,child(a))
            end if
          end do
          call repair_assignment(child,cost,budget,profit_or_loss,maximize,state)
          u=assignment_value(child,profit_or_loss)
          bestp=worst_index(fit,maximize)
          if(better(u,fit(bestp),maximize)) then; pop(bestp,:)=child; fit(bestp)=u; end if
        end do
      end do
      bestp=best_index(fit,maximize)
      if(better(fit(bestp),bestfit,maximize)) then; bestfit=fit(bestp); best=pop(bestp,:); end if
      deallocate(pop,fit,child)
    end do
    call assignment_costs(best,cost,acost)
    allocate(r%assignment(ntask),r%agent_cost(nagent)); r%assignment=best; r%agent_cost=acost
    r%total_profit_or_loss=bestfit; r%feasible=all(acost<=budget+1e-10_dp)
  end function aux_gap_ga

  subroutine gap_exact(cost,profit,budget,maximize,limit,r)
    real(dp),intent(in)::cost(:,:),profit(:,:),budget(:),limit
    logical,intent(in)::maximize
    type(gap_result),intent(out)::r
    integer::na,nt,j,k
    integer,allocatable::assign(:),bestassign(:),order(:)
    real(dp),allocatable::used(:),remain_bound(:)
    real(dp)::cur,best
    type(timer_type)::timer

    na=size(cost,1); nt=size(cost,2)
    if(size(profit,1)/=na .or. size(profit,2)/=nt .or. size(budget)/=na) error stop "gap_exact: dimension mismatch"
    allocate(assign(nt),bestassign(nt),used(na),order(nt),remain_bound(0:nt))
    assign=0; bestassign=0; used=0.0_dp
    do j=1,nt; order(j)=j; end do
    call sort_tasks_by_regret(order,profit,maximize)
    remain_bound(nt)=0.0_dp
    do k=nt-1,0,-1
      j=order(k+1)
      if(maximize) then
        remain_bound(k)=remain_bound(k+1)+maxval(profit(:,j))
      else
        remain_bound(k)=remain_bound(k+1)+minval(profit(:,j))
      end if
    end do
    best=merge(-huge(1.0_dp),huge(1.0_dp),maximize); cur=0.0_dp
    call timer%start(limit)
    call dfs(1)
    if(r%feasible) then
      allocate(r%assignment(nt),r%agent_cost(na)); r%assignment=bestassign
      call assignment_costs(bestassign,cost,r%agent_cost)
      r%total_profit_or_loss=best
    else
      allocate(r%assignment(0),r%agent_cost(na)); r%agent_cost=0.0_dp
    end if

  contains
    recursive subroutine dfs(pos)
      integer,intent(in)::pos
      integer::jj,aa,q
      integer,allocatable::ao(:)
      real(dp)::oldcur
      r%nodes=r%nodes+1_i8
      if(iand(r%nodes,1023_i8)==0_i8) then
        if(timer%expired()) then; r%timed_out=.true.; return; end if
      end if
      if(pos>nt) then
        if(.not.r%feasible .or. better(cur,best,maximize)) then
          best=cur; bestassign=assign; r%feasible=.true.
        end if
        return
      end if
      if(r%feasible) then
        if(maximize) then
          if(cur+remain_bound(pos-1)<=best) return
        else
          if(cur+remain_bound(pos-1)>=best) return
        end if
      end if
      jj=order(pos); allocate(ao(na)); ao=[(q,q=1,na)]
      call sort_agents(ao,profit(:,jj),maximize)
      oldcur=cur
      do q=1,na
        aa=ao(q); r%bkp_solved=r%bkp_solved+1_i8
        if(used(aa)+cost(aa,jj)<=budget(aa)+1e-12_dp) then
          used(aa)=used(aa)+cost(aa,jj); assign(jj)=aa; cur=oldcur+profit(aa,jj)
          call dfs(pos+1)
          used(aa)=used(aa)-cost(aa,jj); assign(jj)=0; cur=oldcur
          if(r%timed_out) return
        end if
      end do
    end subroutine dfs
  end subroutine gap_exact

  subroutine sort_tasks_by_regret(order,profit,maximize)
    integer,intent(inout)::order(:)
    real(dp),intent(in)::profit(:,:)
    logical,intent(in)::maximize
    integer::i,j,key
    real(dp)::rk,rj
    do i=2,size(order)
      key=order(i); rk=task_regret(profit(:,key),maximize); j=i-1
      do while(j>=1)
        rj=task_regret(profit(:,order(j)),maximize)
        if(rj>=rk) exit
        order(j+1)=order(j); j=j-1
      end do
      order(j+1)=key
    end do
  end subroutine sort_tasks_by_regret

  real(dp) function task_regret(x,maximize) result(rg)
    real(dp),intent(in)::x(:); logical,intent(in)::maximize
    real(dp)::b1,b2
    integer::i
    if(size(x)<2) then; rg=huge(1.0_dp); return; end if
    if(maximize) then
      b1=-huge(1.0_dp); b2=b1
      do i=1,size(x)
        if(x(i)>b1) then; b2=b1;b1=x(i); else if(x(i)>b2) then; b2=x(i); end if
      end do
      rg=b1-b2
    else
      b1=huge(1.0_dp); b2=b1
      do i=1,size(x)
        if(x(i)<b1) then; b2=b1;b1=x(i); else if(x(i)<b2) then; b2=x(i); end if
      end do
      rg=b2-b1
    end if
  end function task_regret

  subroutine sort_agents(ord,x,maximize)
    integer,intent(inout)::ord(:); real(dp),intent(in)::x(:); logical,intent(in)::maximize
    integer::i,j,key
    do i=2,size(ord)
      key=ord(i); j=i-1
      do while(j>=1)
        if(maximize) then
          if(x(ord(j))>=x(key)) exit
        else
          if(x(ord(j))<=x(key)) exit
        end if
        ord(j+1)=ord(j);j=j-1
      end do
      ord(j+1)=key
    end do
  end subroutine sort_agents

  logical pure function better(x,y,maximize)
    real(dp),intent(in)::x,y;logical,intent(in)::maximize
    if(maximize) then; better=x>y; else; better=x<y; end if
  end function better

  subroutine assignment_costs(a,cost,ac)
    integer,intent(in)::a(:); real(dp),intent(in)::cost(:,:); real(dp),intent(out)::ac(:)
    integer::j
    ac=0.0_dp
    do j=1,size(a); if(a(j)>0) ac(a(j))=ac(a(j))+cost(a(j),j); end do
  end subroutine assignment_costs

  real(dp) function assignment_value(a,p) result(v)
    integer,intent(in)::a(:); real(dp),intent(in)::p(:,:)
    integer::j
    v=0.0_dp; do j=1,size(a); v=v+p(a(j),j); end do
  end function assignment_value

  subroutine repair_assignment(a,cost,budget,p,maximize,state)
    integer,intent(inout)::a(:); real(dp),intent(in)::cost(:,:),budget(:),p(:,:)
    logical,intent(in)::maximize; integer(i8),intent(inout)::state
    real(dp),allocatable::ac(:)
    integer::bad,j,old,newa,tries
    allocate(ac(size(budget))); call assignment_costs(a,cost,ac)
    do tries=1,10000
      bad=0
      do j=1,size(ac); if(ac(j)>budget(j)+1e-12_dp) then; bad=j; exit; end if; end do
      if(bad==0) return
      j=1+int(rand_u(state)*real(size(a),dp)); j=min(size(a),j)
      if(a(j)/=bad) cycle
      old=a(j); newa=best_feasible_agent(j,a,cost,budget,ac,p,maximize)
      if(newa==0) newa=1+mod(old,size(budget))
      ac(old)=ac(old)-cost(old,j); a(j)=newa; ac(newa)=ac(newa)+cost(newa,j)
    end do
  end subroutine repair_assignment

  integer function best_feasible_agent(task,a,cost,budget,ac,p,maximize) result(b)
    integer,intent(in)::task,a(:);real(dp),intent(in)::cost(:,:),budget(:),ac(:),p(:,:);logical,intent(in)::maximize
    integer::q;real(dp)::bv
    b=0;bv=merge(-huge(1.0_dp),huge(1.0_dp),maximize)
    do q=1,size(budget)
      if(q==a(task)) cycle
      if(ac(q)+cost(q,task)<=budget(q)+1e-12_dp) then
        if(b==0 .or. better(p(q,task),bv,maximize)) then;b=q;bv=p(q,task);end if
      end if
    end do
  end function best_feasible_agent

  integer function tournament(fit,maximize,state) result(k)
    real(dp),intent(in)::fit(:);logical,intent(in)::maximize;integer(i8),intent(inout)::state
    integer::a,b
    a=1+int(rand_u(state)*size(fit));a=min(size(fit),a)
    b=1+int(rand_u(state)*size(fit));b=min(size(fit),b)
    if(better(fit(a),fit(b),maximize)) then;k=a;else;k=b;end if
  end function tournament

  integer function best_index(fit,maximize) result(k)
    real(dp),intent(in)::fit(:);logical,intent(in)::maximize;integer::i
    k=1;do i=2,size(fit);if(better(fit(i),fit(k),maximize))k=i;end do
  end function best_index
  integer function worst_index(fit,maximize) result(k)
    real(dp),intent(in)::fit(:);logical,intent(in)::maximize
    k=best_index(fit,.not.maximize)
  end function worst_index

  real(dp) function rand_u(state) result(u)
    integer(i8),intent(inout)::state
    integer(i8)::x
    x=state; x=ieor(x,shiftl(x,13)); x=ieor(x,shiftr(x,7)); x=ieor(x,shiftl(x,17)); state=x
    u=real(iand(x,int(z'7FFFFFFFFFFFFFFF',i8)),dp)/real(huge(1_i8),dp)
    if(u>=1.0_dp)u=nearest(1.0_dp,-1.0_dp)
  end function rand_u

end module flsss_gap
