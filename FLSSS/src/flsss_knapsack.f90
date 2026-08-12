module flsss_knapsack
  use flsss_kinds, only : dp, i8
  use flsss_types, only : knapsack_result, knapsack_multi_result
  use flsss_util, only : timer_type, argsort_real, top_k_sum
  use flsss_integerize, only : integerize_problem
  implicit none
  private
  public :: mm_knapsack, mm_knapsack_integerized, aux_knapsack01bb, aux_knapsack01dp

contains

  function mm_knapsack(len, items_profits, items_costs, capacities, heuristic, tlimit) result(r)
    integer, intent(in) :: len
    real(dp), intent(in) :: items_profits(:), items_costs(:,:), capacities(:)
    logical, intent(in), optional :: heuristic
    real(dp), intent(in), optional :: tlimit
    type(knapsack_result) :: r
    real(dp) :: lim
    if (present(heuristic)) then
      if (heuristic .and. len < 0) error stop "mm_knapsack: invalid len"
    end if
    lim = huge(1.0_dp); if (present(tlimit)) lim = tlimit
    call solve_md_knapsack(len, items_profits, items_costs, capacities, lim, r)
  end function mm_knapsack

  function mm_knapsack_integerized(len, items_profits, items_costs, capacities, precision_level, &
                                    return_before_mining, tlimit) result(r)
    integer, intent(in) :: len
    real(dp), intent(in) :: items_profits(:), items_costs(:,:), capacities(:)
    integer, intent(in), optional :: precision_level(:)
    logical, intent(in), optional :: return_before_mining
    real(dp), intent(in), optional :: tlimit
    type(knapsack_result) :: r
    real(dp), allocatable :: target(:), me(:)
    logical :: rbm
    integer :: k
    k = len
    if (k == 0) k = size(items_profits)
    allocate(target(size(capacities)), me(size(capacities)))
    target = capacities / 2.0_dp; me = capacities / 2.0_dp
    r%integerized = integerize_problem(k, items_costs, target, max(me, tiny(1.0_dp)), precision_level)
    rbm = .false.; if (present(return_before_mining)) rbm = return_before_mining
    if (rbm) then
      allocate(r%solution(0), r%selection_costs(size(capacities)), r%budgets(size(capacities)))
      r%selection_costs = 0.0_dp; r%budgets = capacities
      return
    end if
    r = mm_knapsack(len, items_profits, items_costs, capacities, tlimit=tlimit)
    r%integerized = integerize_problem(k, items_costs, target, max(me, tiny(1.0_dp)), precision_level)
  end function mm_knapsack_integerized

  subroutine solve_md_knapsack(len, profits, costs, caps, limit, r)
    integer, intent(in) :: len
    real(dp), intent(in) :: profits(:), costs(:,:), caps(:), limit
    type(knapsack_result), intent(out) :: r
    integer :: n,d,i,k, fixed_k
    integer, allocatable :: ord(:), sel(:), bestsel(:)
    real(dp), allocatable :: csum(:)
    real(dp) :: best, pcur, ub
    type(timer_type) :: timer

    n=size(profits); d=size(costs,2)
    if (size(costs,1)/=n .or. size(caps)/=d) error stop "mm_knapsack: dimension mismatch"
    ord=argsort_real(profits,.true.)
    allocate(sel(n),bestsel(n),csum(d)); sel=0; bestsel=0; csum=0.0_dp
    best=-huge(1.0_dp); pcur=0.0_dp
    fixed_k=len
    call timer%start(limit)
    call dfs(1,0)
    allocate(r%budgets(d),r%selection_costs(d)); r%budgets=caps
    if (r%feasible) then
      allocate(r%solution(count(bestsel==1)))
      k=0
      do i=1,n
        if(bestsel(i)==1) then; k=k+1; r%solution(k)=i; end if
      end do
      r%selection_costs=0.0_dp
      if(size(r%solution)>0) r%selection_costs=sum(costs(r%solution,:),dim=1)
      r%selection_profit=best
    else
      allocate(r%solution(0)); r%selection_costs=0.0_dp
    end if
    if (len>0) then
      r%unconstrained_max_profit=top_k_sum(profits,len)
    else
      r%unconstrained_max_profit=sum(max(profits,0.0_dp))
    end if

  contains
    recursive subroutine dfs(pos,nsel)
      integer,intent(in)::pos,nsel
      integer::item,need
      real(dp)::savep
      real(dp),allocatable::savec(:)
      r%nodes=r%nodes+1_i8
      if(iand(r%nodes,1023_i8)==0_i8) then
        if(timer%expired()) then; r%timed_out=.true.; return; end if
      end if
      if(pos>n) then
        if((fixed_k==0 .or. nsel==fixed_k) .and. pcur>best) then
          best=pcur; bestsel=sel; r%feasible=.true.
        end if
        return
      end if
      if(fixed_k>0) then
        need=fixed_k-nsel
        if(need<0 .or. need>n-pos+1) return
        if(need==0) then
          if(pcur>best) then; best=pcur; bestsel=sel; r%feasible=.true.; end if
          return
        end if
        ub=pcur+sum(profits(ord(pos:min(n,pos+need-1))))
      else
        ub=pcur+sum(max(profits(ord(pos:n)),0.0_dp))
      end if
      if(ub<=best) return
      item=ord(pos)
      allocate(savec(d)); savec=csum; savep=pcur
      if(all(csum+costs(item,:)<=caps+1.0e-12_dp)) then
        sel(item)=1; csum=csum+costs(item,:); pcur=pcur+profits(item)
        call dfs(pos+1,nsel+1)
        if(r%timed_out) return
        sel(item)=0; csum=savec; pcur=savep
      end if
      call dfs(pos+1,nsel)
    end subroutine dfs
  end subroutine solve_md_knapsack

  function aux_knapsack01bb(weight,value,caps,item_ncaps,tlimit) result(r)
    real(dp),intent(in)::weight(:),value(:),caps(:)
    integer,intent(in),optional::item_ncaps(:)
    real(dp),intent(in),optional::tlimit
    type(knapsack_multi_result)::r
    integer::j,kmax
    type(knapsack_result)::one
    real(dp),allocatable::c(:,:)
    real(dp)::lim
    lim=huge(1.0_dp); if(present(tlimit)) lim=tlimit
    allocate(r%max_value(size(caps)),r%selection(size(caps)),c(size(weight),1)); c(:,1)=weight
    do j=1,size(caps)
      kmax=0
      if(present(item_ncaps)) then
        if(size(item_ncaps)>0) kmax=item_ncaps(j)
      end if
      if(kmax>0) then
        call solve_md_knapsack_at_most(kmax,value,c,[caps(j)],lim,one)
      else
        call solve_md_knapsack(0,value,c,[caps(j)],lim,one)
      end if
      r%max_value(j)=one%selection_profit
      allocate(r%selection(j)%idx(size(one%solution))); r%selection(j)%idx=one%solution
      r%nodes=r%nodes+one%nodes; r%timed_out=r%timed_out .or. one%timed_out
    end do
  end function aux_knapsack01bb

  subroutine solve_md_knapsack_at_most(kmax,profits,costs,caps,limit,r)
    integer,intent(in)::kmax
    real(dp),intent(in)::profits(:),costs(:,:),caps(:),limit
    type(knapsack_result),intent(out)::r
    integer::k
    type(knapsack_result)::tmp
    r%selection_profit=-huge(1.0_dp)
    do k=0,kmax
      call solve_md_knapsack(k,profits,costs,caps,limit,tmp)
      if(tmp%feasible .and. tmp%selection_profit>r%selection_profit) r=tmp
    end do
  end subroutine solve_md_knapsack_at_most

  function aux_knapsack01dp(weight,value,caps,tlimit) result(r)
    integer,intent(in)::weight(:),caps(:)
    real(dp),intent(in)::value(:)
    real(dp),intent(in),optional::tlimit
    type(knapsack_multi_result)::r
    integer::n,maxcap,i,c,j,w
    real(dp),allocatable::dpv(:,:)
    logical,allocatable::take(:,:)
    type(timer_type)::timer
    real(dp)::lim
    n=size(weight); maxcap=maxval(caps)
    if(size(value)/=n) error stop "aux_knapsack01dp: size mismatch"
    allocate(dpv(0:n,0:maxcap),take(n,0:maxcap)); dpv=0.0_dp; take=.false.
    lim=huge(1.0_dp); if(present(tlimit)) lim=tlimit; call timer%start(lim)
    do i=1,n
      dpv(i,:)=dpv(i-1,:)
      w=weight(i)
      do c=max(0,w),maxcap
        if(dpv(i-1,c-w)+value(i)>dpv(i,c)) then
          dpv(i,c)=dpv(i-1,c-w)+value(i); take(i,c)=.true.
        end if
      end do
      if(timer%expired()) then; r%timed_out=.true.; exit; end if
    end do
    allocate(r%max_value(size(caps)),r%selection(size(caps)),r%lookup_table(0:n,0:maxcap))
    r%lookup_table=dpv
    do j=1,size(caps)
      r%max_value(j)=dpv(n,caps(j)); c=caps(j)
      allocate(r%selection(j)%idx(0))
      do i=n,1,-1
        if (take(i,c)) then
          if (c >= weight(i)) then
            if (abs(dpv(i,c) - (dpv(i-1,c-weight(i)) + value(i))) < 1.0e-10_dp) then
              call append_index(r%selection(j)%idx,i)
              c=c-weight(i)
            end if
          end if
        end if
      end do
      if(size(r%selection(j)%idx)>1) r%selection(j)%idx=r%selection(j)%idx(size(r%selection(j)%idx):1:-1)
    end do
  end function aux_knapsack01dp

  subroutine append_index(a,x)
    integer,allocatable,intent(inout)::a(:)
    integer,intent(in)::x
    integer,allocatable::b(:)
    allocate(b(size(a)+1)); if(size(a)>0)b(:size(a))=a; b(size(b))=x; call move_alloc(b,a)
  end subroutine append_index

end module flsss_knapsack
