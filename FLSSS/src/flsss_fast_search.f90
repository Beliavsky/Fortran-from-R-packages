module flsss_fast_search
  use flsss_kinds, only : dp, i8
  use flsss_types, only : subset_solutions
  use flsss_util, only : timer_type, append_solution
  use flsss_search, only : search_md_i8, build_envelope_md_i8
  use flsss_packed, only : packed_plan, make_packed_plan, pack_matrix, pack_vector
  use flsss_packed, only : packed_qualified, packed_completion_possible
  use flsss_packed, only : packed_meets_lower, packed_meets_upper
  use flsss_packed, only : packed_add_rows, is_comonotonic_i8
  implicit none
  private

  public :: search_md_i8_packed, search_md_i8_pat

  integer(i8), parameter :: iinf = 576460752303423487_i8

contains

  subroutine search_md_i8_packed(v, len, target, me, result, solution_need, tlimit, lb, ub, dl, du, prefix)
    integer(i8), intent(in) :: v(:,:), target(:), me(:)
    integer, intent(in) :: len
    type(subset_solutions), intent(out) :: result
    integer, intent(in), optional :: solution_need
    real(dp), intent(in), optional :: tlimit
    integer, intent(in), optional :: lb(:), ub(:), dl, du, prefix
    type(packed_plan) :: plan
    integer :: n, d, need, nlow, nup, pfx, pos, st
    integer, allocatable :: lower(:), upper(:), chosen(:)
    integer(i8), allocatable :: pv(:,:), sums(:)
    integer(i8), allocatable :: mn(:,:,:), mx(:,:,:), mnp(:,:,:), mxp(:,:,:)
    logical, allocatable :: state_ok(:,:)
    type(timer_type) :: timer

    n = size(v,1)
    d = size(v,2)
    need = 1
    if (present(solution_need)) need = max(1,solution_need)
    nlow = d
    if (present(dl)) nlow = max(0,min(d,dl))
    nup = d
    if (present(du)) nup = max(0,min(d,du))
    pfx = 0
    if (present(prefix)) pfx = prefix

    allocate(result%sol(0))
    result%engine = 'packed-dfs'
    if (len < 0 .or. len > n) return
    if (size(target) /= d .or. size(me) /= d) error stop "search_md_i8_packed: dimension mismatch"

    allocate(lower(len), upper(len), chosen(len))
    call initialize_bounds(n, len, lower, upper, lb, ub, "search_md_i8_packed")
    plan = make_packed_plan(v, len, target, me, nlow, nup)
    if (.not. plan%valid) then
      call fallback()
      return
    end if
    result%packed_lanes = plan%nlane
    if (plan%impossible) then
      result%pruned = 1_i8
      return
    end if

    call pack_matrix(plan, v, pv)
    allocate(sums(plan%nlane))
    sums = 0_i8
    call build_envelope_md_i8(v, lower, upper, mn, mx)
    allocate(mnp(len+1,n+2,plan%nlane), mxp(len+1,n+2,plan%nlane))
    allocate(state_ok(len+1,n+2))
    mnp = 0_i8
    mxp = 0_i8
    state_ok = .false.
    do pos = 1, len + 1
      do st = 1, n + 2
        if (pos == len + 1) then
          state_ok(pos,st) = .true.
        else if (mn(pos,st,1) < iinf) then
          state_ok(pos,st) = .true.
        end if
        if (state_ok(pos,st)) then
          call pack_vector(plan, mn(pos,st,:), mnp(pos,st,:))
          call pack_vector(plan, mx(pos,st,:), mxp(pos,st,:))
        end if
      end do
    end do
    deallocate(mn,mx)

    if (present(tlimit)) then
      call timer%start(tlimit)
    else
      call timer%start(huge(1.0_dp))
    end if

    if (len == 0) then
      if (packed_qualified(plan,sums)) call append_solution(result,[integer ::])
      return
    end if

    if (pfx > 0) then
      if (pfx < lower(1) .or. pfx > upper(1)) return
      chosen(1) = pfx
      sums = pv(pfx,:)
      if (.not. completion_possible(2,pfx+1)) then
        result%pruned = result%pruned + 1_i8
        return
      end if
      call dfs(2,pfx+1)
    else
      call dfs(1,1)
    end if

  contains

    recursive subroutine dfs(ipos,start)
      integer, intent(in) :: ipos, start
      integer :: q, lo, hi, rem
      integer(i8) :: old_local(plan%nlane)
      result%nodes = result%nodes + 1_i8
      if (iand(result%nodes,2047_i8) == 0_i8) then
        if (timer%expired()) then
          result%timed_out = .true.
          return
        end if
      end if
      if (result%size() >= need) return
      if (ipos > len) then
        if (packed_qualified(plan,sums)) call append_solution(result,chosen)
        return
      end if
      if (.not. completion_possible(ipos,start)) then
        result%pruned = result%pruned + 1_i8
        return
      end if
      rem = len - ipos
      lo = max(start,lower(ipos))
      hi = min(upper(ipos),n-rem)
      if (lo > hi) return
      old_local = sums
      do q = lo, hi
        chosen(ipos) = q
        sums = old_local + pv(q,:)
        if (.not. completion_possible(ipos+1,q+1)) then
          result%pruned = result%pruned + 1_i8
          cycle
        end if
        call dfs(ipos+1,q+1)
        if (result%timed_out .or. result%size() >= need) exit
      end do
      sums = old_local
    end subroutine dfs

    logical function completion_possible(ipos,start) result(ok)
      integer, intent(in) :: ipos,start
      integer :: s
      if (ipos > len) then
        ok = packed_qualified(plan,sums)
        return
      end if
      s = max(1,min(n+1,start))
      if (.not. state_ok(ipos,s)) then
        ok = .false.
      else
        ok = packed_completion_possible(plan,sums,mnp(ipos,s,:),mxp(ipos,s,:))
      end if
    end function completion_possible

    subroutine fallback()
      if (present(lb) .and. present(ub)) then
        call search_md_i8(v,len,target,me,result,need,tlimit,lb,ub,nlow,nup,pfx)
      else if (present(lb)) then
        call search_md_i8(v,len,target,me,result,need,tlimit,lb=lb,dl=nlow,du=nup,prefix=pfx)
      else if (present(ub)) then
        call search_md_i8(v,len,target,me,result,need,tlimit,ub=ub,dl=nlow,du=nup,prefix=pfx)
      else
        call search_md_i8(v,len,target,me,result,need,tlimit,dl=nlow,du=nup,prefix=pfx)
      end if
      result%engine = 'dfs-fallback'
    end subroutine fallback

  end subroutine search_md_i8_packed

  subroutine search_md_i8_pat(v, len, target, me, result, solution_need, tlimit, lb, ub, dl, du)
    integer(i8), intent(in) :: v(:,:), target(:), me(:)
    integer, intent(in) :: len
    type(subset_solutions), intent(out) :: result
    integer, intent(in), optional :: solution_need
    real(dp), intent(in), optional :: tlimit
    integer, intent(in), optional :: lb(:), ub(:), dl, du
    type(packed_plan) :: plan
    integer :: n, d, need, nlow, nup
    integer, allocatable :: lower(:), upper(:)
    integer(i8), allocatable :: pv(:,:)
    type(timer_type) :: timer

    n = size(v,1)
    d = size(v,2)
    need = 1
    if (present(solution_need)) need = max(1,solution_need)
    nlow = d
    if (present(dl)) nlow = max(0,min(d,dl))
    nup = d
    if (present(du)) nup = max(0,min(d,du))
    allocate(result%sol(0))
    result%engine = 'pat-packed'
    if (len < 0 .or. len > n) return
    if (size(target) /= d .or. size(me) /= d) error stop "search_md_i8_pat: dimension mismatch"
    if (.not. is_comonotonic_i8(v)) then
      call packed_fallback()
      return
    end if
    allocate(lower(len),upper(len))
    call initialize_bounds(n,len,lower,upper,lb,ub,"search_md_i8_pat")
    plan = make_packed_plan(v,len,target,me,nlow,nup)
    if (.not. plan%valid) then
      call packed_fallback()
      return
    end if
    result%packed_lanes = plan%nlane
    if (plan%impossible) then
      result%pruned = 1_i8
      return
    end if
    call pack_matrix(plan,v,pv)
    if (present(tlimit)) then
      call timer%start(tlimit)
    else
      call timer%start(huge(1.0_dp))
    end if
    call visit(lower,upper)

  contains

    recursive subroutine visit(lbi,ubi)
      integer, intent(in) :: lbi(:), ubi(:)
      integer :: lbw(len), ubw(len)
      integer :: lleft(len), uleft(len), lright(len), uright(len)
      integer(i8) :: smin(plan%nlane), smax(plan%nlane)
      integer :: p, i, gap, bestgap, mid
      logical :: feasible

      result%nodes = result%nodes + 1_i8
      result%bound_states = result%bound_states + 1_i8
      if (iand(result%nodes,1023_i8) == 0_i8) then
        if (timer%expired()) then
          result%timed_out = .true.
          return
        end if
      end if
      if (result%size() >= need) return
      lbw=lbi;ubw=ubi
      call tighten_state(lbw,ubw,feasible)
      if (.not.feasible) then
        result%pruned = result%pruned + 1_i8
        return
      end if
      call packed_add_rows(pv,lbw,smin)
      call packed_add_rows(pv,ubw,smax)
      if (.not. packed_completion_possible(plan,0_i8*smin,smin,smax)) then
        result%pruned = result%pruned + 1_i8
        return
      end if
      if (all(lbw == ubw)) then
        if (packed_qualified(plan,smin)) call append_solution(result,lbw)
        return
      end if

      p = 0
      bestgap = huge(1)
      do i = 1, len
        gap = ubw(i) - lbw(i)
        if (gap > 0 .and. gap < bestgap) then
          bestgap = gap
          p = i
        end if
      end do
      if (p == 0) return
      mid = lbw(p) + (ubw(p)-lbw(p))/2

      lleft = lbw
      uleft = ubw
      uleft(p) = mid
      do i = p-1, 1, -1
        uleft(i) = min(uleft(i),uleft(i+1)-1)
      end do
      if (all(lleft <= uleft)) call visit(lleft,uleft)
      if (result%timed_out .or. result%size() >= need) return

      lright = lbw
      uright = ubw
      lright(p) = mid + 1
      do i = p+1, len
        lright(i) = max(lright(i),lright(i-1)+1)
      end do
      if (all(lright <= uright)) call visit(lright,uright)
    end subroutine visit

    subroutine tighten_state(lbx,ubx,feasible)
      integer,intent(inout)::lbx(:),ubx(:)
      logical,intent(out)::feasible
      integer::p,i,loq,hiq,midq,newq,pass
      logical::changed
      feasible=.true.
      do pass=1,2*max(1,len)
        changed=.false.
        do p=1,len
          loq=lbx(p);hiq=ubx(p)
          if(loq<hiq) then
            if(lower_candidate_ok(p,loq,lbx,ubx)) cycle
            if(.not.lower_candidate_ok(p,hiq,lbx,ubx)) then
              feasible=.false.;return
            end if
            do while(loq<hiq)
              midq=loq+(hiq-loq)/2
              if(lower_candidate_ok(p,midq,lbx,ubx)) then
                hiq=midq
              else
                loq=midq+1
              end if
            end do
            newq=loq
            if(newq>lbx(p)) then
              lbx(p)=newq;changed=.true.
              do i=p+1,len
                lbx(i)=max(lbx(i),lbx(i-1)+1)
              end do
              if(any(lbx>ubx)) then;feasible=.false.;return;end if
            end if
          end if
        end do
        do p=len,1,-1
          loq=lbx(p);hiq=ubx(p)
          if(loq<hiq) then
            if(upper_candidate_ok(p,hiq,lbx,ubx)) cycle
            if(.not.upper_candidate_ok(p,loq,lbx,ubx)) then
              feasible=.false.;return
            end if
            do while(loq<hiq)
              midq=loq+(hiq-loq+1)/2
              if(upper_candidate_ok(p,midq,lbx,ubx)) then
                loq=midq
              else
                hiq=midq-1
              end if
            end do
            newq=loq
            if(newq<ubx(p)) then
              ubx(p)=newq;changed=.true.
              do i=p-1,1,-1
                ubx(i)=min(ubx(i),ubx(i+1)-1)
              end do
              if(any(lbx>ubx)) then;feasible=.false.;return;end if
            end if
          end if
        end do
        if(.not.changed)exit
      end do
    end subroutine tighten_state

    logical function lower_candidate_ok(p,q,lbx,ubx) result(ok)
      integer,intent(in)::p,q,lbx(:),ubx(:)
      integer::ix(len),i
      integer(i8)::s(plan%nlane)
      ix=ubx
      ix(p)=q
      do i=p-1,1,-1
        ix(i)=min(ix(i),ix(i+1)-1)
      end do
      if(any(ix<lbx)) then;ok=.false.;return;end if
      call packed_add_rows(pv,ix,s)
      ok=packed_meets_lower(plan,s)
    end function lower_candidate_ok

    logical function upper_candidate_ok(p,q,lbx,ubx) result(ok)
      integer,intent(in)::p,q,lbx(:),ubx(:)
      integer::ix(len),i
      integer(i8)::s(plan%nlane)
      ix=lbx
      ix(p)=q
      do i=p+1,len
        ix(i)=max(ix(i),ix(i-1)+1)
      end do
      if(any(ix>ubx)) then;ok=.false.;return;end if
      call packed_add_rows(pv,ix,s)
      ok=packed_meets_upper(plan,s)
    end function upper_candidate_ok

    subroutine packed_fallback()
      if (present(lb) .and. present(ub)) then
        call search_md_i8_packed(v,len,target,me,result,need,tlimit,lb,ub,nlow,nup)
      else if (present(lb)) then
        call search_md_i8_packed(v,len,target,me,result,need,tlimit,lb=lb,dl=nlow,du=nup)
      else if (present(ub)) then
        call search_md_i8_packed(v,len,target,me,result,need,tlimit,ub=ub,dl=nlow,du=nup)
      else
        call search_md_i8_packed(v,len,target,me,result,need,tlimit,dl=nlow,du=nup)
      end if
      if (trim(result%engine) == 'packed-dfs') result%engine = 'packed-fallback'
    end subroutine packed_fallback

  end subroutine search_md_i8_pat

  subroutine initialize_bounds(n,len,lower,upper,lb,ub,where)
    integer, intent(in) :: n,len
    integer, intent(out) :: lower(:),upper(:)
    integer, intent(in), optional :: lb(:),ub(:)
    character(len=*), intent(in) :: where
    integer :: i
    if (len == 0) return
    lower = [(i,i=1,len)]
    upper = [(n-len+i,i=1,len)]
    if (present(lb)) then
      if (size(lb) /= len) error stop trim(where)//": lb has wrong size"
      lower = lb
    end if
    if (present(ub)) then
      if (size(ub) /= len) error stop trim(where)//": ub has wrong size"
      upper = ub
    end if
    do i = 1, len
      if (lower(i) < 1 .or. upper(i) > n .or. lower(i) > upper(i)) then
        error stop trim(where)//": invalid bounds"
      end if
    end do
    do i = 2, len
      if (lower(i) <= lower(i-1) .or. upper(i) <= upper(i-1)) then
        error stop trim(where)//": bounds must increase by position"
      end if
    end do
  end subroutine initialize_bounds

end module flsss_fast_search
