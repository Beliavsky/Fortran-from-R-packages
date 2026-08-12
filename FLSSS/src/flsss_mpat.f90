module flsss_mpat
  use flsss_kinds, only : dp, i8
  use flsss_types, only : subset_solutions
  use flsss_util, only : timer_type, append_solution
  use flsss_fast_search, only : search_md_i8_packed
  use flsss_packed, only : packed_plan, make_packed_plan, pack_matrix
  use flsss_packed, only : packed_qualified, packed_completion_possible
  use flsss_packed, only : packed_meets_lower, packed_meets_upper
  use flsss_packed, only : packed_add_rows, is_comonotonic_i8
  implicit none
  private

  type :: triangular_packed_sums
    integer :: n = 0
    integer :: maxlen = 0
    integer :: nlane = 0
    integer :: nentry = 0
    integer, allocatable :: offset(:)
    integer(i8), allocatable :: value(:,:)
  contains
    procedure :: get => triangular_get
  end type triangular_packed_sums

  public :: search_md_i8_mpat

contains

  subroutine search_md_i8_mpat(v, len, target, me, result, solution_need, tlimit, lb, ub, dl, du)
    integer(i8), intent(in) :: v(:,:), target(:), me(:)
    integer, intent(in) :: len
    type(subset_solutions), intent(out) :: result
    integer, intent(in), optional :: solution_need
    real(dp), intent(in), optional :: tlimit
    integer, intent(in), optional :: lb(:), ub(:), dl, du

    type(packed_plan) :: plan
    type(triangular_packed_sums) :: tri
    type(timer_type) :: timer
    integer(i8), allocatable :: pv(:,:)
    integer, allocatable :: lower(:), upper(:)
    integer :: n, d, need, nlow, nup
    integer(i8), allocatable :: work_lb(:,:), work_ub(:,:)
    integer(i8), allocatable :: initial_sumlb(:), initial_sumub(:)

    n = size(v,1)
    d = size(v,2)
    need = 1
    if (present(solution_need)) need = max(1,solution_need)
    nlow = d
    if (present(dl)) nlow = max(0,min(d,dl))
    nup = d
    if (present(du)) nup = max(0,min(d,du))

    allocate(result%sol(0))
    result%engine = 'mpat-packed'
    if (len < 0 .or. len > n) return
    if (size(target) /= d .or. size(me) /= d) error stop 'search_md_i8_mpat: dimension mismatch'
    if (len == 0) then
      plan = make_packed_plan(v,len,target,me,nlow,nup)
      if (plan%valid .and. .not. plan%impossible) then
        allocate(initial_sumlb(plan%nlane))
        initial_sumlb=0_i8
        if(packed_qualified(plan,initial_sumlb)) call append_solution(result,[integer ::])
      end if
      return
    end if
    if (.not. is_comonotonic_i8(v)) then
      call fallback()
      return
    end if

    allocate(lower(len),upper(len))
    call initialize_bounds(n,len,lower,upper,lb,ub)
    plan = make_packed_plan(v,len,target,me,nlow,nup)
    if (.not. plan%valid) then
      call fallback()
      return
    end if
    result%packed_lanes = plan%nlane
    if (plan%impossible) then
      result%pruned = 1_i8
      return
    end if

    call pack_matrix(plan,v,pv)
    call build_triangular(pv,len,tri)
    result%tri_entries = int(tri%nentry,i8)

    allocate(work_lb(plan%nlane,len+1),work_ub(plan%nlane,len+1))
    if (present(tlimit)) then
      call timer%start(tlimit)
    else
      call timer%start(huge(1.0_dp))
    end if
    allocate(initial_sumlb(plan%nlane),initial_sumub(plan%nlane))
    call packed_add_rows(pv,lower,initial_sumlb)
    call packed_add_rows(pv,upper,initial_sumub)
    call visit(lower,upper,initial_sumlb,initial_sumub)

  contains

    recursive subroutine visit(lbin,ubin,slbin,subin)
      integer, intent(in) :: lbin(:),ubin(:)
      integer(i8), intent(in) :: slbin(:),subin(:)
      integer :: lbw(len),ubw(len),lleft(len),uleft(len),lright(len),uright(len)
      integer(i8) :: sumlb(plan%nlane),sumub(plan%nlane)
      integer(i8) :: slleft(plan%nlane),suleft(plan%nlane)
      integer(i8) :: slright(plan%nlane),suright(plan%nlane)
      integer :: ipos,cut
      logical :: ok,left_ok,right_ok

      if (result%size() >= need .or. result%timed_out) return
      result%nodes = result%nodes + 1_i8
      result%bound_states = result%bound_states + 1_i8
      if (iand(result%nodes,1023_i8) == 0_i8) then
        if (timer%expired()) then
          result%timed_out = .true.
          return
        end if
      end if

      lbw=lbin;ubw=ubin
      sumlb=slbin;sumub=subin
      call tighten_state(lbw,ubw,sumlb,sumub,ok)
      if(.not.ok)then
        result%pruned=result%pruned+1_i8
        return
      end if
      if(.not.packed_completion_possible(plan,0_i8*sumlb,sumlb,sumub))then
        result%pruned=result%pruned+1_i8
        return
      end if
      if(all(lbw==ubw))then
        if(packed_qualified(plan,sumlb))call append_solution(result,lbw)
        return
      end if

      ipos=smallest_gap(lbw,ubw)
      if(ipos==0)return
      cut=lbw(ipos)+(ubw(ipos)-lbw(ipos))/2
      result%mpat_splits=result%mpat_splits+1_i8

      lleft=lbw;uleft=ubw;slleft=sumlb;suleft=sumub
      call lower_upper(uleft,lleft,suleft,ipos,cut,left_ok)
      lright=lbw;uright=ubw;slright=sumlb;suright=sumub
      call raise_lower(lright,uright,slright,ipos,cut+1,right_ok)

      if(left_ok)call visit(lleft,uleft,slleft,suleft)
      if(result%size()>=need .or. result%timed_out)return
      if(right_ok)call visit(lright,uright,slright,suright)
    end subroutine visit

    subroutine tighten_state(lbx,ubx,sumlb,sumub,ok)
      integer, intent(inout) :: lbx(:),ubx(:)
      integer(i8), intent(inout) :: sumlb(:),sumub(:)
      logical, intent(out) :: ok
      integer :: pass, q, loq, hiq, midq, newq
      logical :: changed

      ok = .true.
      do pass = 1, 2*max(1,len)
        changed = .false.
        call build_bound_prefix(ubx,work_ub)
        do q = 1, len
          loq = lbx(q)
          hiq = ubx(q)
          if (loq >= hiq) cycle
          if (.not. lower_candidate_ok(q,loq,ubx,sumub)) then
            if (.not. lower_candidate_ok(q,hiq,ubx,sumub)) then
              ok = .false.
              return
            end if
            do while (loq < hiq)
              midq = loq + (hiq-loq)/2
              if (lower_candidate_ok(q,midq,ubx,sumub)) then
                hiq = midq
              else
                loq = midq + 1
              end if
            end do
            newq = loq
            if (newq > lbx(q)) then
              call raise_lower(lbx,ubx,sumlb,q,newq,ok)
              if (.not. ok) return
              result%bound_updates = result%bound_updates + 1_i8
              changed = .true.
            end if
          end if
        end do

        call build_bound_prefix(lbx,work_lb)
        do q = len, 1, -1
          loq = lbx(q)
          hiq = ubx(q)
          if (loq >= hiq) cycle
          if (.not. upper_candidate_ok(q,hiq,lbx,sumlb)) then
            if (.not. upper_candidate_ok(q,loq,lbx,sumlb)) then
              ok = .false.
              return
            end if
            do while (loq < hiq)
              midq = loq + (hiq-loq+1)/2
              if (upper_candidate_ok(q,midq,lbx,sumlb)) then
                loq = midq
              else
                hiq = midq - 1
              end if
            end do
            newq = loq
            if (newq < ubx(q)) then
              call lower_upper(ubx,lbx,sumub,q,newq,ok)
              if (.not. ok) return
              result%bound_updates = result%bound_updates + 1_i8
              changed = .true.
            end if
          end if
        end do
        if (.not. changed) exit
      end do
    end subroutine tighten_state

    subroutine build_bound_prefix(idx,pref)
      integer, intent(in) :: idx(:)
      integer(i8), intent(out) :: pref(:,:)
      integer :: i
      pref(:,1) = 0_i8
      do i = 1, len
        pref(:,i+1) = pref(:,i) + pv(idx(i),:)
      end do
    end subroutine build_bound_prefix

    logical function lower_candidate_ok(ipos,q,ubx,sumub) result(ok)
      integer, intent(in) :: ipos,q,ubx(:)
      integer(i8), intent(in) :: sumub(:)
      integer :: first, j, k, qstart
      integer(i8) :: candidate(plan%nlane), oldsum(plan%nlane), newsum(plan%nlane)

      first = ipos
      do j = ipos-1, 1, -1
        if (ubx(j) <= q-(ipos-j)) exit
        first = j
      end do
      k = ipos-first+1
      qstart = q-k+1
      oldsum = work_ub(:,ipos+1) - work_ub(:,first)
      result%tri_lookups=result%tri_lookups+1_i8
      call tri%get(k,qstart,newsum)
      candidate = sumub - oldsum + newsum
      ok = packed_meets_lower(plan,candidate)
    end function lower_candidate_ok

    logical function upper_candidate_ok(ipos,q,lbx,sumlb) result(ok)
      integer, intent(in) :: ipos,q,lbx(:)
      integer(i8), intent(in) :: sumlb(:)
      integer :: last, j, k
      integer(i8) :: candidate(plan%nlane), oldsum(plan%nlane), newsum(plan%nlane)

      last = ipos
      do j = ipos+1, len
        if (lbx(j) >= q+(j-ipos)) exit
        last = j
      end do
      k = last-ipos+1
      oldsum = work_lb(:,last+1) - work_lb(:,ipos)
      result%tri_lookups=result%tri_lookups+1_i8
      call tri%get(k,q,newsum)
      candidate = sumlb - oldsum + newsum
      ok = packed_meets_upper(plan,candidate)
    end function upper_candidate_ok

    subroutine raise_lower(lbx,ubx,sumlb,ipos,newq,ok)
      integer, intent(inout) :: lbx(:)
      integer, intent(in) :: ubx(:),ipos,newq
      integer(i8), intent(inout) :: sumlb(:)
      logical, intent(out) :: ok
      integer :: i, q, old
      ok = .true.
      q = newq
      do i = ipos, len
        if (i > ipos) q = max(lbx(i),q+1)
        if (q > ubx(i)) then
          ok = .false.
          return
        end if
        old = lbx(i)
        if (q <= old) exit
        sumlb = sumlb - pv(old,:) + pv(q,:)
        lbx(i) = q
      end do
    end subroutine raise_lower

    subroutine lower_upper(ubx,lbx,sumub,ipos,newq,ok)
      integer, intent(inout) :: ubx(:)
      integer, intent(in) :: lbx(:),ipos,newq
      integer(i8), intent(inout) :: sumub(:)
      logical, intent(out) :: ok
      integer :: i, q, old
      ok = .true.
      q = newq
      do i = ipos, 1, -1
        if (i < ipos) q = min(ubx(i),q-1)
        if (q < lbx(i)) then
          ok = .false.
          return
        end if
        old = ubx(i)
        if (q >= old) exit
        sumub = sumub - pv(old,:) + pv(q,:)
        ubx(i) = q
      end do
    end subroutine lower_upper

    integer function smallest_gap(lbx,ubx) result(ipos)
      integer, intent(in) :: lbx(:),ubx(:)
      integer :: i,g,best
      ipos = 0
      best = huge(1)
      do i = 1, len
        g = ubx(i)-lbx(i)
        if (g > 0 .and. g < best) then
          best = g
          ipos = i
        end if
      end do
    end function smallest_gap

    subroutine fallback()
      if (present(lb) .and. present(ub)) then
        call search_md_i8_packed(v,len,target,me,result,need,tlimit,lb,ub,nlow,nup)
      else if (present(lb)) then
        call search_md_i8_packed(v,len,target,me,result,need,tlimit,lb=lb,dl=nlow,du=nup)
      else if (present(ub)) then
        call search_md_i8_packed(v,len,target,me,result,need,tlimit,ub=ub,dl=nlow,du=nup)
      else
        call search_md_i8_packed(v,len,target,me,result,need,tlimit,dl=nlow,du=nup)
      end if
      if (trim(result%engine) == 'packed-dfs') result%engine = 'mpat-fallback'
    end subroutine fallback

  end subroutine search_md_i8_mpat

  subroutine build_triangular(pv,maxlen,tri)
    integer(i8), intent(in) :: pv(:,:)
    integer, intent(in) :: maxlen
    type(triangular_packed_sums), intent(out) :: tri
    integer :: n,k,s,e,total
    n = size(pv,1)
    tri%n=n
    tri%maxlen=maxlen
    tri%nlane=size(pv,2)
    allocate(tri%offset(maxlen+1))
    total=0
    do k=1,maxlen
      tri%offset(k)=total+1
      total=total+n-k+1
    end do
    tri%offset(maxlen+1)=total+1
    tri%nentry=total
    allocate(tri%value(tri%nlane,total))
    if(maxlen<=0)return
    do s=1,n
      e=tri%offset(1)+s-1
      tri%value(:,e)=pv(s,:)
    end do
    do k=2,maxlen
      do s=1,n-k+1
        e=tri%offset(k)+s-1
        tri%value(:,e)=tri%value(:,tri%offset(k-1)+s-1)+pv(s+k-1,:)
      end do
    end do
  end subroutine build_triangular

  subroutine triangular_get(self,k,start,s)
    class(triangular_packed_sums), intent(in) :: self
    integer, intent(in) :: k,start
    integer(i8), intent(out) :: s(:)
    integer :: e
    if(k<1 .or. k>self%maxlen) error stop 'triangular_get: invalid length'
    if(start<1 .or. start>self%n-k+1) error stop 'triangular_get: invalid start'
    if(size(s)/=self%nlane) error stop 'triangular_get: lane mismatch'
    e=self%offset(k)+start-1
    s=self%value(:,e)
  end subroutine triangular_get

  subroutine initialize_bounds(n,len,lower,upper,lb,ub)
    integer, intent(in) :: n,len
    integer, intent(out) :: lower(:),upper(:)
    integer, intent(in), optional :: lb(:),ub(:)
    integer :: i
    lower=[(i,i=1,len)]
    upper=[(n-len+i,i=1,len)]
    if(present(lb))then
      if(size(lb)/=len)error stop 'search_md_i8_mpat: lb has wrong size'
      lower=lb
    end if
    if(present(ub))then
      if(size(ub)/=len)error stop 'search_md_i8_mpat: ub has wrong size'
      upper=ub
    end if
    do i=1,len
      if(lower(i)<1 .or. upper(i)>n .or. lower(i)>upper(i)) &
        error stop 'search_md_i8_mpat: invalid bounds'
    end do
    do i=2,len
      if(lower(i)<=lower(i-1) .or. upper(i)<=upper(i-1)) &
        error stop 'search_md_i8_mpat: bounds must increase by position'
    end do
  end subroutine initialize_bounds


end module flsss_mpat
