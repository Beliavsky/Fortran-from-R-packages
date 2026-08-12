module flsss_search
  use flsss_kinds, only : dp, i8
  use flsss_types, only : subset_solutions
  use flsss_util, only : timer_type, append_solution
  implicit none
  private
  public :: search_1d, search_md, search_md_i8, build_envelope_md_i8

  real(dp), parameter :: rinf = huge(1.0_dp) / 16.0_dp
  integer(i8), parameter :: iinf = 576460752303423487_i8

contains

  subroutine search_1d(v, len, target, me, result, solution_need, tlimit, lb, ub, prefix)
    real(dp), intent(in) :: v(:), target, me
    integer, intent(in) :: len
    type(subset_solutions), intent(out) :: result
    integer, intent(in), optional :: solution_need
    real(dp), intent(in), optional :: tlimit
    integer, intent(in), optional :: lb(:), ub(:)
    integer, intent(in), optional :: prefix
    integer :: need, n, pfx, i
    integer, allocatable :: lower(:), upper(:), chosen(:)
    real(dp), allocatable :: minrem(:,:), maxrem(:,:)
    type(timer_type) :: timer

    n = size(v)
    need = 1
    if (present(solution_need)) need = max(1, solution_need)
    pfx = 0
    if (present(prefix)) pfx = prefix
    allocate(result%sol(0))
    if (len < 0 .or. len > n) return

    allocate(lower(len), upper(len), chosen(len))
    if (len > 0) then
      lower = [(i, i=1,len)]
      upper = [(n - len + i, i=1,len)]
      if (present(lb)) then
        if (size(lb) /= len) error stop "search_1d: lb has wrong size"
        lower = lb
      end if
      if (present(ub)) then
        if (size(ub) /= len) error stop "search_1d: ub has wrong size"
        upper = ub
      end if
      call validate_bounds(lower, upper, n, "search_1d")
    end if
    if (present(tlimit)) then
      call timer%start(tlimit)
    else
      call timer%start(huge(1.0_dp))
    end if

    if (len == 0) then
      if (abs(target) <= me) call append_solution(result, [integer ::])
      return
    end if

    ! v0.2: exact completion envelopes for every DFS state.  Unlike the
    ! v0.1 monotone-only bound, this is valid for arbitrary value order and
    ! position-specific lower/upper index bounds.
    call build_envelope_1d(v, lower, upper, minrem, maxrem)

    if (pfx > 0) then
      if (pfx < lower(1) .or. pfx > upper(1)) return
      if (.not. completion_possible_1d(2, pfx + 1, v(pfx))) then
        result%pruned = result%pruned + 1_i8
        return
      end if
      chosen(1) = pfx
      call dfs(2, pfx + 1, v(pfx))
    else
      call dfs(1, 1, 0.0_dp)
    end if

  contains

    recursive subroutine dfs(pos, start, s)
      integer, intent(in) :: pos, start
      real(dp), intent(in) :: s
      integer :: q, lo, hi, rem
      result%nodes = result%nodes + 1_i8
      if (iand(result%nodes, 1023_i8) == 0_i8) then
        if (timer%expired()) then
          result%timed_out = .true.
          return
        end if
      end if
      if (size(result%sol) >= need) return
      if (pos > len) then
        if (abs(s - target) <= me) call append_solution(result, chosen)
        return
      end if
      if (.not. completion_possible_1d(pos, start, s)) then
        result%pruned = result%pruned + 1_i8
        return
      end if
      rem = len - pos
      lo = max(start, lower(pos))
      hi = min(upper(pos), n - rem)
      if (lo > hi) return
      do q = lo, hi
        if (.not. completion_possible_1d(pos + 1, q + 1, s + v(q))) then
          result%pruned = result%pruned + 1_i8
          cycle
        end if
        chosen(pos) = q
        call dfs(pos + 1, q + 1, s + v(q))
        if (result%timed_out .or. size(result%sol) >= need) return
      end do
    end subroutine dfs

    logical function completion_possible_1d(pos, start, s) result(ok)
      integer, intent(in) :: pos, start
      real(dp), intent(in) :: s
      integer :: st
      if (pos > len) then
        ok = abs(s - target) <= me
        return
      end if
      st = max(1, min(n + 1, start))
      if (minrem(pos,st) >= 0.5_dp*rinf) then
        ok = .false.
      else
        ok = (s + minrem(pos,st) <= target + me) .and. &
             (s + maxrem(pos,st) >= target - me)
      end if
    end function completion_possible_1d

  end subroutine search_1d

  subroutine search_md(v, len, target, me, result, solution_need, tlimit, lb, ub, dl, du, prefix)
    real(dp), intent(in) :: v(:,:), target(:), me(:)
    integer, intent(in) :: len
    type(subset_solutions), intent(out) :: result
    integer, intent(in), optional :: solution_need
    real(dp), intent(in), optional :: tlimit
    integer, intent(in), optional :: lb(:), ub(:), dl, du, prefix
    integer :: need, n, d, nlow, nup, pfx, i
    integer, allocatable :: lower(:), upper(:), chosen(:)
    real(dp), allocatable :: sums(:), minrem(:,:,:), maxrem(:,:,:)
    type(timer_type) :: timer

    n = size(v,1)
    d = size(v,2)
    allocate(result%sol(0))
    if (size(target) /= d .or. size(me) /= d) error stop "search_md: dimension mismatch"
    if (len < 0 .or. len > n) return
    need = 1
    if (present(solution_need)) need = max(1, solution_need)
    nlow = d
    if (present(dl)) nlow = max(0, min(d, dl))
    nup = d
    if (present(du)) nup = max(0, min(d, du))
    pfx = 0
    if (present(prefix)) pfx = prefix
    allocate(lower(len), upper(len), chosen(len), sums(d))
    if (len > 0) then
      lower = [(i, i=1,len)]
      upper = [(n-len+i, i=1,len)]
      if (present(lb)) then
        if (size(lb) /= len) error stop "search_md: lb has wrong size"
        lower = lb
      end if
      if (present(ub)) then
        if (size(ub) /= len) error stop "search_md: ub has wrong size"
        upper = ub
      end if
      call validate_bounds(lower, upper, n, "search_md")
    end if
    sums = 0.0_dp
    if (present(tlimit)) then
      call timer%start(tlimit)
    else
      call timer%start(huge(1.0_dp))
    end if

    if (len == 0) then
      if (qualified_real(sums, target, me, nlow, nup)) call append_solution(result, [integer ::])
      return
    end if

    call build_envelope_md_real(v, lower, upper, minrem, maxrem)

    if (pfx > 0) then
      if (pfx < lower(1) .or. pfx > upper(1)) return
      chosen(1) = pfx
      sums = v(pfx,:)
      if (.not. completion_possible(2, pfx + 1)) then
        result%pruned = result%pruned + 1_i8
        return
      end if
      call dfs(2, pfx + 1)
    else
      call dfs(1, 1)
    end if

  contains

    recursive subroutine dfs(pos, start)
      integer, intent(in) :: pos, start
      integer :: q, lo, hi, rem
      real(dp), allocatable :: old(:)
      result%nodes = result%nodes + 1_i8
      if (iand(result%nodes, 1023_i8) == 0_i8) then
        if (timer%expired()) then
          result%timed_out = .true.
          return
        end if
      end if
      if (size(result%sol) >= need) return
      if (pos > len) then
        if (qualified_real(sums, target, me, nlow, nup)) call append_solution(result, chosen)
        return
      end if
      if (.not. completion_possible(pos, start)) then
        result%pruned = result%pruned + 1_i8
        return
      end if
      rem = len - pos
      lo = max(start, lower(pos))
      hi = min(upper(pos), n-rem)
      if (lo > hi) return
      allocate(old(d))
      old = sums
      do q = lo, hi
        chosen(pos) = q
        sums = old + v(q,:)
        if (.not. completion_possible(pos + 1, q + 1)) then
          result%pruned = result%pruned + 1_i8
          cycle
        end if
        call dfs(pos + 1, q + 1)
        if (result%timed_out .or. size(result%sol) >= need) exit
      end do
      sums = old
    end subroutine dfs

    logical function completion_possible(pos, start) result(ok)
      integer, intent(in) :: pos, start
      integer :: j, st
      if (pos > len) then
        ok = qualified_real(sums, target, me, nlow, nup)
        return
      end if
      st = max(1, min(n + 1, start))
      if (minrem(pos,st,1) >= 0.5_dp*rinf) then
        ok = .false.
        return
      end if
      ok = .true.
      do j = 1, d
        if (j <= nlow) then
          if (sums(j) + maxrem(pos,st,j) < target(j) - me(j)) then
            ok = .false.
            return
          end if
        end if
        if (j > d - nup) then
          if (sums(j) + minrem(pos,st,j) > target(j) + me(j)) then
            ok = .false.
            return
          end if
        end if
      end do
    end function completion_possible

  end subroutine search_md

  subroutine search_md_i8(v, len, target, me, result, solution_need, tlimit, lb, ub, dl, du, prefix)
    integer(i8), intent(in) :: v(:,:), target(:), me(:)
    integer, intent(in) :: len
    type(subset_solutions), intent(out) :: result
    integer, intent(in), optional :: solution_need
    real(dp), intent(in), optional :: tlimit
    integer, intent(in), optional :: lb(:), ub(:), dl, du, prefix
    integer :: need, n, d, nlow, nup, pfx, i
    integer, allocatable :: lower(:), upper(:), chosen(:)
    integer(i8), allocatable :: sums(:), minrem(:,:,:), maxrem(:,:,:)
    type(timer_type) :: timer

    n = size(v,1)
    d = size(v,2)
    allocate(result%sol(0))
    if (size(target) /= d .or. size(me) /= d) error stop "search_md_i8: dimension mismatch"
    if (len < 0 .or. len > n) return
    need = 1
    if (present(solution_need)) need = max(1, solution_need)
    nlow = d
    if (present(dl)) nlow = max(0, min(d, dl))
    nup = d
    if (present(du)) nup = max(0, min(d, du))
    pfx = 0
    if (present(prefix)) pfx = prefix
    allocate(lower(len), upper(len), chosen(len), sums(d))
    if (len > 0) then
      lower = [(i, i=1,len)]
      upper = [(n-len+i, i=1,len)]
      if (present(lb)) then
        if (size(lb) /= len) error stop "search_md_i8: lb has wrong size"
        lower = lb
      end if
      if (present(ub)) then
        if (size(ub) /= len) error stop "search_md_i8: ub has wrong size"
        upper = ub
      end if
      call validate_bounds(lower, upper, n, "search_md_i8")
    end if
    sums = 0_i8
    if (present(tlimit)) then
      call timer%start(tlimit)
    else
      call timer%start(huge(1.0_dp))
    end if

    if (len == 0) then
      if (qualified_i8(sums, target, me, nlow, nup)) call append_solution(result, [integer ::])
      return
    end if

    call build_envelope_md_i8(v, lower, upper, minrem, maxrem)

    if (pfx > 0) then
      if (pfx < lower(1) .or. pfx > upper(1)) return
      chosen(1) = pfx
      sums = v(pfx,:)
      if (.not. completion_possible(2, pfx + 1)) then
        result%pruned = result%pruned + 1_i8
        return
      end if
      call dfs(2, pfx + 1)
    else
      call dfs(1, 1)
    end if

  contains

    recursive subroutine dfs(pos, start)
      integer, intent(in) :: pos, start
      integer :: q, lo, hi, rem, j
      integer(i8), allocatable :: old(:)
      result%nodes = result%nodes + 1_i8
      if (iand(result%nodes, 1023_i8) == 0_i8) then
        if (timer%expired()) then
          result%timed_out = .true.
          return
        end if
      end if
      if (size(result%sol) >= need) return
      if (pos > len) then
        if (qualified_i8(sums, target, me, nlow, nup)) call append_solution(result, chosen)
        return
      end if
      if (.not. completion_possible(pos, start)) then
        result%pruned = result%pruned + 1_i8
        return
      end if
      rem = len - pos
      lo = max(start, lower(pos))
      hi = min(upper(pos), n-rem)
      if (lo > hi) return
      allocate(old(d))
      old = sums
      do q = lo, hi
        chosen(pos) = q
        do j = 1, d
          sums(j) = sat_add_i8(old(j), v(q,j))
        end do
        if (.not. completion_possible(pos + 1, q + 1)) then
          result%pruned = result%pruned + 1_i8
          cycle
        end if
        call dfs(pos + 1, q + 1)
        if (result%timed_out .or. size(result%sol) >= need) exit
      end do
      sums = old
    end subroutine dfs

    logical function completion_possible(pos, start) result(ok)
      integer, intent(in) :: pos, start
      integer :: j, st
      integer(i8) :: lo_target, hi_target, lo_possible, hi_possible
      if (pos > len) then
        ok = qualified_i8(sums, target, me, nlow, nup)
        return
      end if
      st = max(1, min(n + 1, start))
      if (minrem(pos,st,1) >= iinf) then
        ok = .false.
        return
      end if
      ok = .true.
      do j = 1, d
        lo_target = sat_add_i8(target(j), -me(j))
        hi_target = sat_add_i8(target(j), me(j))
        lo_possible = sat_add_i8(sums(j), minrem(pos,st,j))
        hi_possible = sat_add_i8(sums(j), maxrem(pos,st,j))
        if (j <= nlow) then
          if (hi_possible < lo_target) then
            ok = .false.
            return
          end if
        end if
        if (j > d - nup) then
          if (lo_possible > hi_target) then
            ok = .false.
            return
          end if
        end if
      end do
    end function completion_possible

  end subroutine search_md_i8

  subroutine build_envelope_1d(v, lower, upper, mn, mx)
    real(dp), intent(in) :: v(:)
    integer, intent(in) :: lower(:), upper(:)
    real(dp), allocatable, intent(out) :: mn(:,:), mx(:,:)
    integer :: n, len, pos, start, q
    real(dp) :: bestmn, bestmx, cmn, cmx
    n = size(v)
    len = size(lower)
    allocate(mn(len+1,n+2), mx(len+1,n+2))
    mn = rinf
    mx = -rinf
    mn(len+1,:) = 0.0_dp
    mx(len+1,:) = 0.0_dp
    do pos = len, 1, -1
      bestmn = rinf
      bestmx = -rinf
      do start = n, 1, -1
        q = start
        if (q >= lower(pos) .and. q <= upper(pos)) then
          if (mn(pos+1,q+1) < 0.5_dp*rinf) then
            cmn = v(q) + mn(pos+1,q+1)
            cmx = v(q) + mx(pos+1,q+1)
            bestmn = min(bestmn, cmn)
            bestmx = max(bestmx, cmx)
          end if
        end if
        mn(pos,start) = bestmn
        mx(pos,start) = bestmx
      end do
      mn(pos,n+1) = rinf
      mx(pos,n+1) = -rinf
    end do
  end subroutine build_envelope_1d

  subroutine build_envelope_md_real(v, lower, upper, mn, mx)
    real(dp), intent(in) :: v(:,:)
    integer, intent(in) :: lower(:), upper(:)
    real(dp), allocatable, intent(out) :: mn(:,:,:), mx(:,:,:)
    integer :: n, d, len, pos, start, q, j
    real(dp), allocatable :: bestmn(:), bestmx(:)
    n = size(v,1)
    d = size(v,2)
    len = size(lower)
    allocate(mn(len+1,n+2,d), mx(len+1,n+2,d), bestmn(d), bestmx(d))
    mn = rinf
    mx = -rinf
    mn(len+1,:,:) = 0.0_dp
    mx(len+1,:,:) = 0.0_dp
    do pos = len, 1, -1
      bestmn = rinf
      bestmx = -rinf
      do start = n, 1, -1
        q = start
        if (q >= lower(pos) .and. q <= upper(pos)) then
          if (mn(pos+1,q+1,1) < 0.5_dp*rinf) then
            do j = 1, d
              bestmn(j) = min(bestmn(j), v(q,j) + mn(pos+1,q+1,j))
              bestmx(j) = max(bestmx(j), v(q,j) + mx(pos+1,q+1,j))
            end do
          end if
        end if
        mn(pos,start,:) = bestmn
        mx(pos,start,:) = bestmx
      end do
      mn(pos,n+1,:) = rinf
      mx(pos,n+1,:) = -rinf
    end do
  end subroutine build_envelope_md_real

  subroutine build_envelope_md_i8(v, lower, upper, mn, mx)
    integer(i8), intent(in) :: v(:,:)
    integer, intent(in) :: lower(:), upper(:)
    integer(i8), allocatable, intent(out) :: mn(:,:,:), mx(:,:,:)
    integer :: n, d, len, pos, start, q, j
    integer(i8), allocatable :: bestmn(:), bestmx(:)
    integer(i8) :: cmn, cmx
    n = size(v,1)
    d = size(v,2)
    len = size(lower)
    allocate(mn(len+1,n+2,d), mx(len+1,n+2,d), bestmn(d), bestmx(d))
    mn = iinf
    mx = -iinf
    mn(len+1,:,:) = 0_i8
    mx(len+1,:,:) = 0_i8
    do pos = len, 1, -1
      bestmn = iinf
      bestmx = -iinf
      do start = n, 1, -1
        q = start
        if (q >= lower(pos) .and. q <= upper(pos)) then
          if (mn(pos+1,q+1,1) < iinf) then
            do j = 1, d
              cmn = sat_add_i8(v(q,j), mn(pos+1,q+1,j))
              cmx = sat_add_i8(v(q,j), mx(pos+1,q+1,j))
              bestmn(j) = min(bestmn(j), cmn)
              bestmx(j) = max(bestmx(j), cmx)
            end do
          end if
        end if
        mn(pos,start,:) = bestmn
        mx(pos,start,:) = bestmx
      end do
      mn(pos,n+1,:) = iinf
      mx(pos,n+1,:) = -iinf
    end do
  end subroutine build_envelope_md_i8

  logical function qualified_real(sums, target, me, nlow, nup) result(ok)
    real(dp), intent(in) :: sums(:), target(:), me(:)
    integer, intent(in) :: nlow, nup
    integer :: j, d
    d = size(sums)
    ok = .true.
    do j = 1, d
      if (j <= nlow) then
        if (sums(j) < target(j) - me(j)) then
          ok = .false.
          return
        end if
      end if
      if (j > d - nup) then
        if (sums(j) > target(j) + me(j)) then
          ok = .false.
          return
        end if
      end if
    end do
  end function qualified_real

  logical function qualified_i8(sums, target, me, nlow, nup) result(ok)
    integer(i8), intent(in) :: sums(:), target(:), me(:)
    integer, intent(in) :: nlow, nup
    integer :: j, d
    integer(i8) :: lo_target, hi_target
    d = size(sums)
    ok = .true.
    do j = 1, d
      lo_target = sat_add_i8(target(j), -me(j))
      hi_target = sat_add_i8(target(j), me(j))
      if (j <= nlow) then
        if (sums(j) < lo_target) then
          ok = .false.
          return
        end if
      end if
      if (j > d - nup) then
        if (sums(j) > hi_target) then
          ok = .false.
          return
        end if
      end if
    end do
  end function qualified_i8

  subroutine validate_bounds(lower, upper, n, where)
    integer, intent(in) :: lower(:), upper(:), n
    character(len=*), intent(in) :: where
    integer :: q
    do q = 1, size(lower)
      if (lower(q) < 1 .or. upper(q) > n .or. lower(q) > upper(q)) then
        error stop trim(where)//": invalid bounds"
      end if
    end do
    do q = 2, size(lower)
      if (lower(q) <= lower(q-1) .or. upper(q) <= upper(q-1)) then
        error stop trim(where)//": bounds must be strictly increasing by position"
      end if
    end do
  end subroutine validate_bounds

  pure integer(i8) function sat_add_i8(a, b) result(c)
    integer(i8), intent(in) :: a, b
    if (b > 0_i8) then
      if (a > iinf - b) then
        c = iinf
      else
        c = a + b
      end if
    else if (b < 0_i8) then
      if (a < -iinf - b) then
        c = -iinf
      else
        c = a + b
      end if
    else
      c = a
    end if
  end function sat_add_i8

end module flsss_search
