module flsss_packed
  use flsss_kinds, only : i8
  implicit none
  private

  integer, parameter :: max_lane_bits = 62
  integer(i8), parameter :: safe_i8 = 1152921504606846975_i8

  type, public :: packed_plan
    logical :: valid = .false.
    logical :: impossible = .false.
    integer :: d = 0
    integer :: nlane = 0
    integer, allocatable :: lane(:)
    integer, allocatable :: shift(:)
    integer, allocatable :: width(:)
    integer(i8), allocatable :: guard(:)
    integer(i8), allocatable :: lower_guard(:)
    integer(i8), allocatable :: upper_guard(:)
    integer(i8), allocatable :: lower_target(:)
    integer(i8), allocatable :: upper_target(:)
    integer(i8), allocatable :: maxsum(:)
  end type packed_plan

  public :: make_packed_plan, pack_matrix, pack_vector, packed_qualified
  public :: packed_meets_lower, packed_meets_upper
  public :: packed_completion_possible, packed_add_rows, is_comonotonic_i8
  public :: zero_minimum_i8, sort_rows_for_comonotonicity

contains

  function make_packed_plan(v, len, target, me, dl, du) result(plan)
    integer(i8), intent(in) :: v(:,:), target(:), me(:)
    integer, intent(in) :: len, dl, du
    type(packed_plan) :: plan
    integer :: d, j, lane, used, w, data_bits
    integer(i8) :: vmax, msum, lo, hi, bit, fieldmask
    logical :: lower_active, upper_active

    d = size(v,2)
    plan%d = d
    if (size(target) /= d .or. size(me) /= d) return
    if (len < 0 .or. len > size(v,1)) return
    if (any(v < 0_i8) .or. any(me < 0_i8)) return

    allocate(plan%lane(d), plan%shift(d), plan%width(d), plan%maxsum(d))
    lane = 1
    used = 0
    do j = 1, d
      vmax = maxval(v(:,j))
      if (.not. safe_mul_nonnegative(vmax, int(len,i8), msum)) return
      plan%maxsum(j) = msum
      data_bits = bits_needed(msum)
      w = data_bits + 1
      if (w > max_lane_bits) return
      if (used + w > max_lane_bits) then
        lane = lane + 1
        used = 0
      end if
      plan%lane(j) = lane
      plan%shift(j) = used
      plan%width(j) = w
      used = used + w
    end do
    plan%nlane = lane
    allocate(plan%guard(lane), plan%lower_guard(lane), plan%upper_guard(lane))
    allocate(plan%lower_target(lane), plan%upper_target(lane))
    plan%guard = 0_i8
    plan%lower_guard = 0_i8
    plan%upper_guard = 0_i8
    plan%lower_target = 0_i8
    plan%upper_target = 0_i8

    do j = 1, d
      lane = plan%lane(j)
      bit = shiftl(1_i8, plan%shift(j) + plan%width(j) - 1)
      plan%guard(lane) = ior(plan%guard(lane), bit)
      fieldmask = bit - 1_i8

      lo = sat_sub(target(j), me(j))
      hi = sat_add(target(j), me(j))
      lower_active = (j <= max(0,min(d,dl)))
      upper_active = (j > d - max(0,min(d,du)))

      if (lower_active) then
        if (lo > plan%maxsum(j)) then
          plan%impossible = .true.
          plan%valid = .true.
          return
        end if
        if (lo > 0_i8) then
          plan%lower_guard(lane) = ior(plan%lower_guard(lane), bit)
          plan%lower_target(lane) = ior(plan%lower_target(lane), &
            shiftl(iand(lo,fieldmask), plan%shift(j)))
        end if
      end if
      if (upper_active) then
        if (hi < 0_i8) then
          plan%impossible = .true.
          plan%valid = .true.
          return
        end if
        if (hi < plan%maxsum(j)) then
          plan%upper_guard(lane) = ior(plan%upper_guard(lane), bit)
          plan%upper_target(lane) = ior(plan%upper_target(lane), &
            shiftl(iand(max(0_i8,hi),fieldmask), plan%shift(j)))
        end if
      end if
    end do
    plan%valid = .true.
  end function make_packed_plan

  subroutine pack_matrix(plan, v, pv)
    type(packed_plan), intent(in) :: plan
    integer(i8), intent(in) :: v(:,:)
    integer(i8), allocatable, intent(out) :: pv(:,:)
    integer :: i, j, l
    integer(i8) :: x
    if (.not. plan%valid) error stop "pack_matrix: invalid plan"
    if (size(v,2) /= plan%d) error stop "pack_matrix: dimension mismatch"
    allocate(pv(size(v,1), plan%nlane))
    pv = 0_i8
    do j = 1, plan%d
      l = plan%lane(j)
      do i = 1, size(v,1)
        x = v(i,j)
        pv(i,l) = ior(pv(i,l), shiftl(x, plan%shift(j)))
      end do
    end do
  end subroutine pack_matrix

  subroutine pack_vector(plan, x, px)
    type(packed_plan), intent(in) :: plan
    integer(i8), intent(in) :: x(:)
    integer(i8), intent(out) :: px(:)
    integer :: j, l
    if (size(x) /= plan%d .or. size(px) /= plan%nlane) then
      error stop "pack_vector: dimension mismatch"
    end if
    px = 0_i8
    do j = 1, plan%d
      l = plan%lane(j)
      px(l) = ior(px(l), shiftl(x(j), plan%shift(j)))
    end do
  end subroutine pack_vector

  logical function packed_qualified(plan, sums) result(ok)
    type(packed_plan), intent(in) :: plan
    integer(i8), intent(in) :: sums(:)
    integer :: l
    integer(i8) :: z
    ok = .false.
    if (.not. plan%valid .or. plan%impossible) return
    if (size(sums) /= plan%nlane) return
    do l = 1, plan%nlane
      if (plan%lower_guard(l) /= 0_i8) then
        z = ior(sums(l), plan%guard(l)) - plan%lower_target(l)
        if (iand(z, plan%lower_guard(l)) /= plan%lower_guard(l)) return
      end if
      if (plan%upper_guard(l) /= 0_i8) then
        z = ior(plan%upper_target(l), plan%guard(l)) - sums(l)
        if (iand(z, plan%upper_guard(l)) /= plan%upper_guard(l)) return
      end if
    end do
    ok = .true.
  end function packed_qualified

  logical function packed_meets_lower(plan, sums) result(ok)
    type(packed_plan), intent(in) :: plan
    integer(i8), intent(in) :: sums(:)
    integer :: l
    integer(i8) :: z
    ok = .false.
    if (.not. plan%valid .or. plan%impossible) return
    do l = 1, plan%nlane
      if (plan%lower_guard(l) /= 0_i8) then
        z = ior(sums(l),plan%guard(l)) - plan%lower_target(l)
        if (iand(z,plan%lower_guard(l)) /= plan%lower_guard(l)) return
      end if
    end do
    ok = .true.
  end function packed_meets_lower

  logical function packed_meets_upper(plan, sums) result(ok)
    type(packed_plan), intent(in) :: plan
    integer(i8), intent(in) :: sums(:)
    integer :: l
    integer(i8) :: z
    ok = .false.
    if (.not. plan%valid .or. plan%impossible) return
    do l = 1, plan%nlane
      if (plan%upper_guard(l) /= 0_i8) then
        z = ior(plan%upper_target(l),plan%guard(l)) - sums(l)
        if (iand(z,plan%upper_guard(l)) /= plan%upper_guard(l)) return
      end if
    end do
    ok = .true.
  end function packed_meets_upper

  logical function packed_completion_possible(plan, sums, mn, mx) result(ok)
    type(packed_plan), intent(in) :: plan
    integer(i8), intent(in) :: sums(:), mn(:), mx(:)
    integer :: l
    integer(i8) :: z, possible
    ok = .false.
    if (.not. plan%valid .or. plan%impossible) return
    do l = 1, plan%nlane
      if (plan%lower_guard(l) /= 0_i8) then
        possible = sums(l) + mx(l)
        z = ior(possible, plan%guard(l)) - plan%lower_target(l)
        if (iand(z, plan%lower_guard(l)) /= plan%lower_guard(l)) return
      end if
      if (plan%upper_guard(l) /= 0_i8) then
        possible = sums(l) + mn(l)
        z = ior(plan%upper_target(l), plan%guard(l)) - possible
        if (iand(z, plan%upper_guard(l)) /= plan%upper_guard(l)) return
      end if
    end do
    ok = .true.
  end function packed_completion_possible

  subroutine packed_add_rows(pv, idx, sums)
    integer(i8), intent(in) :: pv(:,:)
    integer, intent(in) :: idx(:)
    integer(i8), intent(out) :: sums(:)
    integer :: i
    if (size(sums) /= size(pv,2)) error stop "packed_add_rows: dimension mismatch"
    sums = 0_i8
    do i = 1, size(idx)
      sums = sums + pv(idx(i),:)
    end do
  end subroutine packed_add_rows

  logical function is_comonotonic_i8(v) result(ok)
    integer(i8), intent(in) :: v(:,:)
    integer :: i, j
    ok = .true.
    do j = 1, size(v,2)
      do i = 2, size(v,1)
        if (v(i,j) < v(i-1,j)) then
          ok = .false.
          return
        end if
      end do
    end do
  end function is_comonotonic_i8

  subroutine zero_minimum_i8(v, len, target, shifted, shifted_target, ok)
    integer(i8), intent(in) :: v(:,:), target(:)
    integer, intent(in) :: len
    integer(i8), allocatable, intent(out) :: shifted(:,:), shifted_target(:)
    logical, intent(out) :: ok
    integer :: j
    integer(i8) :: mn, adj
    allocate(shifted(size(v,1),size(v,2)), shifted_target(size(target)))
    shifted = v
    shifted_target = target
    ok = .true.
    do j = 1, size(v,2)
      mn = minval(v(:,j))
      if (.not. safe_mul_signed(mn, int(len,i8), adj)) then
        ok = .false.
        return
      end if
      if (.not. safe_sub_signed(target(j), adj, shifted_target(j))) then
        ok = .false.
        return
      end if
      if (mn /= 0_i8) shifted(:,j) = v(:,j) - mn
      if (any(shifted(:,j) < 0_i8)) then
        ok = .false.
        return
      end if
    end do
  end subroutine zero_minimum_i8

  subroutine sort_rows_for_comonotonicity(v, order, sorted, success)
    integer(i8), intent(in) :: v(:,:)
    integer, allocatable, intent(out) :: order(:)
    integer(i8), allocatable, intent(out) :: sorted(:,:)
    logical, intent(out) :: success
    integer :: n, d, i, j, tmp, bestcol
    integer(i8) :: span, bestspan

    n = size(v,1)
    d = size(v,2)
    allocate(order(n), sorted(n,d))
    order = [(i, i=1,n)]
    bestcol = 1
    bestspan = -1_i8
    do j = 1, d
      span = maxval(v(:,j)) - minval(v(:,j))
      if (span > bestspan) then
        bestspan = span
        bestcol = j
      end if
    end do

    do i = 2, n
      tmp = order(i)
      j = i - 1
      do while (j >= 1)
        if (.not. row_less(v, tmp, order(j), bestcol)) exit
        order(j+1) = order(j)
        j = j - 1
      end do
      order(j+1) = tmp
    end do
    do i = 1, n
      sorted(i,:) = v(order(i),:)
    end do
    success = is_comonotonic_i8(sorted)

  contains
    logical function row_less(a, ia, ib, lead) result(less)
      integer(i8), intent(in) :: a(:,:)
      integer, intent(in) :: ia, ib, lead
      integer :: q
      if (a(ia,lead) < a(ib,lead)) then
        less = .true.
        return
      else if (a(ia,lead) > a(ib,lead)) then
        less = .false.
        return
      end if
      do q = 1, size(a,2)
        if (a(ia,q) < a(ib,q)) then
          less = .true.
          return
        else if (a(ia,q) > a(ib,q)) then
          less = .false.
          return
        end if
      end do
      less = ia < ib
    end function row_less
  end subroutine sort_rows_for_comonotonicity

  integer function bits_needed(x) result(n)
    integer(i8), intent(in) :: x
    integer(i8) :: y
    if (x <= 0_i8) then
      n = 1
      return
    end if
    n = 0
    y = x
    do while (y > 0_i8)
      n = n + 1
      y = shiftr(y,1)
    end do
  end function bits_needed

  logical function safe_mul_nonnegative(a, b, c) result(ok)
    integer(i8), intent(in) :: a, b
    integer(i8), intent(out) :: c
    ok = .false.
    c = 0_i8
    if (a < 0_i8 .or. b < 0_i8) return
    if (a /= 0_i8) then
      if (b > safe_i8 / a) return
    end if
    c = a * b
    ok = .true.
  end function safe_mul_nonnegative

  logical function safe_mul_signed(a, b, c) result(ok)
    integer(i8), intent(in) :: a, b
    integer(i8), intent(out) :: c
    ok = .false.
    c = 0_i8
    if (a == 0_i8 .or. b == 0_i8) then
      ok = .true.
      return
    end if
    if (a > 0_i8) then
      if (b > 0_i8) then
        if (a > safe_i8/b) return
      else
        if (b < -safe_i8/a) return
      end if
    else
      if (b > 0_i8) then
        if (a < -safe_i8/b) return
      else
        if (a < safe_i8/b) return
      end if
    end if
    c = a*b
    ok = .true.
  end function safe_mul_signed

  logical function safe_sub_signed(a, b, c) result(ok)
    integer(i8), intent(in) :: a, b
    integer(i8), intent(out) :: c
    if (b > 0_i8) then
      if (a < -safe_i8 + b) then
        ok = .false.; c = 0_i8; return
      end if
    else if (b < 0_i8) then
      if (a > safe_i8 + b) then
        ok = .false.; c = 0_i8; return
      end if
    end if
    c = a - b
    ok = .true.
  end function safe_sub_signed

  pure integer(i8) function sat_add(a,b) result(c)
    integer(i8), intent(in) :: a,b
    if (b > 0_i8 .and. a > safe_i8-b) then
      c = safe_i8
    else if (b < 0_i8 .and. a < -safe_i8-b) then
      c = -safe_i8
    else
      c = a+b
    end if
  end function sat_add

  pure integer(i8) function sat_sub(a,b) result(c)
    integer(i8), intent(in) :: a,b
    c = sat_add(a,-b)
  end function sat_sub

end module flsss_packed
