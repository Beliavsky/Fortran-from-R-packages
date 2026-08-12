! SPDX-License-Identifier: MIT
module optmatch_feasibility
   use optmatch_kinds, only : dp
   use optmatch_types, only : distance_spec, match_result
   use optmatch_matching, only : fullmatch
   implicit none
   private
   public :: caliper_size, caliper_upper_bound, max_caliper
   public :: max_controls_cap, min_controls_cap

contains

integer function caliper_size(scores, z, width, block) result(narcs)
   real(dp), intent(in) :: scores(:), width
   logical, intent(in) :: z(:)
   integer, intent(in), optional :: block(:)
   integer :: i, j
   if (size(scores) /= size(z)) error stop 'optmatch: scores/z length mismatch'
   if (width <= 0.0_dp) error stop 'optmatch: caliper width must be positive'
   if (present(block)) then
      if (size(block) /= size(z)) error stop 'optmatch: block/z length mismatch'
   end if
   narcs = 0
   do i = 1, size(scores)
      if (.not. z(i)) cycle
      do j = 1, size(scores)
         if (z(j)) cycle
         if (present(block)) then
            if (block(i) /= block(j)) cycle
         end if
         if (abs(scores(i) - scores(j)) <= width) narcs = narcs + 1
      end do
   end do
end function caliper_size

integer function caliper_upper_bound(scores, z, width, block) result(narcs)
   real(dp), intent(in) :: scores(:), width
   logical, intent(in) :: z(:)
   integer, intent(in), optional :: block(:)
   ! The R implementation uses a histogram-based upper bound.  Returning the
   ! exact count is also a valid (and sharper) upper bound, trading some speed
   ! for deterministic behavior and a simpler array API.
   if (present(block)) then
      narcs = caliper_size(scores, z, width, block)
   else
      narcs = caliper_size(scores, z, width)
   end if
end function caliper_upper_bound

real(dp) function max_caliper(scores, z, widths, max_problem_size, block, exact) result(width_out)
   real(dp), intent(in) :: scores(:), widths(:)
   logical, intent(in) :: z(:)
   integer, intent(in) :: max_problem_size
   integer, intent(in), optional :: block(:)
   logical, intent(in), optional :: exact
   real(dp), allocatable :: w(:)
   integer :: i, n
   logical :: use_exact
   if (size(widths) == 0) error stop 'optmatch: widths must not be empty'
   use_exact = .true.
   if (present(exact)) use_exact = exact
   w = widths
   call sort_descending(w)
   do i = 1, size(w)
      if (present(block)) then
         if (use_exact) then
            n = caliper_size(scores, z, w(i), block)
         else
            n = caliper_upper_bound(scores, z, w(i), block)
         end if
      else
         if (use_exact) then
            n = caliper_size(scores, z, w(i))
         else
            n = caliper_upper_bound(scores, z, w(i))
         end if
      end if
      if (n <= max_problem_size) then
         width_out = w(i)
         return
      end if
   end do
   error stop 'optmatch: no sufficiently small caliper width found'
end function max_caliper

real(dp) function max_controls_cap(distance, min_controls) result(cap)
   type(distance_spec), intent(in) :: distance
   real(dp), intent(in), optional :: min_controls
   real(dp), allocatable :: candidates(:)
   type(match_result) :: m
   real(dp) :: mn
   integer :: i
   mn = 0.0_dp
   if (present(min_controls)) mn = min_controls
   call ratio_candidates(size(distance%value, 1), size(distance%value, 2), candidates)
   cap = huge(1.0_dp)
   do i = 1, size(candidates)
      if (candidates(i) + sqrt(epsilon(1.0_dp)) < mn) cycle
      m = fullmatch(distance, min_controls=mn, max_controls=candidates(i))
      if (m%feasible) then
         cap = candidates(i)
         return
      end if
   end do
   m = fullmatch(distance, min_controls=mn)
   if (.not. m%feasible) cap = huge(1.0_dp)
end function max_controls_cap

real(dp) function min_controls_cap(distance, max_controls) result(cap)
   type(distance_spec), intent(in) :: distance
   real(dp), intent(in), optional :: max_controls
   real(dp), allocatable :: candidates(:)
   type(match_result) :: m
   real(dp) :: mx
   integer :: i
   mx = huge(1.0_dp)
   if (present(max_controls)) mx = max_controls
   call ratio_candidates(size(distance%value, 1), size(distance%value, 2), candidates)
   cap = 0.0_dp
   do i = size(candidates), 1, -1
      if (candidates(i) > mx + sqrt(epsilon(1.0_dp))) cycle
      m = fullmatch(distance, min_controls=candidates(i), max_controls=mx)
      if (m%feasible) then
         cap = candidates(i)
         return
      end if
   end do
end function min_controls_cap

subroutine ratio_candidates(nt, nc, values)
   integer, intent(in) :: nt, nc
   real(dp), allocatable, intent(out) :: values(:)
   real(dp), allocatable :: temp(:)
   integer :: k, n, i, j
   allocate(temp(max(1, nt) + max(1, nc)))
   n = 0
   do k = max(1, nt), 1, -1
      n = n + 1
      temp(n) = 1.0_dp / real(k, dp)
   end do
   do k = 2, max(1, nc)
      n = n + 1
      temp(n) = real(k, dp)
   end do
   call sort_ascending(temp(1:n))
   j = 0
   do i = 1, n
      if (j == 0) then
         j = 1
         temp(j) = temp(i)
      else if (abs(temp(i) - temp(j)) > 32.0_dp * epsilon(1.0_dp)) then
         j = j + 1
         temp(j) = temp(i)
      end if
   end do
   allocate(values(j))
   values = temp(1:j)
end subroutine ratio_candidates

subroutine sort_ascending(x)
   real(dp), intent(inout) :: x(:)
   integer :: i, j
   real(dp) :: key
   do i = 2, size(x)
      key = x(i)
      j = i - 1
      do while (j >= 1)
         if (x(j) <= key) exit
         x(j + 1) = x(j)
         j = j - 1
      end do
      x(j + 1) = key
   end do
end subroutine sort_ascending

subroutine sort_descending(x)
   real(dp), intent(inout) :: x(:)
   integer :: i, j
   real(dp) :: key
   do i = 2, size(x)
      key = x(i)
      j = i - 1
      do while (j >= 1)
         if (x(j) >= key) exit
         x(j + 1) = x(j)
         j = j - 1
      end do
      x(j + 1) = key
   end do
end subroutine sort_descending

end module optmatch_feasibility
