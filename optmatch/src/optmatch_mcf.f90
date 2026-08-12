! SPDX-License-Identifier: MIT
module optmatch_mcf
   use optmatch_kinds, only : dp, optmatch_inf
   implicit none
   private
   public :: min_cost_flow_integer

contains

subroutine min_cost_flow_integer(source, target, capacity, cost, supply, flow, total_cost, feasible)
   integer, intent(in) :: source(:), target(:), capacity(:), supply(:)
   real(dp), intent(in) :: cost(:)
   integer, allocatable, intent(out) :: flow(:)
   real(dp), intent(out) :: total_cost
   logical, intent(out) :: feasible
   integer :: n, m, extra, ss, tt, narcs, e, i, required, sent, delta
   integer, allocatable :: u(:), v(:), cap(:), rev(:), orig(:), parent(:)
   real(dp), allocatable :: c(:), dist(:)
   logical :: found
   if (size(source) /= size(target) .or. size(source) /= size(capacity) .or. size(source) /= size(cost)) &
      error stop 'optmatch: inconsistent min-cost-flow arc arrays'
   if (any(capacity < 0)) error stop 'optmatch: negative arc capacity'
   if (sum(supply) /= 0) error stop 'optmatch: min-cost-flow supplies must sum to zero'
   n = size(supply)
   m = size(source)
   if (m > 0) then
      if (minval(source) < 1 .or. minval(target) < 1 .or. maxval(source) > n .or. maxval(target) > n) &
         error stop 'optmatch: min-cost-flow node index out of range'
   end if
   extra = count(supply > 0) + count(supply < 0)
   ss = n + 1
   tt = n + 2
   narcs = 2 * (m + extra)
   allocate(u(narcs), v(narcs), cap(narcs), c(narcs), rev(narcs), orig(narcs))
   e = 0
   do i = 1, m
      call add_residual_arc(source(i), target(i), capacity(i), cost(i), i, e, u, v, cap, c, rev, orig)
   end do
   required = 0
   do i = 1, n
      if (supply(i) > 0) then
         call add_residual_arc(ss, i, supply(i), 0.0_dp, 0, e, u, v, cap, c, rev, orig)
         required = required + supply(i)
      else if (supply(i) < 0) then
         call add_residual_arc(i, tt, -supply(i), 0.0_dp, 0, e, u, v, cap, c, rev, orig)
      end if
   end do
   allocate(flow(m), parent(n + 2), dist(n + 2))
   flow = 0
   sent = 0
   total_cost = 0.0_dp
   do while (sent < required)
      call shortest_residual_path(ss, tt, u, v, cap, c, n + 2, dist, parent, found)
      if (.not. found) exit
      delta = required - sent
      i = tt
      do while (i /= ss)
         e = parent(i)
         if (e == 0) exit
         delta = min(delta, cap(e))
         i = u(e)
      end do
      if (i /= ss .or. delta <= 0) exit
      i = tt
      do while (i /= ss)
         e = parent(i)
         cap(e) = cap(e) - delta
         cap(rev(e)) = cap(rev(e)) + delta
         if (orig(e) > 0) then
            flow(orig(e)) = flow(orig(e)) + delta
         else if (orig(rev(e)) > 0) then
            flow(orig(rev(e))) = flow(orig(rev(e))) - delta
         end if
         total_cost = total_cost + real(delta, dp) * c(e)
         i = u(e)
      end do
      sent = sent + delta
   end do
   feasible = sent == required
   if (.not. feasible) then
      flow = 0
      total_cost = 0.0_dp
   end if
end subroutine min_cost_flow_integer

subroutine add_residual_arc(a, b, capacity, cost, original_id, e, u, v, cap, c, rev, orig)
   integer, intent(in) :: a, b, capacity, original_id
   real(dp), intent(in) :: cost
   integer, intent(inout) :: e
   integer, intent(inout) :: u(:), v(:), cap(:), rev(:), orig(:)
   real(dp), intent(inout) :: c(:)
   integer :: f, r
   f = e + 1
   r = e + 2
   u(f) = a
   v(f) = b
   cap(f) = capacity
   c(f) = cost
   rev(f) = r
   orig(f) = original_id
   u(r) = b
   v(r) = a
   cap(r) = 0
   c(r) = -cost
   rev(r) = f
   orig(r) = 0
   e = r
end subroutine add_residual_arc

subroutine shortest_residual_path(source, sink, u, v, cap, cost, n, dist, parent, found)
   integer, intent(in) :: source, sink, u(:), v(:), cap(:), n
   real(dp), intent(in) :: cost(:)
   real(dp), intent(out) :: dist(:)
   integer, intent(out) :: parent(:)
   logical, intent(out) :: found
   integer :: iter, e
   logical :: changed
   real(dp), parameter :: eps = 1.0e-12_dp
   dist = optmatch_inf
   parent = 0
   dist(source) = 0.0_dp
   do iter = 1, n - 1
      changed = .false.
      do e = 1, size(u)
         if (cap(e) <= 0) cycle
         if (dist(u(e)) >= optmatch_inf / 2.0_dp) cycle
         if (dist(v(e)) > dist(u(e)) + cost(e) + eps) then
            dist(v(e)) = dist(u(e)) + cost(e)
            parent(v(e)) = e
            changed = .true.
         end if
      end do
      if (.not. changed) exit
   end do
   found = parent(sink) /= 0 .or. source == sink
end subroutine shortest_residual_path

end module optmatch_mcf
