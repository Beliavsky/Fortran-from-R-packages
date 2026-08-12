! SPDX-License-Identifier: MIT
module optmatch_matching
   use optmatch_kinds, only : dp
   use optmatch_types, only : distance_spec, match_result
   use optmatch_distance, only : find_subproblems
   use optmatch_mcf, only : min_cost_flow_integer
   implicit none
   private
   public :: fullmatch, pairmatch, fmatch_core, matching_cost

contains

function fullmatch(x, min_controls, max_controls, omit_fraction, mean_controls) result(out)
   type(distance_spec), intent(in) :: x
   real(dp), intent(in), optional :: min_controls, max_controls, omit_fraction, mean_controls
   type(match_result) :: out
   integer, allocatable :: tc(:), cc(:), ti(:), ci(:)
   type(distance_spec) :: sub
   type(match_result) :: part
   integer :: ncomp, k, i, j, nt, nc, ngroup
   real(dp) :: mn, mx, omf, mean_c
   logical :: has_omf, has_mean, all_ok

   nt = size(x%value, 1)
   nc = size(x%value, 2)
   call init_result(out, nt, nc)
   if (nt == 0 .or. nc == 0) return
   if (any(x%value < 0.0_dp .and. x%allowed)) error stop 'optmatch: distances must be nonnegative'

   mn = 0.0_dp
   mx = huge(1.0_dp)
   if (present(min_controls)) mn = min_controls
   if (present(max_controls)) mx = max_controls
   if (mn < 0.0_dp .or. mx <= 0.0_dp .or. mn > mx) error stop 'optmatch: invalid control ratio bounds'
   has_omf = present(omit_fraction)
   has_mean = present(mean_controls)
   if (has_omf .and. has_mean) error stop 'optmatch: omit_fraction and mean_controls are mutually exclusive'
   omf = 0.0_dp
   if (has_omf) then
      omf = omit_fraction
      if (abs(omf) > 1.0_dp) error stop 'optmatch: omit_fraction must be in [-1,1]'
   end if
   mean_c = 0.0_dp
   if (has_mean) then
      mean_c = mean_controls
      if (mean_c <= 0.0_dp) error stop 'optmatch: mean_controls must be positive'
      if (mean_c < mn - sqrt(epsilon(1.0_dp)) .or. mean_c > mx + sqrt(epsilon(1.0_dp))) &
         error stop 'optmatch: mean_controls outside min/max controls'
   end if

   call find_subproblems(x, tc, cc, ncomp)
   ngroup = 0
   all_ok = .true.
   do k = 1, ncomp
      ti = pack([(i, i = 1, nt)], tc == k)
      ci = pack([(j, j = 1, nc)], cc == k)
      if (size(ti) == 0 .or. size(ci) == 0) cycle
      call extract_distance(x, ti, ci, sub)
      if (has_mean) then
         omf = 1.0_dp - mean_c * real(size(ti), dp) / real(size(ci), dp)
         part = fullmatch_single(sub, mn, mx, .true., omf, size(ti), size(ci))
      else
         part = fullmatch_single(sub, mn, mx, has_omf, omf, size(ti), size(ci))
      end if
      if (.not. part%feasible) all_ok = .false.
      call merge_part(out, part, ti, ci, ngroup)
   end do
   out%feasible = all_ok .and. ncomp > 0
   out%n_selected = count(out%selected)
   out%objective = matching_cost(x, out)
end function fullmatch

function pairmatch(x, controls, remove_unmatchables) result(out)
   type(distance_spec), intent(in) :: x
   integer, intent(in), optional :: controls
   logical, intent(in), optional :: remove_unmatchables
   type(match_result) :: out
   integer, allocatable :: tc(:), cc(:), ti(:), ci(:)
   type(distance_spec) :: sub
   type(match_result) :: part
   integer :: kctl, nt, nc, ncomp, k, i, j, ngroup
   real(dp) :: omf, control_overage, treatment_overage
   logical :: remove_um, all_ok
   nt = size(x%value, 1)
   nc = size(x%value, 2)
   call init_result(out, nt, nc)
   kctl = 1
   if (present(controls)) kctl = controls
   if (kctl < 1) error stop 'optmatch: pairmatch controls must be >= 1'
   remove_um = .false.
   if (present(remove_unmatchables)) remove_um = remove_unmatchables
   if (.not. remove_um) then
      do i = 1, nt
         if (.not. any(x%allowed(i, :))) return
      end do
   end if
   call find_subproblems(x, tc, cc, ncomp)
   ngroup = 0
   all_ok = .true.
   do k = 1, ncomp
      ti = pack([(i, i = 1, nt)], tc == k)
      ci = pack([(j, j = 1, nc)], cc == k)
      if (size(ti) == 0 .or. size(ci) == 0) cycle
      call extract_distance(x, ti, ci, sub)
      control_overage = real(size(ci) - kctl * size(ti), dp)
      treatment_overage = real(size(ti), dp) - real(size(ci), dp) / real(kctl, dp)
      if (control_overage >= 0.0_dp) then
         omf = control_overage / real(size(ci), dp)
      else
         omf = -treatment_overage / real(size(ti), dp)
      end if
      part = fullmatch_single(sub, real(kctl, dp), real(kctl, dp), .true., omf, size(ti), size(ci))
      if (.not. part%feasible) all_ok = .false.
      call merge_part(out, part, ti, ci, ngroup)
   end do
   out%feasible = all_ok .and. ncomp > 0
   out%n_selected = count(out%selected)
   out%objective = matching_cost(x, out)
end function pairmatch

function fullmatch_single(x, min_controls, max_controls, has_omit, omit_fraction, &
                          original_nt, original_nc) result(out)
   type(distance_spec), intent(in) :: x
   real(dp), intent(in) :: min_controls, max_controls, omit_fraction
   logical, intent(in) :: has_omit
   integer, intent(in) :: original_nt, original_nc
   type(match_result) :: out
   integer, allocatable :: ti(:), ci(:)
   type(distance_spec) :: reduced, oriented
   type(match_result) :: core
   integer :: i, j, ntf, ncf
   real(dp) :: mn, mx, omf, adjusted
   logical :: flip, use_omit

   call init_result(out, size(x%value, 1), size(x%value, 2))
   ti = pack([(i, i = 1, size(x%value, 1))], [(any(x%allowed(i, :)), i = 1, size(x%value, 1))])
   ci = pack([(j, j = 1, size(x%value, 2))], [(any(x%allowed(:, j)), j = 1, size(x%value, 2))])
   if (size(ti) == 0 .or. size(ci) == 0) return
   call extract_distance(x, ti, ci, reduced)

   flip = .false.
   if (has_omit) then
      flip = omit_fraction < 0.0_dp
   else
      flip = max_controls <= 0.5_dp
   end if

   if (.not. flip) then
      mx = min(max_controls, real(size(reduced%value, 2), dp))
      mn = max(min_controls, 1.0_dp / real(size(reduced%value, 1), dp))
      oriented = reduced
      ntf = size(reduced%value, 1)
      ncf = size(reduced%value, 2)
      omf = omit_fraction
      use_omit = has_omit
      if (use_omit .and. omf > 0.0_dp .and. ncf < original_nc) then
         adjusted = (omf * real(original_nc, dp) - real(original_nc - ncf, dp)) / real(ncf, dp)
         if (adjusted <= 0.0_dp) then
            use_omit = .false.
         else
            omf = adjusted
         end if
      end if
      core = fmatch_core(oriented, mn, mx, use_omit, max(0.0_dp, omf))
      if (core%feasible) then
         do j = 1, size(ci)
            do i = 1, size(ti)
               if (core%selected(i, j)) out%selected(ti(i), ci(j)) = .true.
            end do
         end do
      end if
   else
      mx = min(reciprocal_or_huge(min_controls), real(size(reduced%value, 1), dp))
      mn = max(1.0_dp / max_controls, 1.0_dp / real(size(reduced%value, 2), dp))
      call transpose_distance(reduced, oriented)
      ntf = size(oriented%value, 1)
      ncf = size(oriented%value, 2)
      omf = -omit_fraction
      use_omit = has_omit
      if (use_omit .and. omf > 0.0_dp .and. ncf < original_nt) then
         adjusted = (omf * real(original_nt, dp) - real(original_nt - ncf, dp)) / real(ncf, dp)
         if (adjusted <= 0.0_dp) then
            use_omit = .false.
         else
            omf = adjusted
         end if
      end if
      core = fmatch_core(oriented, mn, mx, use_omit, max(0.0_dp, omf))
      if (core%feasible) then
         do j = 1, size(ti)
            do i = 1, size(ci)
               if (core%selected(i, j)) out%selected(ti(j), ci(i)) = .true.
            end do
         end do
      end if
   end if

   if (core%feasible) then
      out%feasible = .true.
      call groups_from_selected(out%selected, out%treatment_group, out%control_group)
      out%n_selected = count(out%selected)
      out%objective = matching_cost(x, out)
   end if
end function fullmatch_single

function fmatch_core(distance, min_cpt, max_cpt, has_omit, omit_fraction) result(out)
   type(distance_spec), intent(in) :: distance
   real(dp), intent(in) :: min_cpt, max_cpt, omit_fraction
   logical, intent(in) :: has_omit
   type(match_result) :: out
   integer :: nt, nc, mxc, mnc, mxr, nmc, narcs, m, e, i, j, end_id, sink_id
   integer, allocatable :: source(:), target(:), capacity(:), supply(:), flow(:)
   real(dp), allocatable :: cost(:)
   real(dp) :: total_cost, f
   logical :: feasible

   nt = size(distance%value, 1)
   nc = size(distance%value, 2)
   call init_result(out, nt, nc)
   if (nt == 0 .or. nc == 0) return
   if (min_cpt <= 0.0_dp .or. max_cpt <= 0.0_dp) error stop 'optmatch: min/max cpt must be positive'
   if (max_cpt < min_cpt) error stop 'optmatch: min cpt exceeds max cpt'

   mxc = ceiling(max_cpt)
   mnc = max(1, floor(min_cpt))
   mxr = ceiling(1.0_dp / min_cpt)
   if (mnc > 1) mxr = 1
   if (mxc < mnc .or. min(mxc, mnc, mxr) < 1) error stop 'optmatch: invalid integral matching bounds'
   f = 1.0_dp
   if (has_omit) then
      if (omit_fraction < 0.0_dp .or. omit_fraction > 1.0_dp) error stop 'optmatch: oriented omit fraction outside [0,1]'
      f = 1.0_dp - omit_fraction
   end if
   nmc = r_round_nonnegative(real(nc, dp) * f)

   if ((mxr > 1 .and. real(nt, dp) / real(mxr, dp) > real(nmc, dp)) .or. &
       (mxr == 1 .and. nt * mnc > nmc) .or. nt * mxc < nmc) return

   narcs = count(distance%allowed)
   m = narcs + nt + 2 * nc
   allocate(source(m), target(m), capacity(m), cost(m))
   allocate(supply(nt + nc + 2))
   end_id = nt + nc + 1
   sink_id = end_id + 1
   e = 0
   do j = 1, nc
      do i = 1, nt
         if (.not. distance%allowed(i, j)) cycle
         e = e + 1
         source(e) = i
         target(e) = nt + j
         capacity(e) = 1
         cost(e) = distance%value(i, j)
      end do
   end do
   do i = 1, nt
      e = e + 1
      source(e) = i
      target(e) = end_id
      capacity(e) = mxc - mnc
      cost(e) = 0.0_dp
   end do
   do j = 1, nc
      e = e + 1
      source(e) = nt + j
      target(e) = end_id
      capacity(e) = mxr - 1
      cost(e) = 0.0_dp
   end do
   do j = 1, nc
      e = e + 1
      source(e) = nt + j
      target(e) = sink_id
      capacity(e) = 1
      cost(e) = 0.0_dp
   end do
   if (e /= m) error stop 'optmatch: internal arc count mismatch'
   supply = 0
   supply(1:nt) = mxc
   supply(end_id) = -(mxc * nt - nmc)
   supply(sink_id) = -nmc
   call min_cost_flow_integer(source, target, capacity, cost, supply, flow, total_cost, feasible)
   if (.not. feasible) return

   e = 0
   do j = 1, nc
      do i = 1, nt
         if (.not. distance%allowed(i, j)) cycle
         e = e + 1
         if (flow(e) > 0) out%selected(i, j) = .true.
      end do
   end do
   out%feasible = .true.
   out%objective = total_cost
   out%n_selected = count(out%selected)
   call groups_from_selected(out%selected, out%treatment_group, out%control_group)
end function fmatch_core

real(dp) function matching_cost(x, match) result(ans)
   type(distance_spec), intent(in) :: x
   type(match_result), intent(in) :: match
   if (.not. allocated(match%selected)) then
      ans = 0.0_dp
   else
      ans = sum(x%value, mask=match%selected)
   end if
end function matching_cost

subroutine init_result(out, nt, nc)
   type(match_result), intent(out) :: out
   integer, intent(in) :: nt, nc
   allocate(out%treatment_group(nt), out%control_group(nc), out%selected(nt, nc))
   out%treatment_group = 0
   out%control_group = 0
   out%selected = .false.
   out%feasible = .false.
   out%objective = 0.0_dp
   out%n_selected = 0
end subroutine init_result

subroutine extract_distance(x, ti, ci, sub)
   type(distance_spec), intent(in) :: x
   integer, intent(in) :: ti(:), ci(:)
   type(distance_spec), intent(out) :: sub
   integer :: i, j
   allocate(sub%value(size(ti), size(ci)), sub%allowed(size(ti), size(ci)))
   do j = 1, size(ci)
      do i = 1, size(ti)
         sub%value(i, j) = x%value(ti(i), ci(j))
         sub%allowed(i, j) = x%allowed(ti(i), ci(j))
      end do
   end do
end subroutine extract_distance

subroutine transpose_distance(x, xt)
   type(distance_spec), intent(in) :: x
   type(distance_spec), intent(out) :: xt
   allocate(xt%value(size(x%value, 2), size(x%value, 1)))
   allocate(xt%allowed(size(x%allowed, 2), size(x%allowed, 1)))
   xt%value = transpose(x%value)
   xt%allowed = transpose(x%allowed)
end subroutine transpose_distance

subroutine merge_part(out, part, ti, ci, ngroup)
   type(match_result), intent(inout) :: out
   type(match_result), intent(in) :: part
   integer, intent(in) :: ti(:), ci(:)
   integer, intent(inout) :: ngroup
   integer :: i, j, gmax
   if (.not. part%feasible) return
   gmax = 0
   if (size(part%treatment_group) > 0) gmax = max(gmax, maxval(part%treatment_group))
   if (size(part%control_group) > 0) gmax = max(gmax, maxval(part%control_group))
   do i = 1, size(ti)
      if (part%treatment_group(i) > 0) out%treatment_group(ti(i)) = ngroup + part%treatment_group(i)
   end do
   do j = 1, size(ci)
      if (part%control_group(j) > 0) out%control_group(ci(j)) = ngroup + part%control_group(j)
   end do
   do j = 1, size(ci)
      do i = 1, size(ti)
         if (part%selected(i, j)) out%selected(ti(i), ci(j)) = .true.
      end do
   end do
   ngroup = ngroup + gmax
end subroutine merge_part

subroutine groups_from_selected(selected, tg, cg)
   logical, intent(in) :: selected(:, :)
   integer, intent(out) :: tg(:), cg(:)
   integer, allocatable :: parent(:), root_to_group(:)
   integer :: nt, nc, i, j, r, ng
   nt = size(selected, 1)
   nc = size(selected, 2)
   allocate(parent(nt + nc), root_to_group(nt + nc))
   parent = [(i, i = 1, nt + nc)]
   root_to_group = 0
   do j = 1, nc
      do i = 1, nt
         if (selected(i, j)) call uf_union(parent, i, nt + j)
      end do
   end do
   ng = 0
   tg = 0
   cg = 0
   do i = 1, nt
      if (.not. any(selected(i, :))) cycle
      r = uf_find(parent, i)
      if (root_to_group(r) == 0) then
         ng = ng + 1
         root_to_group(r) = ng
      end if
      tg(i) = root_to_group(r)
   end do
   do j = 1, nc
      if (.not. any(selected(:, j))) cycle
      r = uf_find(parent, nt + j)
      if (root_to_group(r) == 0) then
         ng = ng + 1
         root_to_group(r) = ng
      end if
      cg(j) = root_to_group(r)
   end do
end subroutine groups_from_selected

recursive integer function uf_find(parent, x) result(r)
   integer, intent(inout) :: parent(:)
   integer, intent(in) :: x
   if (parent(x) == x) then
      r = x
   else
      parent(x) = uf_find(parent, parent(x))
      r = parent(x)
   end if
end function uf_find

subroutine uf_union(parent, a, b)
   integer, intent(inout) :: parent(:)
   integer, intent(in) :: a, b
   integer :: ra, rb
   ra = uf_find(parent, a)
   rb = uf_find(parent, b)
   if (ra /= rb) parent(rb) = ra
end subroutine uf_union

pure real(dp) function reciprocal_or_huge(x) result(ans)
   real(dp), intent(in) :: x
   if (x <= 0.0_dp) then
      ans = huge(1.0_dp)
   else
      ans = 1.0_dp / x
   end if
end function reciprocal_or_huge

integer function r_round_nonnegative(x) result(ans)
   real(dp), intent(in) :: x
   integer :: lo
   real(dp) :: frac, tol
   if (x < 0.0_dp) error stop 'optmatch: internal negative rounding input'
   lo = floor(x)
   frac = x - real(lo, dp)
   tol = 16.0_dp * epsilon(1.0_dp) * max(1.0_dp, abs(x))
   if (frac < 0.5_dp - tol) then
      ans = lo
   else if (frac > 0.5_dp + tol) then
      ans = lo + 1
   else if (mod(lo, 2) == 0) then
      ans = lo
   else
      ans = lo + 1
   end if
end function r_round_nonnegative

end module optmatch_matching
