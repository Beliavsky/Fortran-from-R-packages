! SPDX-License-Identifier: MIT
module optmatch_utilities
   use optmatch_kinds, only : dp
   use optmatch_types, only : distance_spec, match_result, stratum_summary
   use optmatch_stats, only : mad_value, sample_variance
   implicit none
   private
   public :: standardization_scale, effective_sample_size, summarize_strata
   public :: matched_distances, matched_units, integer_digits

contains

real(dp) function standardization_scale(x, z, use_mad, supplied_scale) result(ans)
   real(dp), intent(in) :: x(:)
   logical, intent(in) :: z(:)
   logical, intent(in), optional :: use_mad
   real(dp), intent(in), optional :: supplied_scale
   real(dp) :: st, sc
   integer :: nt, nc, n
   logical :: robust
   if (size(x) /= size(z)) error stop 'optmatch: x/z length mismatch in standardization_scale'
   if (present(supplied_scale)) then
      ans = supplied_scale
      return
   end if
   n = size(x)
   nt = count(z)
   nc = n - nt
   if (n <= 2 .or. nt == 0 .or. nc == 0) error stop 'optmatch: both groups required for standardization scale'
   robust = .true.
   if (present(use_mad)) robust = use_mad
   if (robust) then
      st = mad_value(pack(x, z))
      sc = mad_value(pack(x, .not. z))
   else
      st = sqrt(sample_variance(pack(x, z)))
      sc = sqrt(sample_variance(pack(x, .not. z)))
   end if
   ans = sqrt((real(nt - 1, dp) * st**2 + real(nc - 1, dp) * sc**2) / real(n - 2, dp))
end function standardization_scale

real(dp) function effective_sample_size(match) result(ans)
   type(match_result), intent(in) :: match
   integer :: ng, g, nt, nc
   ans = 0.0_dp
   if (.not. allocated(match%treatment_group) .or. .not. allocated(match%control_group)) return
   ng = 0
   if (size(match%treatment_group) > 0) ng = max(ng, maxval(match%treatment_group))
   if (size(match%control_group) > 0) ng = max(ng, maxval(match%control_group))
   do g = 1, ng
      nt = count(match%treatment_group == g)
      nc = count(match%control_group == g)
      if (nt > 0 .and. nc > 0) ans = ans + 2.0_dp / (1.0_dp / real(nt, dp) + 1.0_dp / real(nc, dp))
   end do
end function effective_sample_size

function summarize_strata(match) result(out)
   type(match_result), intent(in) :: match
   type(stratum_summary) :: out
   integer, allocatable :: ttmp(:), ctmp(:), ftmp(:)
   integer :: ng, g, nt, nc, k, found
   ng = 0
   if (size(match%treatment_group) > 0) ng = max(ng, maxval(match%treatment_group))
   if (size(match%control_group) > 0) ng = max(ng, maxval(match%control_group))
   allocate(ttmp(max(1, ng)), ctmp(max(1, ng)), ftmp(max(1, ng)))
   ttmp = 0
   ctmp = 0
   ftmp = 0
   k = 0
   do g = 1, ng
      nt = count(match%treatment_group == g)
      nc = count(match%control_group == g)
      if (nt == 0 .and. nc == 0) cycle
      found = 0
      if (k > 0) then
         do found = 1, k
            if (ttmp(found) == nt .and. ctmp(found) == nc) exit
         end do
         if (found > k) found = 0
      end if
      if (found == 0) then
         k = k + 1
         ttmp(k) = nt
         ctmp(k) = nc
         ftmp(k) = 1
      else
         ftmp(found) = ftmp(found) + 1
      end if
   end do
   out%n_strata = k
   allocate(out%treatments(k), out%controls(k), out%frequency(k))
   if (k > 0) then
      out%treatments = ttmp(1:k)
      out%controls = ctmp(1:k)
      out%frequency = ftmp(1:k)
   end if
   out%effective_sample_size = effective_sample_size(match)
end function summarize_strata

subroutine matched_distances(x, match, distances, group_id)
   type(distance_spec), intent(in) :: x
   type(match_result), intent(in) :: match
   real(dp), allocatable, intent(out) :: distances(:)
   integer, allocatable, intent(out) :: group_id(:)
   integer :: i, j, n, k, g
   if (size(x%value, 1) /= size(match%treatment_group) .or. &
       size(x%value, 2) /= size(match%control_group)) error stop 'optmatch: match/distance dimension mismatch'
   n = 0
   do j = 1, size(x%value, 2)
      do i = 1, size(x%value, 1)
         g = match%treatment_group(i)
         if (g > 0 .and. match%control_group(j) == g .and. x%allowed(i, j)) n = n + 1
      end do
   end do
   allocate(distances(n), group_id(n))
   k = 0
   do j = 1, size(x%value, 2)
      do i = 1, size(x%value, 1)
         g = match%treatment_group(i)
         if (g <= 0 .or. match%control_group(j) /= g .or. .not. x%allowed(i, j)) cycle
         k = k + 1
         distances(k) = x%value(i, j)
         group_id(k) = g
      end do
   end do
end subroutine matched_distances

subroutine matched_units(match, treatment_matched, control_matched)
   type(match_result), intent(in) :: match
   logical, allocatable, intent(out) :: treatment_matched(:), control_matched(:)
   allocate(treatment_matched(size(match%treatment_group)), control_matched(size(match%control_group)))
   treatment_matched = match%treatment_group > 0
   control_matched = match%control_group > 0
end subroutine matched_units

integer function integer_digits(n) result(ans)
   integer, intent(in) :: n
   integer :: x
   if (n == 0) then
      ans = 1
      return
   end if
   if (n == -huge(n)-1) then
      ! Avoid overflow from abs(minimum integer).
      ans = range(n) + 2
      return
   end if
   x = abs(n)
   ans = 0
   do while (x > 0)
      ans = ans + 1
      x = x / 10
   end do
   if (n < 0) ans = ans + 1
end function integer_digits

end module optmatch_utilities
