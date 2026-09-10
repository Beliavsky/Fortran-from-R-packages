! SPDX-License-Identifier: GPL-2.0-or-later
! Irregular-grid wavelet methods translated from wavethresh 4.7.3.
module wavethresh_irregular
   use r_kinds, only : dp
   use wavethresh_transform_1d, only : wd
   use wavethresh_types, only : grid_data_t, irregular_wavelet_t, wd_t
   implicit none
   private

   public :: irregwd, accessc, threshold_irregwd

contains

   function irregwd(grid, filter_number, family, boundary) result(object)
      !! Decomposes interpolated irregular observations and propagates their unit-noise variances.
      type(grid_data_t), intent(in) :: grid !! Valid interpolation map and regular-grid observations from makegrid.
      real(dp), intent(in), optional :: filter_number !! Wavethresh filter number; default is two.
      character(len=*), intent(in), optional :: family !! Real filter family; default is DaubExPhase.
      character(len=*), intent(in), optional :: boundary !! Boundary rule; the current dense backend supports periodic.
      type(irregular_wavelet_t) :: object
      type(wd_t) :: influence_transform
      real(dp), allocatable :: influence(:)
      real(dp) :: fnum
      real(dp) :: weight
      character(len=24) :: fam
      character(len=16) :: bc
      integer :: n_observations
      integer :: observation
      integer :: grid_index
      integer :: left_index
      integer :: level

      if (.not. grid%ok .or. .not. allocated(grid%grid_y) .or. &
         .not. allocated(grid%left_index) .or. .not. allocated(grid%weight_left)) then
         object%message = "grid must be a valid makegrid result"
         return
      end if
      if (size(grid%left_index) /= size(grid%grid_y) .or. &
         size(grid%weight_left) /= size(grid%grid_y)) then
         object%message = "grid interpolation arrays must conform"
         return
      end if
      fnum = 2.0_dp
      if (present(filter_number)) fnum = filter_number
      fam = "DaubExPhase"
      if (present(family)) fam = family
      bc = "periodic"
      if (present(boundary)) bc = boundary
      object%transform = wd(grid%grid_y, fnum, trim(fam), trim(bc))
      if (.not. object%transform%ok) then
         object%message = object%transform%message
         return
      end if

      n_observations = grid%n_observations
      if (n_observations < 1) n_observations = maxval(grid%left_index) + 1
      if (any(grid%left_index < 1) .or. any(grid%left_index >= n_observations)) then
         object%message = "grid interpolation indices are invalid"
         return
      end if
      allocate(object%coefficient_variance(0:object%transform%nlevels - 1))
      do level = 0, object%transform%nlevels - 1
         allocate(object%coefficient_variance(level)%values( &
            size(object%transform%detail(level)%values)), source=0.0_dp)
      end do

      allocate(influence(size(grid%grid_y)))
      do observation = 1, n_observations
         influence = 0.0_dp
         do grid_index = 1, size(influence)
            left_index = grid%left_index(grid_index)
            weight = grid%weight_left(grid_index)
            if (left_index == observation) influence(grid_index) = influence(grid_index) + weight
            if (left_index + 1 == observation) influence(grid_index) = influence(grid_index) + 1.0_dp - weight
         end do
         influence_transform = wd(influence, fnum, trim(fam), trim(bc))
         if (.not. influence_transform%ok) then
            object%message = influence_transform%message
            return
         end if
         do level = 0, object%transform%nlevels - 1
            object%coefficient_variance(level)%values = object%coefficient_variance(level)%values + &
               influence_transform%detail(level)%values**2
         end do
      end do
      object%ok = .true.
      object%message = "ok"
   end function irregwd

   pure function accessc(object, level, boundary) result(values)
      !! Returns coefficient noise-variance multipliers for one irregular-transform level.
      type(irregular_wavelet_t), intent(in) :: object !! Irregular wavelet decomposition.
      integer, intent(in) :: level !! R-style zero-based detail level.
      logical, intent(in), optional :: boundary !! Include boundary coefficients; retained for R API correspondence.
      real(dp), allocatable :: values(:)
      logical :: include_boundary

      include_boundary = .false.
      if (present(boundary)) include_boundary = boundary
      if (.not. object%ok .or. .not. allocated(object%coefficient_variance) .or. &
         level < 0 .or. level >= object%transform%nlevels) then
         allocate(values(0))
         return
      end if
      values = object%coefficient_variance(level)%values
      if (include_boundary) values = object%coefficient_variance(level)%values
   end function accessc

   pure function threshold_irregwd(input, levels, threshold_type, policy, value, by_level) result(output)
      !! Applies variance-adjusted manual or universal shrinkage to an irregular wavelet transform.
      type(irregular_wavelet_t), intent(in) :: input !! Irregular decomposition whose details are thresholded.
      integer, intent(in), optional :: levels(:) !! R-style detail levels; default is levels three through finest.
      character(len=*), intent(in), optional :: threshold_type !! Hard or soft shrinkage; default is hard.
      character(len=*), intent(in), optional :: policy !! Manual or universal threshold policy; default is universal.
      real(dp), intent(in), optional :: value !! Manual normalized threshold; default is zero.
      logical, intent(in), optional :: by_level !! Estimate noise and universal multipliers separately by level.
      type(irregular_wavelet_t) :: output
      integer, allocatable :: selected(:)
      real(dp), allocatable :: normalized(:)
      real(dp), allocatable :: variance(:)
      real(dp), allocatable :: cutoff(:)
      real(dp), allocatable :: multipliers(:)
      real(dp), allocatable :: noise(:)
      character(len=16) :: kind
      character(len=16) :: rule
      logical :: per_level
      integer :: first
      integer :: i
      integer :: level
      integer :: total

      output = input
      if (.not. input%ok) return
      kind = "hard"
      if (present(threshold_type)) kind = threshold_type
      rule = "universal"
      if (present(policy)) rule = policy
      if (trim(kind) /= "hard" .and. trim(kind) /= "soft") then
         output%ok = .false.
         output%message = "threshold type must be hard or soft"
         return
      end if
      if (trim(rule) /= "manual" .and. trim(rule) /= "universal") then
         output%ok = .false.
         output%message = "translated irregular threshold policies are manual and universal"
         return
      end if
      if (present(levels)) then
         selected = levels
      else
         first = min(3, input%transform%nlevels)
         if (first >= input%transform%nlevels) then
            allocate(selected(0))
         else
            allocate(selected(input%transform%nlevels - first))
            do i = 1, size(selected)
               selected(i) = first + i - 1
            end do
         end if
      end if
      if (any(selected < 0) .or. any(selected >= input%transform%nlevels)) then
         output%ok = .false.
         output%message = "threshold level is outside the detail range"
         return
      end if
      if (size(selected) == 0) return

      per_level = .false.
      if (present(by_level)) per_level = by_level
      allocate(noise(size(selected)), multipliers(size(selected)))
      if (per_level) then
         do i = 1, size(selected)
            noise(i) = level_noise(input, selected(i))
         end do
      else
         noise = pooled_noise(input, selected)
      end if
      if (trim(rule) == "manual") then
         multipliers = 0.0_dp
         if (present(value)) multipliers = max(value, 0.0_dp)
      else if (per_level) then
         do i = 1, size(selected)
            total = size(input%transform%detail(selected(i))%values)
            multipliers(i) = sqrt(2.0_dp * log(real(max(total, 2), dp)))
         end do
      else
         total = 0
         do i = 1, size(selected)
            total = total + size(input%transform%detail(selected(i))%values)
         end do
         multipliers = sqrt(2.0_dp * log(real(max(total, 2), dp)))
      end if

      do i = 1, size(selected)
         level = selected(i)
         variance = input%coefficient_variance(level)%values
         cutoff = multipliers(i) * noise(i) * sqrt(max(variance, 0.0_dp))
         normalized = output%transform%detail(level)%values
         if (trim(kind) == "hard") then
            where (abs(normalized) <= cutoff) normalized = 0.0_dp
         else
            normalized = sign(max(abs(normalized) - cutoff, 0.0_dp), normalized)
         end if
         output%transform%detail(level)%values = normalized
      end do
      output%message = "ok"
   end function threshold_irregwd

   pure function level_noise(object, level) result(noise)
      !! Estimates unit-observation noise from variance-standardized coefficients at one level.
      type(irregular_wavelet_t), intent(in) :: object !! Irregular decomposition supplying details and variances.
      integer, intent(in) :: level !! Valid zero-based detail level.
      real(dp) :: noise
      real(dp), allocatable :: values(:)
      real(dp), allocatable :: variance(:)
      logical, allocatable :: usable(:)

      variance = object%coefficient_variance(level)%values
      usable = variance > 1.0e-5_dp
      values = pack(object%transform%detail(level)%values / sqrt(max(variance, tiny(1.0_dp))), usable)
      noise = sqrt(max(sample_variance(values), 0.0_dp))
   end function level_noise

   pure function pooled_noise(object, levels) result(noise)
      !! Estimates one unit-observation noise scale across selected detail levels.
      type(irregular_wavelet_t), intent(in) :: object !! Irregular decomposition supplying details and variances.
      integer, intent(in) :: levels(:) !! Valid zero-based detail levels to pool.
      real(dp) :: noise
      real(dp), allocatable :: values(:)
      real(dp), allocatable :: variance(:)
      real(dp), allocatable :: level_values(:)
      integer :: total_count
      integer :: i
      integer :: position

      total_count = 0
      do i = 1, size(levels)
         total_count = total_count + count(object%coefficient_variance(levels(i))%values > 1.0e-5_dp)
      end do
      allocate(values(total_count))
      position = 1
      do i = 1, size(levels)
         variance = object%coefficient_variance(levels(i))%values
         level_values = pack(object%transform%detail(levels(i))%values / &
            sqrt(max(variance, tiny(1.0_dp))), variance > 1.0e-5_dp)
         if (size(level_values) > 0) then
            values(position:position + size(level_values) - 1) = level_values
            position = position + size(level_values)
         end if
      end do
      noise = sqrt(max(sample_variance(values), 0.0_dp))
   end function pooled_noise

   pure function sample_variance(values) result(variance)
      !! Computes the unbiased sample variance used by the upstream default dev=var path.
      real(dp), intent(in) :: values(:) !! Values whose sample variance is requested.
      real(dp) :: variance
      real(dp) :: center

      if (size(values) <= 1) then
         variance = 0.0_dp
         return
      end if
      center = sum(values) / real(size(values), dp)
      variance = sum((values - center)**2) / real(size(values) - 1, dp)
   end function sample_variance

end module wavethresh_irregular
