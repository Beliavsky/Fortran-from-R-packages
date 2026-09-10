! SPDX-License-Identifier: GPL-2.0-or-later
! Wavelet shrinkage policies translated from wavethresh 4.7.3.
module wavethresh_threshold
   use wavethresh_types, only : dp, wd_t, wp_t, mwd_t, imwd_t, wd3d_t
   use wavethresh_stats, only : mad_value, madmad, robust_covariance, sure, shrink_soft, shrink_hard
   implicit none
   private

   public :: threshold_wd, threshold_wd_levels, threshold_wp, threshold_mwd, threshold_imwd, threshold_wd3d
   public :: threshold_value, universal_threshold

contains

   pure function universal_threshold(coefficients) result(threshold)
      real(dp), intent(in) :: coefficients(:) !! Coefficients used to estimate noise and determine a universal threshold.
      real(dp) :: threshold
      real(dp) :: noise
      if (size(coefficients) <= 1) then
         threshold = 0.0_dp
         return
      end if
      noise = sqrt(madmad(coefficients))
      threshold = sqrt(2.0_dp * log(real(size(coefficients), dp))) * noise
   end function universal_threshold

   elemental function threshold_value(coefficient, threshold, threshold_type) result(value)
      real(dp), intent(in) :: coefficient !! Wavelet coefficient to shrink.
      real(dp), intent(in) :: threshold !! Nonnegative threshold magnitude.
      character(len=*), intent(in) :: threshold_type !! hard or soft thresholding rule.
      real(dp) :: value
      select case (trim(threshold_type))
      case ("hard")
         value = shrink_hard(coefficient, threshold)
      case default
         value = shrink_soft(coefficient, threshold)
      end select
   end function threshold_value

   function threshold_wd(input, levels, threshold_type, policy, value, by_level) result(output)
      type(wd_t), intent(in) :: input !! One-dimensional wd or wst object to threshold.
      integer, intent(in), optional :: levels(:) !! R-style zero-based detail levels; default is levels 3 through finest.
      character(len=*), intent(in), optional :: threshold_type !! hard or soft; default is soft.
      character(len=*), intent(in), optional :: policy !! manual, universal, LSuniversal, or sure; default is sure.
      real(dp), intent(in), optional :: value !! Manual scalar threshold; default is zero.
      logical, intent(in), optional :: by_level !! Compute automatic thresholds separately by level when true.
      type(wd_t) :: output
      integer, allocatable :: selected(:)
      real(dp), allocatable :: thresholds(:)
      real(dp), allocatable :: pooled(:)
      real(dp), allocatable :: data(:)
      character(len=16) :: kind
      character(len=20) :: rule
      logical :: per_level
      real(dp) :: manual
      integer :: i
      integer :: level
      integer :: start
      integer :: total
      integer :: pos

      output = input
      if (.not. input%ok .or. input%nlevels <= 0) return
      kind = "soft"
      if (present(threshold_type)) kind = threshold_type
      rule = "sure"
      if (present(policy)) rule = policy
      manual = 0.0_dp
      if (present(value)) manual = value
      per_level = .false.
      if (present(by_level)) per_level = by_level
      if (present(levels)) then
         selected = levels
      else
         start = min(3, input%nlevels)
         if (start >= input%nlevels) then
            allocate(selected(0))
         else
            allocate(selected(input%nlevels - start))
            selected = [(i, i = start, input%nlevels - 1)]
         end if
      end if
      if (size(selected) == 0) return
      allocate(thresholds(size(selected)))
      select case (trim(rule))
      case ("manual", "mannum")
         thresholds = max(manual, 0.0_dp)
      case ("universal", "LSuniversal", "sure")
         if (per_level) then
            do i = 1, size(selected)
               level = selected(i)
               if (level < 0 .or. level >= input%nlevels) cycle
               data = input%detail(level)%values
               thresholds(i) = automatic_threshold(data, trim(rule))
            end do
         else
            total = 0
            do i = 1, size(selected)
               level = selected(i)
               if (level >= 0 .and. level < input%nlevels) total = total + size(input%detail(level)%values)
            end do
            allocate(pooled(total))
            pos = 1
            do i = 1, size(selected)
               level = selected(i)
               if (level < 0 .or. level >= input%nlevels) cycle
               data = input%detail(level)%values
               pooled(pos:pos + size(data) - 1) = data
               pos = pos + size(data)
            end do
            thresholds = automatic_threshold(pooled, trim(rule))
         end if
      case default
         output%message = "threshold policy not translated; object left unchanged"
         return
      end select
      do i = 1, size(selected)
         level = selected(i)
         if (level < 0 .or. level >= output%nlevels) cycle
         output%detail(level)%values = threshold_value(output%detail(level)%values, thresholds(i), trim(kind))
      end do
      output%message = "ok"
   end function threshold_wd

   function threshold_wd_levels(input, levels, values, threshold_type) result(output)
      type(wd_t), intent(in) :: input !! One-dimensional wd or wst object whose selected detail levels are thresholded.
      integer, intent(in) :: levels(:) !! Zero-based detail levels paired elementwise with values.
      real(dp), intent(in) :: values(:) !! Nonnegative manual threshold value for each selected detail level.
      character(len=*), intent(in), optional :: threshold_type !! hard or soft threshold rule; default is soft.
      type(wd_t) :: output
      character(len=16) :: kind
      integer :: i
      integer :: level

      output = input
      if (.not. input%ok) return
      if (size(levels) /= size(values) .or. any(values < 0.0_dp)) then
         output%ok = .false.
         output%message = "manual level thresholds require one nonnegative value per selected level"
         return
      end if
      kind = "soft"
      if (present(threshold_type)) kind = threshold_type
      if (trim(kind) /= "soft" .and. trim(kind) /= "hard") then
         output%ok = .false.
         output%message = "threshold type must be hard or soft"
         return
      end if
      do i = 1, size(levels)
         level = levels(i)
         if (level < 0 .or. level >= output%nlevels) then
            output%ok = .false.
            output%message = "manual level threshold index is outside the detail range"
            return
         end if
         output%detail(level)%values = threshold_value(output%detail(level)%values, values(i), trim(kind))
      end do
      output%message = "ok"
   end function threshold_wd_levels

   function threshold_wp(input, minimum_level, threshold, threshold_type) result(output)
      type(wp_t), intent(in) :: input !! Wavelet-packet tree whose packet coefficients are to be thresholded.
      integer, intent(in), optional :: minimum_level !! First packet-tree depth to shrink; default is one.
      real(dp), intent(in), optional :: threshold !! Manual threshold magnitude; default is zero.
      character(len=*), intent(in), optional :: threshold_type !! hard or soft; default is soft.
      type(wp_t) :: output
      integer :: first
      integer :: depth
      integer :: node
      real(dp) :: cutoff
      character(len=16) :: kind
      output = input
      if (.not. input%ok) return
      first = 1
      if (present(minimum_level)) first = minimum_level
      cutoff = 0.0_dp
      if (present(threshold)) cutoff = max(threshold, 0.0_dp)
      kind = "soft"
      if (present(threshold_type)) kind = threshold_type
      do depth = max(first, 1), input%nlevels
         do node = 0, 2**depth - 1
            output%packet(depth, node)%values = &
               threshold_value(output%packet(depth, node)%values, cutoff, trim(kind))
         end do
      end do
   end function threshold_wp

   function threshold_mwd(input, levels, threshold_type, policy, value, robust, bivariate, covtol) result(output)
      type(mwd_t), intent(in) :: input !! Multiple-wavelet decomposition whose selected detail levels are to be shrunk.
      integer, intent(in), optional :: levels(:) !! R-style zero-based detail levels; default is levels 3 through finest.
      character(len=*), intent(in), optional :: threshold_type !! hard or soft thresholding rule; default is hard.
      character(len=*), intent(in), optional :: policy !! manual, universal, or visushrink policy; default is universal.
      real(dp), intent(in), optional :: value(:) !! Manual threshold values; one scalar or one value per wavelet component.
      logical, intent(in), optional :: robust !! Use robust MAD/covariance estimates when true; default is true.
      logical, intent(in), optional :: bivariate !! Use joint chi-square shrinkage across wavelet components; default is true.
      real(dp), intent(in), optional :: covtol !! Covariance pivot tolerance for joint shrinkage; default is 1e-9.
      type(mwd_t) :: output
      integer, allocatable :: selected(:)
      real(dp), allocatable :: combined(:,:)
      real(dp), allocatable :: covariance(:,:)
      real(dp), allocatable :: inverse(:,:)
      real(dp), allocatable :: thresholds(:)
      real(dp) :: chi_square
      real(dp) :: cutoff
      real(dp) :: noise
      real(dp) :: pivot_tol
      real(dp) :: shrink
      character(len=16) :: kind
      character(len=20) :: rule
      logical :: joint
      logical :: robust_estimate
      logical :: invert_ok
      integer :: i
      integer :: j
      integer :: level
      integer :: n_total
      integer :: n_thresholdable
      integer :: position
      integer :: ncols
      integer :: npsi

      output = input
      if (.not. input%ok .or. input%nlevels <= 0) return
      npsi = input%filter%npsi
      if (npsi <= 0) then
         output%message = "threshold.mwd requires a valid multiple-wavelet filter"
         return
      end if
      kind = "hard"
      if (present(threshold_type)) kind = threshold_type
      rule = "universal"
      if (present(policy)) rule = policy
      robust_estimate = .true.
      if (present(robust)) robust_estimate = robust
      joint = .true.
      if (present(bivariate)) joint = bivariate
      pivot_tol = 1.0e-9_dp
      if (present(covtol)) pivot_tol = max(covtol, 0.0_dp)

      if (trim(kind) /= "hard" .and. trim(kind) /= "soft") then
         output%message = "threshold.mwd threshold_type must be hard or soft"
         return
      end if
      if (trim(rule) /= "manual" .and. trim(rule) /= "universal" .and. trim(rule) /= "visushrink") then
         output%message = "threshold.mwd policy must be manual, universal, or visushrink"
         return
      end if
      if (present(levels)) then
         selected = levels
      else if (input%nlevels > 3) then
         allocate(selected(input%nlevels - 3))
         selected = [(i, i = 3, input%nlevels - 1)]
      else
         allocate(selected(0))
      end if
      if (size(selected) == 0) return
      do i = 1, size(selected)
         if (selected(i) < 0 .or. selected(i) >= input%nlevels) then
            output%message = "threshold.mwd level is outside the detail range"
            return
         end if
      end do
      if (trim(rule) == "manual") then
         if (.not. present(value) .or. size(value) == 0 .or. any(value <= 0.0_dp)) then
            output%message = "manual threshold.mwd requires positive threshold values"
            return
         end if
      end if

      if (.not. joint) then
         n_total = 0
         do i = 1, size(selected)
            n_total = n_total + size(input%detail(selected(i))%values, 2)
         end do
         allocate(combined(npsi, n_total), source=0.0_dp)
         position = 1
         do i = 1, size(selected)
            level = selected(i)
            ncols = size(input%detail(level)%values, 2)
            combined(:, position:position + ncols - 1) = input%detail(level)%values
            position = position + ncols
         end do
         allocate(thresholds(npsi), source=0.0_dp)
         do j = 1, npsi
            if (trim(rule) == "manual") then
               if (size(value) == 1) then
                  thresholds(j) = value(1)
               else if (size(value) >= npsi) then
                  thresholds(j) = value(j)
               else
                  output%message = "manual threshold.mwd needs one threshold or one per wavelet component"
                  return
               end if
            else
               if (robust_estimate) then
                  noise = sqrt(max(mad_value(combined(j, :)), 0.0_dp))
               else
                  noise = sqrt(max(sample_variance(combined(j, :)), 0.0_dp))
               end if
               thresholds(j) = sqrt(2.0_dp * log(real(max(n_total, 1), dp))) * noise
               if (trim(rule) == "visushrink" .and. n_total > 0) then
                  thresholds(j) = thresholds(j) / sqrt(real(n_total, dp))
               end if
            end if
            combined(j, :) = threshold_value(combined(j, :), thresholds(j), trim(kind))
         end do
         position = 1
         do i = 1, size(selected)
            level = selected(i)
            ncols = size(output%detail(level)%values, 2)
            output%detail(level)%values = combined(:, position:position + ncols - 1)
            position = position + ncols
         end do
         output%message = "ok"
         return
      end if

      if (trim(rule) == "visushrink") then
         output%message = "bivariate threshold.mwd does not define the visushrink policy"
         return
      end if
      n_thresholdable = 0
      do i = 1, size(selected)
         level = selected(i)
         if (robust_estimate) then
            covariance = robust_covariance(input%detail(level)%values)
         else
            covariance = covariance_rows(input%detail(level)%values)
         end if
         call invert_small_matrix(covariance, pivot_tol, inverse, invert_ok)
         if (invert_ok) n_thresholdable = n_thresholdable + size(input%detail(level)%values, 2)
      end do
      if (trim(rule) == "manual") then
         cutoff = value(1)
      else if (n_thresholdable > 0) then
         cutoff = 2.0_dp * log(real(n_thresholdable, dp))
      else
         output%message = "no nonsingular multiple-wavelet levels were available for thresholding"
         return
      end if

      do i = 1, size(selected)
         level = selected(i)
         if (robust_estimate) then
            covariance = robust_covariance(input%detail(level)%values)
         else
            covariance = covariance_rows(input%detail(level)%values)
         end if
         call invert_small_matrix(covariance, pivot_tol, inverse, invert_ok)
         if (.not. invert_ok) cycle
         do j = 1, size(output%detail(level)%values, 2)
            chi_square = dot_product(output%detail(level)%values(:, j), &
               matmul(inverse, output%detail(level)%values(:, j)))
            if (trim(kind) == "hard") then
               if (chi_square < cutoff) output%detail(level)%values(:, j) = 0.0_dp
            else
               if (chi_square > 0.0_dp) then
                  shrink = max(chi_square - cutoff, 0.0_dp) / chi_square
               else
                  shrink = 0.0_dp
               end if
               output%detail(level)%values(:, j) = output%detail(level)%values(:, j) * shrink
            end if
         end do
      end do
      output%message = "ok"
   end function threshold_mwd

   function threshold_imwd(input, levels, threshold_type, policy, value, by_level) result(output)
      type(imwd_t), intent(in) :: input !! Two-dimensional decimated or stationary transform to threshold.
      integer, intent(in), optional :: levels(:) !! R-style zero-based 2-D detail levels; default is all levels.
      character(len=*), intent(in), optional :: threshold_type !! hard or soft; default is soft.
      character(len=*), intent(in), optional :: policy !! manual or universal; default is universal.
      real(dp), intent(in), optional :: value !! Manual scalar threshold; default is zero.
      logical, intent(in), optional :: by_level !! Compute universal thresholds separately by level when true.
      type(imwd_t) :: output
      integer, allocatable :: selected(:)
      real(dp), allocatable :: pooled(:)
      real(dp), allocatable :: data(:)
      real(dp), allocatable :: thresholds(:)
      character(len=16) :: kind
      character(len=20) :: rule
      logical :: per_level
      real(dp) :: manual
      integer :: i
      integer :: level
      integer :: total
      integer :: pos

      output = input
      if (.not. input%ok) return
      kind = "soft"
      if (present(threshold_type)) kind = threshold_type
      rule = "universal"
      if (present(policy)) rule = policy
      manual = 0.0_dp
      if (present(value)) manual = value
      per_level = .false.
      if (present(by_level)) per_level = by_level
      if (present(levels)) then
         selected = levels
      else
         allocate(selected(input%nlevels))
         selected = [(i, i = 0, input%nlevels - 1)]
      end if
      allocate(thresholds(size(selected)))
      if (trim(rule) == "manual") then
         thresholds = max(manual, 0.0_dp)
      else if (trim(rule) == "universal") then
         if (per_level) then
            do i = 1, size(selected)
               level = selected(i)
               if (level < 0 .or. level >= input%nlevels) cycle
               data = flatten_level(input, level)
               thresholds(i) = universal_threshold(data)
            end do
         else
            total = 0
            do i = 1, size(selected)
               level = selected(i)
               if (level < 0 .or. level >= input%nlevels) cycle
               total = total + 3 * size(input%level(level)%lh)
            end do
            allocate(pooled(total))
            pos = 1
            do i = 1, size(selected)
               level = selected(i)
               if (level < 0 .or. level >= input%nlevels) cycle
               data = flatten_level(input, level)
               pooled(pos:pos + size(data) - 1) = data
               pos = pos + size(data)
            end do
            thresholds = universal_threshold(pooled)
         end if
      else
         output%message = "2-D threshold policy not translated; object left unchanged"
         return
      end if
      do i = 1, size(selected)
         level = selected(i)
         if (level < 0 .or. level >= output%nlevels) cycle
         output%level(level)%lh = threshold_value(output%level(level)%lh, thresholds(i), trim(kind))
         output%level(level)%hl = threshold_value(output%level(level)%hl, thresholds(i), trim(kind))
         output%level(level)%hh = threshold_value(output%level(level)%hh, thresholds(i), trim(kind))
      end do
      output%message = "ok"
   end function threshold_imwd

   function threshold_wd3d(input, levels, threshold_type, policy, value, by_level) result(output)
      type(wd3d_t), intent(in) :: input !! Three-dimensional wavelet decomposition to threshold.
      integer, intent(in), optional :: levels(:) !! R-style zero-based 3-D detail levels; default is all levels.
      character(len=*), intent(in), optional :: threshold_type !! hard or soft; default is soft.
      character(len=*), intent(in), optional :: policy !! manual or universal; default is universal.
      real(dp), intent(in), optional :: value !! Manual scalar threshold; default is zero.
      logical, intent(in), optional :: by_level !! Compute universal thresholds separately by level when true.
      type(wd3d_t) :: output
      integer, allocatable :: selected(:)
      real(dp), allocatable :: thresholds(:)
      real(dp), allocatable :: pooled(:)
      real(dp), allocatable :: data(:)
      character(len=16) :: kind
      character(len=20) :: rule
      logical :: per_level
      real(dp) :: manual
      integer :: i
      integer :: level
      integer :: total
      integer :: pos

      output = input
      if (.not. input%ok) return
      kind = "soft"
      if (present(threshold_type)) kind = threshold_type
      rule = "universal"
      if (present(policy)) rule = policy
      manual = 0.0_dp
      if (present(value)) manual = value
      per_level = .false.
      if (present(by_level)) per_level = by_level
      if (present(levels)) then
         selected = levels
      else
         allocate(selected(input%nlevels))
         selected = [(i, i = 0, input%nlevels - 1)]
      end if
      allocate(thresholds(size(selected)))
      if (trim(rule) == "manual") then
         thresholds = max(manual, 0.0_dp)
      else if (trim(rule) == "universal") then
         if (per_level) then
            do i = 1, size(selected)
               level = selected(i)
               if (level < 0 .or. level >= input%nlevels) cycle
               data = reshape(input%level(level)%band, [size(input%level(level)%band)])
               thresholds(i) = universal_threshold(data)
            end do
         else
            total = 0
            do i = 1, size(selected)
               level = selected(i)
               if (level >= 0 .and. level < input%nlevels) total = total + size(input%level(level)%band)
            end do
            allocate(pooled(total))
            pos = 1
            do i = 1, size(selected)
               level = selected(i)
               if (level < 0 .or. level >= input%nlevels) cycle
               data = reshape(input%level(level)%band, [size(input%level(level)%band)])
               pooled(pos:pos + size(data) - 1) = data
               pos = pos + size(data)
            end do
            thresholds = universal_threshold(pooled)
         end if
      else
         output%message = "3-D threshold policy not translated; object left unchanged"
         return
      end if
      do i = 1, size(selected)
         level = selected(i)
         if (level < 0 .or. level >= output%nlevels) cycle
         output%level(level)%band = threshold_value(output%level(level)%band, thresholds(i), trim(kind))
      end do
      output%message = "ok"
   end function threshold_wd3d

   pure function automatic_threshold(data, policy) result(threshold)
      real(dp), intent(in) :: data(:) !! Coefficients from which an automatic threshold is selected.
      character(len=*), intent(in) :: policy !! universal, LSuniversal, or sure.
      real(dp) :: threshold
      real(dp) :: noise
      noise = sqrt(madmad(data))
      if (noise <= tiny(1.0_dp) .or. size(data) <= 1) then
         threshold = 0.0_dp
         return
      end if
      select case (trim(policy))
      case ("sure")
         threshold = noise * sure(data / noise)
      case ("LSuniversal")
         threshold = log(real(size(data), dp)) * noise
      case default
         threshold = sqrt(2.0_dp * log(real(size(data), dp))) * noise
      end select
   end function automatic_threshold

   pure function flatten_level(object, level) result(values)
      type(imwd_t), intent(in) :: object !! Two-dimensional transform providing the requested detail triplet.
      integer, intent(in) :: level !! Valid R-style zero-based detail level.
      real(dp), allocatable :: values(:)
      real(dp), allocatable :: one(:)
      real(dp), allocatable :: two(:)
      real(dp), allocatable :: three(:)
      one = reshape(object%level(level)%lh, [size(object%level(level)%lh)])
      two = reshape(object%level(level)%hl, [size(object%level(level)%hl)])
      three = reshape(object%level(level)%hh, [size(object%level(level)%hh)])
      values = [one, two, three]
   end function flatten_level

   pure function sample_variance(values) result(variance)
      real(dp), intent(in) :: values(:) !! Sample whose ordinary unbiased variance is required.
      real(dp) :: variance
      real(dp) :: center

      if (size(values) <= 1) then
         variance = 0.0_dp
         return
      end if
      center = sum(values) / real(size(values), dp)
      variance = sum((values - center)**2) / real(size(values) - 1, dp)
   end function sample_variance

   pure function covariance_rows(values) result(covariance)
      real(dp), intent(in) :: values(:,:) !! Matrix with wavelet components in rows and observations in columns.
      real(dp), allocatable :: covariance(:,:)
      real(dp), allocatable :: centered(:,:)
      real(dp), allocatable :: means(:)
      integer :: i

      allocate(covariance(size(values, 1), size(values, 1)), source=0.0_dp)
      if (size(values, 2) <= 1) return
      allocate(means(size(values, 1)))
      do i = 1, size(values, 1)
         means(i) = sum(values(i, :)) / real(size(values, 2), dp)
      end do
      allocate(centered(size(values, 1), size(values, 2)))
      do i = 1, size(values, 1)
         centered(i, :) = values(i, :) - means(i)
      end do
      covariance = matmul(centered, transpose(centered)) / real(size(values, 2) - 1, dp)
   end function covariance_rows

   pure subroutine invert_small_matrix(matrix, tolerance, inverse, ok)
      real(dp), intent(in) :: matrix(:,:) !! Square covariance matrix to invert by pivoted Gauss-Jordan elimination.
      real(dp), intent(in) :: tolerance !! Absolute pivot threshold below which the matrix is treated as singular.
      real(dp), allocatable, intent(out) :: inverse(:,:) !! Matrix inverse when successful; empty on invalid shape.
      logical, intent(out) :: ok !! True when every pivot exceeded the supplied tolerance and inversion succeeded.
      real(dp), allocatable :: work(:,:)
      real(dp), allocatable :: temp_row(:)
      real(dp) :: factor
      real(dp) :: pivot
      integer :: i
      integer :: j
      integer :: n
      integer :: pivot_row

      ok = .false.
      if (size(matrix, 1) /= size(matrix, 2) .or. size(matrix, 1) == 0) then
         allocate(inverse(0, 0))
         return
      end if
      n = size(matrix, 1)
      allocate(work(n, 2 * n), source=0.0_dp)
      work(:, 1:n) = matrix
      do i = 1, n
         work(i, n + i) = 1.0_dp
      end do
      allocate(temp_row(2 * n))
      do i = 1, n
         pivot_row = i
         do j = i + 1, n
            if (abs(work(j, i)) > abs(work(pivot_row, i))) pivot_row = j
         end do
         if (abs(work(pivot_row, i)) <= tolerance) then
            allocate(inverse(0, 0))
            return
         end if
         if (pivot_row /= i) then
            temp_row = work(i, :)
            work(i, :) = work(pivot_row, :)
            work(pivot_row, :) = temp_row
         end if
         pivot = work(i, i)
         work(i, :) = work(i, :) / pivot
         do j = 1, n
            if (j == i) cycle
            factor = work(j, i)
            work(j, :) = work(j, :) - factor * work(i, :)
         end do
      end do
      inverse = work(:, n + 1:2 * n)
      ok = .true.
   end subroutine invert_small_matrix


end module wavethresh_threshold
