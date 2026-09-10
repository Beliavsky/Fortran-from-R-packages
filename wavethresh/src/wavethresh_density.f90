! SPDX-License-Identifier: GPL-2.0-or-later
! Claw distribution utilities translated from wavethresh 4.7.3.
module wavethresh_density
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   use wavethresh_types, only : dp, density_grid_t, density_projection_t, density_wavelet_t
   use wavethresh_types, only : first_last_t, scaling_function_t, support_t, wt_filter_t
   use wavethresh_basis, only : wavelet_support
   use wavethresh_filters, only : filter_select
   use wavethresh_utilities, only : first_last_dh
   implicit none
   private
   public :: chires5, chires6, cwavde, dencvwd, denproj, denwd, evaluate_density
   public :: dclaw, pclaw, rclaw
contains

   pure function chires5(x, tau, resolution_level, filter_number, family, n_iterations) result(projection)
      !! Estimates high-resolution density scaling coefficients using the Daubechies-Lagarias algorithm.
      real(dp), intent(in) :: x(:) !! Finite observations to project onto translated scaling functions.
      real(dp), intent(in), optional :: tau !! Positive resolution multiplier; default is one.
      integer, intent(in) :: resolution_level !! Dyadic resolution level J.
      real(dp), intent(in), optional :: filter_number !! Wavethresh real-filter number; default is 10.
      character(len=*), intent(in), optional :: family !! Wavethresh real-filter family; default is DaubLeAsymm.
      integer, intent(in), optional :: n_iterations !! Binary-product iterations; default is 20.
      type(density_projection_t) :: projection

      projection = project_density(x, tau, resolution_level, filter_number, family, n_iterations, .false.)
   end function chires5

   pure function chires6(x, tau, resolution_level, filter_number, family, n_iterations) result(projection)
      !! Estimates high-resolution density coefficients and their upper covariance bands.
      real(dp), intent(in) :: x(:) !! Finite observations to project onto translated scaling functions.
      real(dp), intent(in), optional :: tau !! Positive resolution multiplier; default is one.
      integer, intent(in) :: resolution_level !! Dyadic resolution level J.
      real(dp), intent(in), optional :: filter_number !! Wavethresh real-filter number; default is 10.
      character(len=*), intent(in), optional :: family !! Wavethresh real-filter family; default is DaubLeAsymm.
      integer, intent(in), optional :: n_iterations !! Binary-product iterations; default is 20.
      type(density_projection_t) :: projection

      projection = project_density(x, tau, resolution_level, filter_number, family, n_iterations, .true.)
   end function chires6

   pure function denproj(x, resolution_level, tau, filter_number, family, covariance, n_iterations) result(projection)
      !! Dispatches density projection with optional covariance-band calculation.
      real(dp), intent(in) :: x(:) !! Finite observations to project onto translated scaling functions.
      integer, intent(in) :: resolution_level !! Dyadic resolution level J.
      real(dp), intent(in), optional :: tau !! Positive resolution multiplier; default is one.
      real(dp), intent(in), optional :: filter_number !! Wavethresh real-filter number; default is 10.
      character(len=*), intent(in), optional :: family !! Wavethresh real-filter family; default is DaubLeAsymm.
      logical, intent(in), optional :: covariance !! Whether to calculate covariance bands; default is false.
      integer, intent(in), optional :: n_iterations !! Binary-product iterations; default is 20.
      type(density_projection_t) :: projection
      logical :: include_covariance

      include_covariance = .false.
      if (present(covariance)) include_covariance = covariance
      projection = project_density(x, tau, resolution_level, filter_number, family, n_iterations, include_covariance)
   end function denproj

   pure function evaluate_density(projection, x, n_iterations) result(grid)
      !! Evaluates a high-resolution density projection at caller-supplied locations.
      type(density_projection_t), intent(in) :: projection !! Valid projection returned by chires5, chires6, or denproj.
      real(dp), intent(in) :: x(:) !! Finite locations at which to evaluate the projected density.
      integer, intent(in), optional :: n_iterations !! Binary-product iterations; defaults to the projection value.
      type(density_grid_t) :: grid
      real(dp), allocatable :: phi_values(:)
      real(dp) :: z
      integer :: iterations
      integer :: local_min
      integer :: local_max
      integer :: point
      integer :: k
      integer :: phi_index

      if (.not. projection%ok) then
         grid%message = "evaluate_density requires a valid density projection"
         return
      end if
      if (any(.not. ieee_is_finite(x))) then
         grid%message = "density evaluation locations must be finite"
         return
      end if
      iterations = projection%iterations
      if (present(n_iterations)) iterations = n_iterations
      if (iterations < 1) then
         grid%message = "n_iterations must be positive"
         return
      end if

      grid%x = x
      allocate(grid%y(size(x)), source=0.0_dp)
      do point = 1, size(x)
         z = projection%primary_resolution * x(point)
         local_min = ceiling(z - real(size(projection%filter%low) - 1, dp))
         local_max = floor(z)
         phi_values = daubechies_lagarias_phi(z, projection%filter%low, iterations)
         do k = max(local_min, projection%k_min), min(local_max, projection%k_max)
            phi_index = k - local_min + 1
            if (phi_index <= size(phi_values)) then
               grid%y(point) = grid%y(point) + projection%coefficients(k - projection%k_min + 1) * &
                  sqrt(projection%primary_resolution) * phi_values(phi_index)
            end if
         end do
      end do
      grid%ok = .true.
      grid%message = "ok"
   end function evaluate_density

   pure function cwavde(x, j_max, sf, wv, threshold, n_output, primary_resolution, filter_number, family) result(grid)
      !! Estimates a density from sampled scaling and wavelet functions using the upstream CWavDE kernel.
      real(dp), intent(in) :: x(:) !! Finite observations from the target density.
      integer, intent(in) :: j_max !! Number of wavelet resolution levels; must be nonnegative.
      type(scaling_function_t), intent(in) :: sf !! Sampled scaling function on an increasing regular grid.
      type(scaling_function_t), intent(in) :: wv !! Sampled wavelet function on an increasing regular grid.
      real(dp), intent(in), optional :: threshold !! Nonnegative hard threshold for wavelet coefficients; default is zero.
      integer, intent(in), optional :: n_output !! Number of density grid points; default is 100 and minimum is two.
      real(dp), intent(in), optional :: primary_resolution !! Positive primary resolution; default is one.
      real(dp), intent(in), optional :: filter_number !! Wavethresh real-filter number used for support; default is 10.
      character(len=*), intent(in), optional :: family !! Wavethresh real-filter family used for support.
      type(density_grid_t) :: grid
      type(support_t) :: support
      real(dp), allocatable :: coefficients(:)
      real(dp) :: selected_threshold
      real(dp) :: resolution
      real(dp) :: selected_filter_number
      real(dp) :: divisor
      real(dp) :: coefficient
      real(dp) :: lower_bound
      real(dp) :: upper_bound
      character(len=24) :: selected_family
      integer :: output_size
      integer :: k_min
      integer :: k_max
      integer :: level
      integer :: k
      integer :: observation
      integer :: point

      if (size(x) == 0 .or. any(.not. ieee_is_finite(x))) then
         grid%message = "CWavDE requires at least one finite observation"
         return
      end if
      if (j_max < 0) then
         grid%message = "j_max must be nonnegative"
         return
      end if
      if (.not. valid_sampled_function(sf) .or. .not. valid_sampled_function(wv)) then
         grid%message = "CWavDE requires valid sampled scaling and wavelet functions"
         return
      end if
      selected_threshold = 0.0_dp
      if (present(threshold)) selected_threshold = threshold
      if (.not. ieee_is_finite(selected_threshold) .or. selected_threshold < 0.0_dp) then
         grid%message = "threshold must be finite and nonnegative"
         return
      end if
      output_size = 100
      if (present(n_output)) output_size = n_output
      if (output_size < 2) then
         grid%message = "n_output must be at least two"
         return
      end if
      resolution = 1.0_dp
      if (present(primary_resolution)) resolution = primary_resolution
      if (.not. ieee_is_finite(resolution) .or. resolution <= 0.0_dp) then
         grid%message = "primary_resolution must be finite and positive"
         return
      end if
      selected_filter_number = 10.0_dp
      if (present(filter_number)) selected_filter_number = filter_number
      selected_family = "DaubLeAsymm"
      if (present(family)) selected_family = family
      support = wavelet_support(selected_filter_number, trim(selected_family))
      if (.not. support%ok) then
         grid%message = support%message
         return
      end if

      k_min = floor(minval(x) - support%phi_right / resolution)
      k_max = ceiling(maxval(x) - support%phi_left / resolution)
      allocate(grid%scaling_indices(k_max - k_min + 1))
      do k = k_min, k_max
         grid%scaling_indices(k - k_min + 1) = k
      end do
      allocate(grid%wavelet_k_min(j_max), grid%wavelet_k_max(j_max))
      lower_bound = real(k_min, dp) + support%phi_left / resolution
      upper_bound = real(k_max, dp) + support%phi_right / resolution
      do level = 1, j_max
         divisor = resolution * 2.0_dp**level
         grid%wavelet_k_min(level) = floor(minval(x) - support%psi_right / divisor)
         grid%wavelet_k_max(level) = ceiling(maxval(x) - support%psi_left / divisor)
         lower_bound = min(lower_bound, real(grid%wavelet_k_min(level), dp) + support%psi_left / divisor)
         upper_bound = max(upper_bound, real(grid%wavelet_k_max(level), dp) + support%psi_right / divisor)
      end do

      allocate(grid%x(output_size), grid%y(output_size), source=0.0_dp)
      do point = 1, output_size
         grid%x(point) = lower_bound + real(point - 1, dp) * (upper_bound - lower_bound) / real(output_size - 1, dp)
      end do
      allocate(coefficients(k_max - k_min + 1), source=0.0_dp)
      do k = k_min, k_max
         do observation = 1, size(x)
            coefficients(k - k_min + 1) = coefficients(k - k_min + 1) + &
               sampled_value(sf, resolution * x(observation) - real(k, dp))
         end do
         coefficients(k - k_min + 1) = sqrt(resolution) * coefficients(k - k_min + 1) / real(size(x), dp)
         do point = 1, output_size
            grid%y(point) = grid%y(point) + sqrt(resolution) * coefficients(k - k_min + 1) * &
               sampled_value(sf, resolution * grid%x(point) - real(k, dp))
         end do
      end do

      do level = 1, j_max
         divisor = resolution * 2.0_dp**level
         do k = grid%wavelet_k_min(level), grid%wavelet_k_max(level)
            coefficient = 0.0_dp
            do observation = 1, size(x)
               coefficient = coefficient + sampled_value(wv, divisor * x(observation) - real(k, dp))
            end do
            coefficient = sqrt(divisor) * coefficient / real(size(x), dp)
            if (abs(coefficient) <= selected_threshold) cycle
            do point = 1, output_size
               grid%y(point) = grid%y(point) + sqrt(divisor) * coefficient * &
                  sampled_value(wv, divisor * grid%x(point) - real(k, dp))
            end do
         end do
      end do
      grid%ok = .true.
      grid%message = "ok"
   end function cwavde

   pure function valid_sampled_function(sampled_function) result(valid)
      !! Checks the shape, finiteness, and regular increasing grid required by the upstream interpolator.
      type(scaling_function_t), intent(in) :: sampled_function !! Sampled function descriptor to validate.
      logical :: valid
      real(dp) :: spacing
      real(dp) :: tolerance
      integer :: i

      valid = sampled_function%ok
      if (.not. valid) return
      valid = allocated(sampled_function%x) .and. allocated(sampled_function%y)
      if (.not. valid) return
      valid = size(sampled_function%x) == size(sampled_function%y) .and. size(sampled_function%x) >= 3
      if (.not. valid) return
      valid = all(ieee_is_finite(sampled_function%x)) .and. all(ieee_is_finite(sampled_function%y))
      if (.not. valid) return
      spacing = sampled_function%x(2) - sampled_function%x(1)
      valid = spacing > 0.0_dp
      if (.not. valid) return
      tolerance = 64.0_dp * epsilon(1.0_dp) * max(1.0_dp, abs(spacing))
      do i = 3, size(sampled_function%x)
         if (abs(sampled_function%x(i) - sampled_function%x(i - 1) - spacing) > tolerance) then
            valid = .false.
            return
         end if
      end do
   end function valid_sampled_function

   pure function sampled_value(sampled_function, x) result(value)
      !! Reproduces the upstream CWavDE linear interpolator, including its length-minus-two index scaling.
      type(scaling_function_t), intent(in) :: sampled_function !! Valid sampled function on a regular grid.
      real(dp), intent(in) :: x !! Coordinate at which to interpolate, returning zero outside the grid.
      real(dp) :: value
      real(dp) :: fractional_index
      real(dp) :: fraction
      real(dp) :: width
      integer :: left

      value = 0.0_dp
      if (x < sampled_function%x(1) .or. x > sampled_function%x(size(sampled_function%x))) return
      width = sampled_function%x(size(sampled_function%x)) - sampled_function%x(1)
      fractional_index = real(size(sampled_function%x) - 2, dp) * (x - sampled_function%x(1)) / width
      left = int(fractional_index) + 1
      fraction = fractional_index - real(left - 1, dp)
      value = (1.0_dp - fraction) * sampled_function%y(left) + fraction * sampled_function%y(left + 1)
   end function sampled_value

   function denwd(projection) result(transform)
      !! Applies the upstream zero-boundary density wavelet decomposition to projection coefficients.
      type(density_projection_t), intent(in) :: projection !! Valid density projection whose coefficients are decomposed.
      type(density_wavelet_t) :: transform
      type(first_last_t) :: bounds
      real(dp) :: coefficient
      integer :: level
      integer :: k
      integer :: m
      integer :: source_index
      integer :: first_input

      if (.not. projection%ok) then
         transform%message = "denwd requires a valid density projection"
         return
      end if
      bounds = first_last_dh(size(projection%filter%low), size(projection%coefficients), &
         transform_type="wavelet", boundary="zero", firstk=[projection%k_min, projection%k_max])
      if (.not. bounds%ok) then
         transform%message = bounds%message
         return
      end if
      transform%nlevels = size(bounds%detail, 1)
      transform%filter = projection%filter
      allocate(transform%scaling(0:transform%nlevels))
      allocate(transform%detail(0:transform%nlevels - 1))
      transform%scaling_first = bounds%scaling(:, 1)
      transform%scaling_last = bounds%scaling(:, 2)
      transform%detail_first = bounds%detail(:, 1)
      transform%detail_last = bounds%detail(:, 2)
      transform%scaling(transform%nlevels)%values = projection%coefficients

      do level = transform%nlevels - 1, 0, -1
         first_input = transform%scaling_first(level + 2)
         allocate(transform%scaling(level)%values(transform%scaling_last(level + 1) - &
            transform%scaling_first(level + 1) + 1), source=0.0_dp)
         allocate(transform%detail(level)%values(transform%detail_last(level + 1) - &
            transform%detail_first(level + 1) + 1), source=0.0_dp)
         do k = transform%scaling_first(level + 1), transform%scaling_last(level + 1)
            coefficient = 0.0_dp
            do m = 0, size(transform%filter%low) - 1
               source_index = m + 2 * k - first_input + 1
               if (source_index < 1 .or. source_index > size(transform%scaling(level + 1)%values)) cycle
               coefficient = coefficient + transform%filter%low(m + 1) * &
                  transform%scaling(level + 1)%values(source_index)
            end do
            transform%scaling(level)%values(k - transform%scaling_first(level + 1) + 1) = coefficient
         end do
         do k = transform%detail_first(level + 1), transform%detail_last(level + 1)
            coefficient = 0.0_dp
            do m = 0, size(transform%filter%low) - 1
               source_index = 2 * k + 1 - m - first_input + 1
               if (source_index < 1 .or. source_index > size(transform%scaling(level + 1)%values)) cycle
               if (mod(m, 2) == 0) then
                  coefficient = coefficient - transform%filter%low(m + 1) * &
                     transform%scaling(level + 1)%values(source_index)
               else
                  coefficient = coefficient + transform%filter%low(m + 1) * &
                     transform%scaling(level + 1)%values(source_index)
               end if
            end do
            transform%detail(level)%values(k - transform%detail_first(level + 1) + 1) = coefficient
         end do
      end do
      transform%ok = .true.
      transform%message = "ok"
   end function denwd

   function dencvwd(projection) result(transform)
      !! Propagates Chires6 covariance bands through the zero-boundary density wavelet transform.
      type(density_projection_t), intent(in) :: projection !! Valid density projection containing covariance bands.
      type(density_wavelet_t) :: transform
      type(first_last_t) :: bounds
      real(dp), allocatable :: covariance(:,:)
      real(dp), allocatable :: next_covariance(:,:)
      real(dp), allocatable :: detail_covariance(:,:)
      real(dp), allocatable :: scaling_operator(:,:)
      real(dp), allocatable :: detail_operator(:,:)
      integer :: n_coefficients
      integer :: n_scaling
      integer :: n_detail
      integer :: band
      integer :: level
      integer :: row
      integer :: k
      integer :: m
      integer :: source_index
      integer :: first_input

      if (.not. projection%ok .or. .not. projection%has_covariance) then
         transform%message = "dencvwd requires a valid Chires6 covariance projection"
         return
      end if
      n_coefficients = size(projection%coefficients)
      if (size(projection%covariance, 1) /= n_coefficients) then
         transform%message = "density covariance row count does not match the coefficient count"
         return
      end if
      bounds = first_last_dh(size(projection%filter%low), n_coefficients, transform_type="wavelet", &
         boundary="zero", firstk=[projection%k_min, projection%k_max])
      if (.not. bounds%ok) then
         transform%message = bounds%message
         return
      end if
      transform%nlevels = size(bounds%detail, 1)
      transform%filter = projection%filter
      allocate(transform%scaling(0:transform%nlevels))
      allocate(transform%detail(0:transform%nlevels - 1))
      transform%scaling_first = bounds%scaling(:, 1)
      transform%scaling_last = bounds%scaling(:, 2)
      transform%detail_first = bounds%detail(:, 1)
      transform%detail_last = bounds%detail(:, 2)

      allocate(covariance(n_coefficients, n_coefficients), source=0.0_dp)
      do band = 1, size(projection%covariance, 2)
         do row = 1, n_coefficients - band + 1
            covariance(row, row + band - 1) = projection%covariance(row, band)
            covariance(row + band - 1, row) = projection%covariance(row, band)
         end do
      end do
      allocate(transform%scaling(transform%nlevels)%values(n_coefficients))
      do row = 1, n_coefficients
         transform%scaling(transform%nlevels)%values(row) = covariance(row, row)
      end do

      do level = transform%nlevels - 1, 0, -1
         first_input = transform%scaling_first(level + 2)
         n_scaling = transform%scaling_last(level + 1) - transform%scaling_first(level + 1) + 1
         n_detail = transform%detail_last(level + 1) - transform%detail_first(level + 1) + 1
         allocate(scaling_operator(n_scaling, size(covariance, 1)), source=0.0_dp)
         allocate(detail_operator(n_detail, size(covariance, 1)), source=0.0_dp)
         do k = transform%scaling_first(level + 1), transform%scaling_last(level + 1)
            row = k - transform%scaling_first(level + 1) + 1
            do m = 0, size(transform%filter%low) - 1
               source_index = m + 2 * k - first_input + 1
               if (source_index < 1 .or. source_index > size(covariance, 1)) cycle
               scaling_operator(row, source_index) = transform%filter%low(m + 1)
            end do
         end do
         do k = transform%detail_first(level + 1), transform%detail_last(level + 1)
            row = k - transform%detail_first(level + 1) + 1
            do m = 0, size(transform%filter%low) - 1
               source_index = 2 * k + 1 - m - first_input + 1
               if (source_index < 1 .or. source_index > size(covariance, 1)) cycle
               detail_operator(row, source_index) = merge(-1.0_dp, 1.0_dp, mod(m, 2) == 0) * &
                  transform%filter%low(m + 1)
            end do
         end do
         next_covariance = matmul(scaling_operator, matmul(covariance, transpose(scaling_operator)))
         detail_covariance = matmul(detail_operator, matmul(covariance, transpose(detail_operator)))
         allocate(transform%scaling(level)%values(n_scaling))
         allocate(transform%detail(level)%values(n_detail))
         do row = 1, n_scaling
            transform%scaling(level)%values(row) = next_covariance(row, row)
         end do
         do row = 1, n_detail
            transform%detail(level)%values(row) = detail_covariance(row, row)
         end do
         call move_alloc(next_covariance, covariance)
         deallocate(scaling_operator, detail_operator, detail_covariance)
      end do
      transform%ok = .true.
      transform%message = "ok"
   end function dencvwd

   pure function project_density(x, tau, resolution_level, filter_number, family, n_iterations, with_covariance) result(projection)
      !! Implements the common Chires5/Chires6 empirical scaling-coefficient calculation.
      real(dp), intent(in) :: x(:) !! Finite observations to project onto translated scaling functions.
      real(dp), intent(in), optional :: tau !! Positive resolution multiplier; default is one.
      integer, intent(in) :: resolution_level !! Dyadic resolution level J.
      real(dp), intent(in), optional :: filter_number !! Wavethresh real-filter number; default is 10.
      character(len=*), intent(in), optional :: family !! Wavethresh real-filter family; default is DaubLeAsymm.
      integer, intent(in), optional :: n_iterations !! Binary-product iterations; default is 20.
      logical, intent(in) :: with_covariance !! Whether to calculate upper covariance bands.
      type(density_projection_t) :: projection
      type(support_t) :: support
      type(wt_filter_t) :: filter
      real(dp), allocatable :: phi_values(:)
      real(dp) :: resolution_multiplier
      real(dp) :: selected_filter_number
      real(dp) :: z
      real(dp) :: phi_k
      real(dp) :: phi_l
      character(len=24) :: selected_family
      integer :: iterations
      integer :: n_filter_minus_one
      integer :: n_coefficients
      integer :: observation
      integer :: local_min
      integer :: local_max
      integer :: k
      integer :: l
      integer :: k_index
      integer :: l_index

      if (size(x) == 0) then
         projection%message = "density projection requires at least one observation"
         return
      end if
      if (any(.not. ieee_is_finite(x))) then
         projection%message = "density projection observations must be finite"
         return
      end if
      resolution_multiplier = 1.0_dp
      if (present(tau)) resolution_multiplier = tau
      if (.not. ieee_is_finite(resolution_multiplier) .or. resolution_multiplier <= 0.0_dp) then
         projection%message = "tau must be finite and positive"
         return
      end if
      selected_filter_number = 10.0_dp
      if (present(filter_number)) selected_filter_number = filter_number
      selected_family = "DaubLeAsymm"
      if (present(family)) selected_family = family
      iterations = 20
      if (present(n_iterations)) iterations = n_iterations
      if (iterations < 1) then
         projection%message = "n_iterations must be positive"
         return
      end if

      filter = filter_select(selected_filter_number, trim(selected_family))
      if (.not. filter%ok .or. filter%is_complex) then
         projection%message = "density projection requires a supported real wavelet filter"
         return
      end if
      support = wavelet_support(selected_filter_number, trim(selected_family))
      if (.not. support%ok) then
         projection%message = support%message
         return
      end if
      projection%primary_resolution = resolution_multiplier * 2.0_dp**resolution_level
      if (.not. ieee_is_finite(projection%primary_resolution) .or. projection%primary_resolution <= 0.0_dp) then
         projection%message = "tau*2**resolution_level is not a finite positive resolution"
         return
      end if
      projection%k_min = ceiling(projection%primary_resolution * minval(x) - support%phi_right)
      projection%k_max = floor(projection%primary_resolution * maxval(x) - support%phi_left)
      n_coefficients = projection%k_max - projection%k_min + 1
      if (n_coefficients < 1) then
         projection%message = "density projection produced an empty translation range"
         return
      end if

      n_filter_minus_one = size(filter%low) - 1
      allocate(projection%coefficients(n_coefficients), source=0.0_dp)
      if (with_covariance) then
         allocate(projection%covariance(n_coefficients, n_filter_minus_one), source=0.0_dp)
      end if
      do observation = 1, size(x)
         z = projection%primary_resolution * x(observation)
         local_min = ceiling(z - support%phi_right)
         local_max = floor(z - support%phi_left)
         phi_values = daubechies_lagarias_phi(z, filter%low, iterations)
         do k = local_min, local_max
            k_index = k - projection%k_min + 1
            if (k_index < 1 .or. k_index > n_coefficients) cycle
            if (k - local_min + 1 > size(phi_values)) cycle
            phi_k = sqrt(projection%primary_resolution) * phi_values(k - local_min + 1)
            projection%coefficients(k_index) = projection%coefficients(k_index) + phi_k / real(size(x), dp)
            if (.not. with_covariance) cycle
            do l = k, min(k + n_filter_minus_one - 1, local_max)
               l_index = l - local_min + 1
               if (l_index > size(phi_values)) cycle
               phi_l = sqrt(projection%primary_resolution) * phi_values(l_index)
               projection%covariance(k_index, l - k + 1) = projection%covariance(k_index, l - k + 1) + &
                  phi_k * phi_l / real(size(x), dp)**2
            end do
         end do
      end do

      projection%filter = filter
      projection%tau = resolution_multiplier
      projection%resolution_level = resolution_level
      projection%iterations = iterations
      projection%sample_size = size(x)
      projection%has_covariance = with_covariance
      projection%ok = .true.
      projection%message = "ok"
   end function project_density

   pure function daubechies_lagarias_phi(y, filter, iterations) result(values)
      !! Evaluates all nonzero integer translates of phi at y by binary transition-matrix products.
      real(dp), intent(in) :: y !! Scaling-function coordinate whose fractional part selects the binary product.
      real(dp), intent(in) :: filter(:) !! Normalized low-pass filter with length n+1.
      integer, intent(in) :: iterations !! Number of binary digits and transition matrices to apply.
      real(dp), allocatable :: values(:)
      real(dp), allocatable :: product_matrix(:,:)
      real(dp), allocatable :: next_matrix(:,:)
      real(dp) :: fractional
      integer :: n
      integer :: digit
      integer :: step
      integer :: row
      integer :: column
      integer :: inner
      integer :: filter_index

      n = size(filter) - 1
      if (n < 1 .or. iterations < 1) then
         allocate(values(0))
         return
      end if
      allocate(product_matrix(n, n), source=0.0_dp)
      allocate(next_matrix(n, n), source=0.0_dp)
      do row = 1, n
         product_matrix(row, row) = 1.0_dp
      end do
      fractional = y - floor(y)
      do step = 1, iterations
         fractional = 2.0_dp * fractional
         digit = floor(fractional)
         fractional = fractional - real(digit, dp)
         next_matrix = 0.0_dp
         do column = 1, n
            do inner = 1, n
               if (digit == 0) then
                  filter_index = 2 * inner - column
               else
                  filter_index = 2 * inner - column + 1
               end if
               if (filter_index < 1 .or. filter_index > size(filter)) cycle
               next_matrix(:, column) = next_matrix(:, column) + &
                  product_matrix(:, inner) * sqrt(2.0_dp) * filter(filter_index)
            end do
         end do
         product_matrix = next_matrix
      end do
      allocate(values(n))
      do row = 1, n
         values(n - row + 1) = sum(product_matrix(row, :)) / real(n, dp)
      end do
   end function daubechies_lagarias_phi

   elemental function dclaw(x) result(density)
      real(dp), intent(in) :: x !! Point at which the claw-mixture density is evaluated.
      real(dp) :: density
      real(dp) :: mean_value
      integer :: i
      density = 0.5_dp * normal_pdf(x, 0.0_dp, 1.0_dp)
      do i = 0, 4
         mean_value = real(i, dp) / 2.0_dp - 1.0_dp
         density = density + 0.1_dp * normal_pdf(x, mean_value, 0.1_dp)
      end do
   end function dclaw

   elemental function pclaw(q) result(probability)
      real(dp), intent(in) :: q !! Point at which the claw-mixture cumulative probability is evaluated.
      real(dp) :: probability
      real(dp) :: mean_value
      integer :: i
      probability = 0.5_dp * normal_cdf(q, 0.0_dp, 1.0_dp)
      do i = 0, 4
         mean_value = real(i, dp) / 2.0_dp - 1.0_dp
         probability = probability + 0.1_dp * normal_cdf(q, mean_value, 0.1_dp)
      end do
   end function pclaw

   function rclaw(n, seed) result(values)
      integer, intent(in) :: n !! Number of random claw-mixture draws requested.
      integer, intent(in), optional :: seed !! Optional deterministic seed for the Fortran intrinsic RNG.
      real(dp), allocatable :: values(:)
      real(dp), allocatable :: uniforms(:)
      integer :: i
      if (n < 0) then
         allocate(values(0))
         return
      end if
      if (present(seed)) call seed_rng(seed)
      values = normal_vector(n)
      allocate(uniforms(n))
      call random_number(uniforms)
      do i = 1, n
         if (uniforms(i) <= 0.5_dp) then
            values(i) = values(i) / 10.0_dp + real(int(uniforms(i) * 10.0_dp), dp) / 2.0_dp - 1.0_dp
         end if
      end do
   end function rclaw

   elemental function normal_pdf(x, mean_value, standard_deviation) result(density)
      real(dp), intent(in) :: x !! Normal variate at which to evaluate the density.
      real(dp), intent(in) :: mean_value !! Normal mean parameter.
      real(dp), intent(in) :: standard_deviation !! Positive normal standard deviation.
      real(dp) :: density
      real(dp) :: z
      if (standard_deviation <= 0.0_dp) then
         density = 0.0_dp
         return
      end if
      z = (x - mean_value) / standard_deviation
      density = exp(-0.5_dp * z * z) / (standard_deviation * sqrt(2.0_dp * acos(-1.0_dp)))
   end function normal_pdf

   elemental function normal_cdf(x, mean_value, standard_deviation) result(probability)
      real(dp), intent(in) :: x !! Normal variate at which to evaluate the CDF.
      real(dp), intent(in) :: mean_value !! Normal mean parameter.
      real(dp), intent(in) :: standard_deviation !! Positive normal standard deviation.
      real(dp) :: probability
      if (standard_deviation <= 0.0_dp) then
         probability = merge(1.0_dp, 0.0_dp, x >= mean_value)
      else
         probability = 0.5_dp * erfc(-(x - mean_value) / (standard_deviation * sqrt(2.0_dp)))
      end if
   end function normal_cdf

   function normal_vector(n) result(z)
      integer, intent(in) :: n !! Number of independent standard-normal variates requested.
      real(dp), allocatable :: z(:)
      real(dp) :: u1
      real(dp) :: u2
      integer :: i
      allocate(z(n))
      i = 1
      do while (i <= n)
         call random_number(u1)
         call random_number(u2)
         u1 = max(u1, tiny(1.0_dp))
         z(i) = sqrt(-2.0_dp * log(u1)) * cos(2.0_dp * acos(-1.0_dp) * u2)
         if (i + 1 <= n) z(i + 1) = sqrt(-2.0_dp * log(u1)) * sin(2.0_dp * acos(-1.0_dp) * u2)
         i = i + 2
      end do
   end function normal_vector

   subroutine seed_rng(seed)
      integer, intent(in) :: seed !! User seed expanded deterministically to the compiler RNG seed vector.
      integer, allocatable :: put(:)
      integer :: nseed
      integer :: i
      call random_seed(size=nseed)
      allocate(put(nseed))
      do i = 1, nseed
         put(i) = modulo(abs(seed) + 130363 * i, huge(1) - 1)
         if (put(i) == 0) put(i) = i
      end do
      call random_seed(put=put)
   end subroutine seed_rng

end module wavethresh_density
