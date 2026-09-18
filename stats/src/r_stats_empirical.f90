! SPDX-License-Identifier: MIT
! SPDX-FileComment: Empirical distributions, histograms, and kernel-density estimates.
module r_stats_empirical
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_is_nan
   use r_descriptive, only: r_sd
   use r_kinds, only: dp
   use r_optional, only: optval
   use r_ordering, only: r_sort_values_in_place
   use r_quantiles, only: r_quantile_type7
   use r_stats_types, only: density_result_t, ecdf_t, histogram_t
   implicit none
   private

   integer, parameter, public :: density_kernel_gaussian = 1
   integer, parameter, public :: density_kernel_rectangular = 2
   integer, parameter, public :: density_kernel_triangular = 3
   integer, parameter, public :: density_kernel_epanechnikov = 4
   real(dp), parameter :: pi = acos(-1.0_dp)

   public :: bw_nrd, bw_nrd0, density, ecdf, histogram, predict_ecdf

   interface histogram
      module procedure histogram_equal, histogram_explicit
   end interface histogram

contains

   pure function ecdf(x, na_rm) result(model)
      !! Constructs a compact empirical cumulative distribution function.
      real(dp), intent(in) :: x(:) !! Sample values.
      logical, intent(in), optional :: na_rm !! Remove NaNs when true; defaults to false.
      type(ecdf_t) :: model !! Sorted ECDF knots, probabilities, and sample count.
      real(dp), allocatable :: sample(:), unique_values(:), cumulative(:)
      integer :: i, number_unique

      sample = cleaned_sample(x, optval(na_rm, .false.), "ecdf")
      model%n = size(sample)
      if (model%n == 0) then
         allocate (model%values(0), model%probabilities(0))
         return
      end if
      call r_sort_values_in_place(sample)
      allocate (unique_values(model%n), cumulative(model%n))
      number_unique = 1
      unique_values(1) = sample(1)
      do i = 2, model%n
         if (sample(i) /= unique_values(number_unique)) then
            cumulative(number_unique) = real(i - 1, dp)/real(model%n, dp)
            number_unique = number_unique + 1
            unique_values(number_unique) = sample(i)
         end if
      end do
      cumulative(number_unique) = 1.0_dp
      model%values = unique_values(:number_unique)
      model%probabilities = cumulative(:number_unique)
   end function ecdf

   pure function predict_ecdf(model, query) result(probabilities)
      !! Evaluates an empirical distribution at arbitrary query values.
      type(ecdf_t), intent(in) :: model !! Previously constructed empirical distribution.
      real(dp), intent(in) :: query(:) !! Query values.
      real(dp), allocatable :: probabilities(:) !! Fractions of observations not exceeding each query.
      integer :: high, i, low, middle

      allocate (probabilities(size(query)))
      do i = 1, size(query)
         if (ieee_is_nan(query(i))) then
            probabilities(i) = query(i)
            cycle
         end if
         low = 1
         high = size(model%values)
         do while (low <= high)
            middle = (low + high)/2
            if (model%values(middle) <= query(i)) then
               low = middle + 1
            else
               high = middle - 1
            end if
         end do
         if (high == 0) then
            probabilities(i) = 0.0_dp
         else
            probabilities(i) = model%probabilities(high)
         end if
      end do
   end function predict_ecdf

   pure function histogram_explicit(x, breaks, right, include_lowest) result(value)
      !! Counts observations in bins defined by explicit increasing break points.
      real(dp), intent(in) :: x(:) !! Finite observations.
      real(dp), intent(in) :: breaks(:) !! Strictly increasing bin boundaries.
      logical, intent(in), optional :: right !! Use right-closed bins; defaults to true.
      logical, intent(in), optional :: include_lowest !! Include the otherwise open endpoint.
      type(histogram_t) :: value !! Histogram counts, densities, breaks, and midpoints.
      real(dp), allocatable :: widths(:)
      logical :: close_right, include_endpoint
      integer :: bin, i, number_bins

      number_bins = size(breaks) - 1
      if (number_bins < 1) error stop "histogram: at least two breaks are required"
      if (.not. all(ieee_is_finite(x)) .or. .not. all(ieee_is_finite(breaks))) then
         error stop "histogram: values and breaks must be finite"
      end if
      widths = breaks(2:) - breaks(:number_bins)
      if (any(widths <= 0.0_dp)) error stop "histogram: breaks must be strictly increasing"
      if (any(x < breaks(1)) .or. any(x > breaks(number_bins + 1))) then
         error stop "histogram: some observations lie outside the breaks"
      end if
      close_right = optval(right, .true.)
      include_endpoint = optval(include_lowest, .true.)
      value%breaks = breaks
      value%mids = 0.5_dp*(breaks(:number_bins) + breaks(2:))
      allocate (value%counts(number_bins), source=0)
      do i = 1, size(x)
         bin = locate_bin(x(i), breaks, close_right, include_endpoint)
         if (bin == 0) error stop "histogram: an endpoint observation was excluded"
         value%counts(bin) = value%counts(bin) + 1
      end do
      allocate (value%density(number_bins), source=0.0_dp)
      if (size(x) > 0) value%density = real(value%counts, dp)/(real(size(x), dp)*widths)
      value%equidistant = maxval(abs(widths - widths(1))) <= &
                          100.0_dp*epsilon(1.0_dp)*max(1.0_dp, abs(widths(1)))
   end function histogram_explicit

   pure function histogram_equal(x, number_bins, right, include_lowest) result(value)
      !! Constructs equal-width bins spanning the finite sample range.
      real(dp), intent(in) :: x(:) !! Nonempty finite observations.
      integer, intent(in) :: number_bins !! Positive number of equal-width bins.
      logical, intent(in), optional :: right !! Use right-closed bins; defaults to true.
      logical, intent(in), optional :: include_lowest !! Include the otherwise open endpoint.
      type(histogram_t) :: value !! Histogram counts, densities, breaks, and midpoints.
      real(dp), allocatable :: breaks(:)
      real(dp) :: lower, upper, width
      integer :: i

      if (size(x) == 0 .or. .not. all(ieee_is_finite(x))) then
         error stop "histogram: a nonempty finite sample is required"
      end if
      if (number_bins < 1) error stop "histogram: number of bins must be positive"
      lower = minval(x)
      upper = maxval(x)
      if (upper == lower) then
         lower = lower - 0.5_dp
         upper = upper + 0.5_dp
      end if
      width = (upper - lower)/real(number_bins, dp)
      allocate (breaks(number_bins + 1))
      do i = 1, number_bins + 1
         breaks(i) = lower + real(i - 1, dp)*width
      end do
      breaks(number_bins + 1) = upper
      value = histogram_explicit(x, breaks, right, include_lowest)
   end function histogram_equal

   pure function bw_nrd0(x, na_rm) result(bandwidth)
      !! Computes R's robust default `bw.nrd0` bandwidth selector.
      real(dp), intent(in) :: x(:) !! Sample values.
      logical, intent(in), optional :: na_rm !! Remove NaNs when true; defaults to false.
      real(dp) :: bandwidth !! Positive normal-reference bandwidth.
      real(dp), allocatable :: sample(:)
      real(dp) :: high_scale, low_scale

      sample = cleaned_sample(x, optval(na_rm, .false.), "bw_nrd0")
      if (size(sample) < 2) error stop "bw_nrd0: at least two observations are required"
      high_scale = r_sd(sample)
      low_scale = min(high_scale, &
                      (r_quantile_type7(sample, 0.75_dp) - r_quantile_type7(sample, 0.25_dp))/1.34_dp)
      if (low_scale == 0.0_dp) low_scale = high_scale
      if (low_scale == 0.0_dp) low_scale = abs(sample(1))
      if (low_scale == 0.0_dp) low_scale = 1.0_dp
      bandwidth = 0.9_dp*low_scale*real(size(sample), dp)**(-0.2_dp)
   end function bw_nrd0

   pure function bw_nrd(x, na_rm) result(bandwidth)
      !! Computes R's `bw.nrd` normal-reference bandwidth selector.
      real(dp), intent(in) :: x(:) !! Sample values.
      logical, intent(in), optional :: na_rm !! Remove NaNs when true; defaults to false.
      real(dp) :: bandwidth !! Normal-reference bandwidth.
      real(dp), allocatable :: sample(:)
      real(dp) :: scale

      sample = cleaned_sample(x, optval(na_rm, .false.), "bw_nrd")
      if (size(sample) < 2) error stop "bw_nrd: at least two observations are required"
      scale = min(r_sd(sample), &
                  (r_quantile_type7(sample, 0.75_dp) - r_quantile_type7(sample, 0.25_dp))/1.34_dp)
      bandwidth = 1.06_dp*scale*real(size(sample), dp)**(-0.2_dp)
   end function bw_nrd

   pure function density(x, bandwidth, weights, adjust, number_points, from, to, kernel) result(value)
      !! Evaluates a direct univariate kernel-density estimate on a regular grid.
      real(dp), intent(in) :: x(:) !! Nonempty finite sample.
      real(dp), intent(in), optional :: bandwidth !! Positive base bandwidth; defaults to `bw_nrd0`.
      real(dp), intent(in), optional :: weights(:) !! Nonnegative observation weights.
      real(dp), intent(in), optional :: adjust !! Positive bandwidth multiplier; defaults to one.
      integer, intent(in), optional :: number_points !! Grid size; defaults to 512.
      real(dp), intent(in), optional :: from !! Lower grid endpoint; defaults to `min(x)-3*bw`.
      real(dp), intent(in), optional :: to !! Upper grid endpoint; defaults to `max(x)+3*bw`.
      integer, intent(in), optional :: kernel !! Kernel identifier; defaults to Gaussian.
      type(density_result_t) :: value !! Grid, density values, and fitting metadata.
      real(dp), allocatable :: actual_weights(:)
      real(dp) :: lower, multiplier, total_weight, upper, u
      integer :: i, j, number_grid, selected_kernel

      if (size(x) < 2 .or. .not. all(ieee_is_finite(x))) then
         error stop "density: at least two finite observations are required"
      end if
      if (present(weights)) then
         if (size(weights) /= size(x) .or. any(weights < 0.0_dp) .or. &
             .not. all(ieee_is_finite(weights))) error stop "density: invalid weights"
         actual_weights = weights
      else
         allocate (actual_weights(size(x)), source=1.0_dp)
      end if
      total_weight = sum(actual_weights)
      if (total_weight <= 0.0_dp) error stop "density: weights must have positive sum"
      multiplier = optval(adjust, 1.0_dp)
      if (multiplier <= 0.0_dp) error stop "density: adjust must be positive"
      if (present(bandwidth)) then
         value%bandwidth = bandwidth*multiplier
      else
         value%bandwidth = bw_nrd0(x)*multiplier
      end if
      if (value%bandwidth <= 0.0_dp .or. .not. ieee_is_finite(value%bandwidth)) then
         error stop "density: bandwidth must be positive and finite"
      end if
      number_grid = optval(number_points, 512)
      if (number_grid < 2) error stop "density: at least two grid points are required"
      lower = minval(x) - 3.0_dp*value%bandwidth
      upper = maxval(x) + 3.0_dp*value%bandwidth
      if (present(from)) lower = from
      if (present(to)) upper = to
      if (.not. ieee_is_finite(lower) .or. .not. ieee_is_finite(upper) .or. upper <= lower) then
         error stop "density: invalid grid endpoints"
      end if
      selected_kernel = optval(kernel, density_kernel_gaussian)
      if (selected_kernel < density_kernel_gaussian .or. &
          selected_kernel > density_kernel_epanechnikov) error stop "density: unsupported kernel"
      allocate (value%x(number_grid), value%y(number_grid), source=0.0_dp)
      do i = 1, number_grid
         value%x(i) = lower + real(i - 1, dp)*(upper - lower)/real(number_grid - 1, dp)
         do j = 1, size(x)
            u = (value%x(i) - x(j))/value%bandwidth
            value%y(i) = value%y(i) + actual_weights(j)*kernel_value(u, selected_kernel)
         end do
      end do
      value%y = value%y/(total_weight*value%bandwidth)
      value%n = size(x)
      value%kernel = selected_kernel
   end function density

   pure function cleaned_sample(x, remove_nan, caller) result(sample)
      !! Validates a sample and optionally removes its NaN values.
      real(dp), intent(in) :: x(:) !! Sample to validate.
      logical, intent(in) :: remove_nan !! Whether NaNs should be removed.
      character(len=*), intent(in) :: caller !! Procedure name used in diagnostics.
      real(dp), allocatable :: sample(:) !! Finite retained observations.

      if (remove_nan) then
         sample = pack(x, .not. ieee_is_nan(x))
      else
         if (any(ieee_is_nan(x))) error stop caller//": sample contains NaN"
         sample = x
      end if
      if (.not. all(ieee_is_finite(sample))) error stop caller//": sample contains infinity"
   end function cleaned_sample

   pure function locate_bin(x, breaks, close_right, include_endpoint) result(bin)
      !! Locates a finite observation in explicitly bounded histogram bins.
      real(dp), intent(in) :: x !! Observation to locate.
      real(dp), intent(in) :: breaks(:) !! Strictly increasing bin boundaries.
      logical, intent(in) :: close_right !! Whether bins are right closed.
      logical, intent(in) :: include_endpoint !! Whether to include the opposite extreme endpoint.
      integer :: bin !! One-based bin index, or zero when excluded.
      integer :: i, number_bins

      bin = 0
      number_bins = size(breaks) - 1
      if (close_right) then
         if (include_endpoint .and. x == breaks(1)) then
            bin = 1
            return
         end if
         do i = 1, number_bins
            if (x > breaks(i) .and. x <= breaks(i + 1)) then
               bin = i
               return
            end if
         end do
      else
         if (include_endpoint .and. x == breaks(number_bins + 1)) then
            bin = number_bins
            return
         end if
         do i = 1, number_bins
            if (x >= breaks(i) .and. x < breaks(i + 1)) then
               bin = i
               return
            end if
         end do
      end if
   end function locate_bin

   pure elemental function kernel_value(x, kernel) result(value)
      !! Evaluates a unit-standard-deviation density kernel.
      real(dp), intent(in) :: x !! Standardized displacement.
      integer, intent(in) :: kernel !! Supported kernel identifier.
      real(dp) :: value !! Kernel density.
      real(dp) :: radius

      select case (kernel)
      case (density_kernel_gaussian)
         value = exp(-0.5_dp*x**2)/sqrt(2.0_dp*pi)
      case (density_kernel_rectangular)
         radius = sqrt(3.0_dp)
         value = merge(0.5_dp/radius, 0.0_dp, abs(x) <= radius)
      case (density_kernel_triangular)
         radius = sqrt(6.0_dp)
         value = merge((1.0_dp - abs(x)/radius)/radius, 0.0_dp, abs(x) <= radius)
      case (density_kernel_epanechnikov)
         radius = sqrt(5.0_dp)
         value = merge(0.75_dp*(1.0_dp - (x/radius)**2)/radius, 0.0_dp, abs(x) <= radius)
      case default
         value = 0.0_dp
      end select
   end function kernel_value

end module r_stats_empirical
