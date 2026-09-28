! SPDX-License-Identifier: GPL-3.0-only
! Derived from R package ftsa 6.7 by Rob Hyndman and Han Lin Shang.
module ftsa_statistics
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_quiet_nan, ieee_value
   use r_kinds, only : dp
   use r_descriptive, only : r_mean, r_variance, r_sd
   use r_quantiles, only : r_median, r_quantile_type7
   use r_vectors, only : r_difference
   implicit none
   private

   integer, parameter, public :: ftsa_ok = 0
   integer, parameter, public :: ftsa_invalid_input = 1
   integer, parameter, public :: ftsa_unsupported = 2

   public :: ftsa_extract_time, ftsa_extract_x, ftsa_residuals, ftsa_difference
   public :: ftsa_mean_coordinate, ftsa_median_coordinate, ftsa_variance_coordinate
   public :: ftsa_sd_coordinate, ftsa_quantiles, ftsa_sd_vector, ftsa_sd_columns
   public :: ftsa_centre_mean, ftsa_centre_variance, ftsa_error_metric, ftsa_error_all
   public :: ftsa_bootstrap_mean

contains

   pure subroutine ftsa_extract_time(y, indices, extracted, status)
      real(dp), intent(in) :: y(:, :) !! Functional observations with grid points in rows and times in columns.
      integer, intent(in) :: indices(:) !! One-based time-column indices to retain, in requested output order.
      real(dp), allocatable, intent(out) :: extracted(:, :) !! Selected grid-by-time matrix.
      integer, intent(out) :: status !! Zero on success; nonzero when an index is outside the input time range.
      integer :: j

      status = ftsa_ok
      if (any(indices < 1) .or. any(indices > size(y, 2))) then
         allocate(extracted(0, 0))
         status = ftsa_invalid_input
         return
      end if
      allocate(extracted(size(y, 1), size(indices)))
      do j = 1, size(indices)
         extracted(:, j) = y(:, indices(j))
      end do
   end subroutine ftsa_extract_time

   pure subroutine ftsa_extract_x(x, y, indices, x_out, y_out, status)
      real(dp), intent(in) :: x(:) !! Functional grid corresponding to the rows of `y`.
      real(dp), intent(in) :: y(:, :) !! Functional observations with grid points in rows and times in columns.
      integer, intent(in) :: indices(:) !! One-based grid-row indices to retain, in requested output order.
      real(dp), allocatable, intent(out) :: x_out(:) !! Selected grid values.
      real(dp), allocatable, intent(out) :: y_out(:, :) !! Selected functional rows across all times.
      integer, intent(out) :: status !! Zero on success; nonzero for nonconforming dimensions or invalid indices.
      integer :: i

      status = ftsa_ok
      if (size(x) /= size(y, 1) .or. any(indices < 1) .or. any(indices > size(x))) then
         allocate(x_out(0), y_out(0, 0))
         status = ftsa_invalid_input
         return
      end if
      allocate(x_out(size(indices)), y_out(size(indices), size(y, 2)))
      do i = 1, size(indices)
         x_out(i) = x(indices(i))
         y_out(i, :) = y(indices(i), :)
      end do
   end subroutine ftsa_extract_x

   pure subroutine ftsa_residuals(y, fitted, residuals, status)
      real(dp), intent(in) :: y(:, :) !! Observed functional data matrix.
      real(dp), intent(in) :: fitted(:, :) !! Fitted functional data matrix with the same shape as `y`.
      real(dp), allocatable, intent(out) :: residuals(:, :) !! Elementwise observed-minus-fitted residual matrix.
      integer, intent(out) :: status !! Zero on success; nonzero when the matrix shapes differ.

      status = ftsa_ok
      if (any(shape(y) /= shape(fitted))) then
         allocate(residuals(0, 0))
         status = ftsa_invalid_input
         return
      end if
      allocate(residuals(size(y, 1), size(y, 2)))
      residuals = y - fitted
   end subroutine ftsa_residuals

   pure subroutine ftsa_difference(y, lag, differences, differenced, status)
      real(dp), intent(in) :: y(:, :) !! Functional observations with time running along the second dimension.
      integer, intent(in), optional :: lag !! Positive lag between differenced observations; default is one.
      integer, intent(in), optional :: differences !! Number of repeated differencing operations; default is one.
      real(dp), allocatable, intent(out) :: differenced(:, :) !! Differenced functional matrix with shortened time dimension.
      integer, intent(out) :: status !! Zero on success; nonzero for invalid lag/difference settings.
      real(dp), allocatable :: temp(:, :)
      integer :: core_status

      call r_difference(transpose(y), temp, lag, differences, core_status)
      if (core_status /= 0) then
         allocate(differenced(0, 0))
         status = ftsa_invalid_input
         return
      end if
      allocate(differenced(size(temp, 2), size(temp, 1)))
      differenced = transpose(temp)
      status = ftsa_ok
   end subroutine ftsa_difference

   pure subroutine ftsa_mean_coordinate(y, mean_curve, na_rm)
      real(dp), intent(in) :: y(:, :) !! Functional observations with replicates in columns.
      real(dp), intent(out) :: mean_curve(:) !! Coordinatewise mean; size must equal the number of grid rows.
      logical, intent(in), optional :: na_rm !! Omit IEEE NaNs when true, matching `mean.fts(..., na.rm=TRUE)`.
      integer :: i

      do i = 1, min(size(mean_curve), size(y, 1))
         mean_curve(i) = r_mean(y(i, :), na_rm)
      end do
   end subroutine ftsa_mean_coordinate

   pure subroutine ftsa_median_coordinate(y, median_curve, na_rm)
      real(dp), intent(in) :: y(:, :) !! Functional observations with replicates in columns.
      real(dp), intent(out) :: median_curve(:) !! Coordinatewise type-7 median curve.
      logical, intent(in), optional :: na_rm !! Omit IEEE NaNs when true.
      integer :: i

      do i = 1, min(size(median_curve), size(y, 1))
         median_curve(i) = r_median(y(i, :), na_rm)
      end do
   end subroutine ftsa_median_coordinate

   pure subroutine ftsa_variance_coordinate(y, variance_curve, na_rm)
      real(dp), intent(in) :: y(:, :) !! Functional observations with replicates in columns.
      real(dp), intent(out) :: variance_curve(:) !! Coordinatewise population variance using divisor `n`.
      logical, intent(in), optional :: na_rm !! Omit IEEE NaNs when true.
      integer :: i

      do i = 1, min(size(variance_curve), size(y, 1))
         variance_curve(i) = r_variance(y(i, :), na_rm=na_rm, ddof=0)
      end do
   end subroutine ftsa_variance_coordinate

   pure subroutine ftsa_sd_coordinate(y, sd_curve, status)
      real(dp), intent(in) :: y(:, :) !! Finite functional observations with replicates in columns.
      real(dp), intent(out) :: sd_curve(:) !! Coordinatewise scale used by upstream `sd.fts` coordinate method.
      integer, intent(out) :: status !! Zero on success; nonzero when values are non-finite or dimensions are invalid.
      integer :: i, n

      status = ftsa_ok
      n = size(y, 2)
      if (n < 2 .or. size(sd_curve) /= size(y, 1) .or. any(.not. ieee_is_finite(y))) then
         sd_curve = ieee_value(0.0_dp, ieee_quiet_nan)
         status = ftsa_invalid_input
         return
      end if
      do i = 1, size(y, 1)
         sd_curve(i) = real(n - 1, dp) * r_sd(y(i, :), ddof=1) / real(n, dp)
      end do
   end subroutine ftsa_sd_coordinate

   pure subroutine ftsa_quantiles(y, probabilities, values, na_rm, status)
      real(dp), intent(in) :: y(:, :) !! Functional observations with replicates in columns.
      real(dp), intent(in) :: probabilities(:) !! Requested probabilities in the closed interval `[0,1]`.
      real(dp), allocatable, intent(out) :: values(:, :) !! Grid-by-probability matrix of R type-7 quantiles.
      logical, intent(in), optional :: na_rm !! Omit IEEE NaNs when true.
      integer, intent(out) :: status !! Zero on success; nonzero when any probability is outside `[0,1]`.
      integer :: i, j

      if (any(probabilities < 0.0_dp) .or. any(probabilities > 1.0_dp)) then
         allocate(values(0, 0))
         status = ftsa_invalid_input
         return
      end if
      allocate(values(size(y, 1), size(probabilities)))
      do j = 1, size(probabilities)
         do i = 1, size(y, 1)
            values(i, j) = r_quantile_type7(y(i, :), probabilities(j), na_rm)
         end do
      end do
      status = ftsa_ok
   end subroutine ftsa_quantiles

   pure real(dp) function ftsa_sd_vector(x, na_rm) result(value)
      real(dp), intent(in) :: x(:) !! Numeric vector whose sample standard deviation is requested.
      logical, intent(in), optional :: na_rm !! Omit IEEE NaNs when true.

      value = r_sd(x, na_rm=na_rm, ddof=1)
   end function ftsa_sd_vector

   pure subroutine ftsa_sd_columns(x, values, na_rm)
      real(dp), intent(in) :: x(:, :) !! Numeric matrix whose columns are treated independently.
      real(dp), intent(out) :: values(:) !! Sample standard deviations, one per input column.
      logical, intent(in), optional :: na_rm !! Omit IEEE NaNs within each column when true.
      integer :: j

      do j = 1, min(size(values), size(x, 2))
         values(j) = r_sd(x(:, j), na_rm=na_rm, ddof=1)
      end do
   end subroutine ftsa_sd_columns

   pure subroutine ftsa_centre_mean(x, center)
      real(dp), intent(in) :: x(:, :) !! Functional matrix with grid points in rows and replicates in columns.
      real(dp), intent(out) :: center(:) !! Coordinatewise arithmetic center.

      call ftsa_mean_coordinate(x, center, na_rm=.false.)
   end subroutine ftsa_centre_mean

   pure subroutine ftsa_centre_variance(x, center)
      real(dp), intent(in) :: x(:, :) !! Functional matrix with grid points in rows and replicates in columns.
      real(dp), intent(out) :: center(:) !! Coordinatewise population variance center.

      call ftsa_variance_coordinate(x, center, na_rm=.false.)
   end subroutine ftsa_centre_variance

   pure subroutine ftsa_error_metric(forecast, true_values, method, value, forecast_bench, insample_true, status)
      real(dp), intent(in) :: forecast(:) !! Forecast values to evaluate.
      real(dp), intent(in) :: true_values(:) !! Realized out-of-sample values corresponding to `forecast`.
      character(len=*), intent(in) :: method !! Error measure name such as `mae`, `mse`, `mase`, or `rmsse`.
      real(dp), intent(out) :: value !! Requested scalar error measure, including upstream rounding.
      real(dp), intent(in), optional :: forecast_bench(:) !! Benchmark forecast required by relative-error measures.
      real(dp), intent(in), optional :: insample_true(:) !! In-sample truth required by scaled-error measures.
      integer, intent(out) :: status !! Zero on success; nonzero for missing required inputs or shape errors.
      real(dp), allocatable :: error(:), ratio(:), scale_errors(:)
      real(dp) :: denominator
      character(len=:), allocatable :: key
      integer :: n

      status = ftsa_ok
      value = ieee_value(0.0_dp, ieee_quiet_nan)
      if (size(forecast) /= size(true_values) .or. size(forecast) == 0) then
         status = ftsa_invalid_input
         return
      end if
      n = size(forecast)
      allocate(error(n))
      error = true_values - forecast
      key = lower_ascii(trim(method))
      select case (key)
      case ('me')
         value = round_places(sum(error) / real(n, dp), 6)
      case ('mpe')
         value = round_places(sum(100.0_dp * error / true_values) / real(n, dp), 6)
      case ('mae')
         value = round_places(sum(abs(error)) / real(n, dp), 6)
      case ('mse')
         value = round_places(sum(error**2) / real(n, dp), 6)
      case ('sse')
         value = round_places(sum(error**2), 6)
      case ('rmse')
         value = round_places(sqrt(sum(error**2) / real(n, dp)), 6)
      case ('mdae')
         value = round_places(r_median(abs(error)), 4)
      case ('mdse')
         value = round_places(r_median(error**2), 4)
      case ('mape')
         value = round_places(sum(100.0_dp * abs(error / true_values)) / real(n, dp), 6)
      case ('mdape')
         value = round_places(r_median(100.0_dp * abs(error / true_values)), 4)
      case ('smape')
         value = round_places(sum(200.0_dp * abs(error) / (true_values + forecast)) / real(n, dp), 6)
      case ('smdape')
         value = round_places(r_median(200.0_dp * abs(error) / (true_values + forecast)), 6)
      case ('rmspe')
         value = round_places(sqrt(sum((100.0_dp * error / true_values)**2) / real(n, dp)), 6)
      case ('rmdspe')
         value = round_places(sqrt(r_median((100.0_dp * error / true_values)**2)), 6)
      case ('mrae', 'mdrae', 'gmrae', 'relmae', 'relmse')
         if (.not. present(forecast_bench)) then
            status = ftsa_invalid_input
            return
         end if
         if (size(forecast_bench) /= n) then
            status = ftsa_invalid_input
            return
         end if
         allocate(ratio(n))
         ratio = error / (true_values - forecast_bench)
         select case (key)
         case ('mrae')
            value = round_places(sum(abs(ratio)) / real(n, dp), 6)
         case ('mdrae')
            value = round_places(r_median(abs(ratio)), 6)
         case ('gmrae')
            if (any(abs(ratio) <= tiny(1.0_dp))) then
               value = 0.0_dp
            else
               value = round_places(exp(sum(log(abs(ratio))) / real(n, dp)), 6)
            end if
         case ('relmae')
            denominator = sum(abs(true_values - forecast_bench)) / real(n, dp)
            value = round_places((sum(abs(error)) / real(n, dp)) / denominator, 6)
         case ('relmse')
            denominator = sum((true_values - forecast_bench)**2) / real(n, dp)
            value = round_places((sum(error**2) / real(n, dp)) / denominator, 6)
         end select
      case ('mase', 'mdase', 'rmsse')
         if (.not. present(insample_true)) then
            status = ftsa_invalid_input
            return
         end if
         if (size(insample_true) < 2) then
            status = ftsa_invalid_input
            return
         end if
         allocate(scale_errors(size(insample_true) - 1))
         scale_errors = abs(insample_true(2:) - insample_true(:size(insample_true) - 1))
         denominator = sum(scale_errors) / real(size(scale_errors), dp)
         select case (key)
         case ('mase')
            value = round_places(sum(abs(error / denominator)) / real(n, dp), 6)
         case ('mdase')
            value = round_places(r_median(abs(error / denominator)), 6)
         case ('rmsse')
            value = round_places(sqrt(sum((error / denominator)**2) / real(n, dp)), 6)
         end select
      case default
         status = ftsa_unsupported
      end select
   end subroutine ftsa_error_metric

   pure subroutine ftsa_error_all(forecast, true_values, forecast_bench, insample_true, values, status)
      real(dp), intent(in) :: forecast(:) !! Forecast vector evaluated by every supported upstream error measure.
      real(dp), intent(in) :: true_values(:) !! Realized out-of-sample values.
      real(dp), intent(in) :: forecast_bench(:) !! Benchmark forecast for relative error measures.
      real(dp), intent(in) :: insample_true(:) !! In-sample realized values for scaled error measures.
      real(dp), intent(out) :: values(22) !! Measures in upstream `giveall=TRUE` order from ME through RMSSE.
      integer, intent(out) :: status !! Zero when all measures are computed successfully.
      character(len=8), parameter :: names(22) = [character(len=8) :: &
         'me', 'mpe', 'mae', 'mse', 'sse', 'rmse', 'mdae', 'mdse', 'mape', 'mdape', 'smape', &
         'smdape', 'rmspe', 'rmdspe', 'mrae', 'mdrae', 'gmrae', 'relmae', 'relmse', 'mase', 'mdase', 'rmsse']
      integer :: i, local_status

      status = ftsa_ok
      do i = 1, size(names)
         call ftsa_error_metric(forecast, true_values, trim(names(i)), values(i), forecast_bench, &
            insample_true, local_status)
         if (local_status /= ftsa_ok) status = local_status
      end do
   end subroutine ftsa_error_all

   subroutine ftsa_bootstrap_mean(data, nb, alpha, estimate, max_dist, rep_dist, resamples, center, status, seed, center_boot)
      real(dp), intent(in) :: data(:, :) !! Functional data with grid points in rows and observed curves in columns.
      integer, intent(in) :: nb !! Number of nonparametric bootstrap resamples; must be positive.
      real(dp), intent(in) :: alpha !! Upper-tail probability excluded from the bootstrap radius, in `[0,1)`.
      real(dp), allocatable, intent(out) :: estimate(:) !! Coordinatewise sample mean of the original curves.
      real(dp), intent(out) :: max_dist !! Empirical `(1-alpha)` radius of Euclidean bootstrap distances.
      real(dp), allocatable, intent(out) :: rep_dist(:) !! Euclidean distance for every bootstrap replicate.
      real(dp), allocatable, intent(out) :: resamples(:, :) !! Bootstrap mean curves, one per column.
      real(dp), allocatable, intent(out) :: center(:) !! Original mean or bootstrap-mean center used for distances.
      integer, intent(out) :: status !! Zero on success; nonzero for invalid dimensions or probability settings.
      integer, intent(in), optional :: seed !! Optional deterministic seed for the processor's Fortran RNG.
      logical, intent(in), optional :: center_boot !! Use the mean of bootstrap estimates as distance center when true.
      real(dp), allocatable :: sorted(:)
      real(dp) :: u
      integer :: b, draw, i, keep, n
      logical :: use_boot_center

      status = ftsa_ok
      max_dist = ieee_value(0.0_dp, ieee_quiet_nan)
      n = size(data, 2)
      if (nb < 1 .or. n < 1 .or. size(data, 1) < 1 .or. alpha < 0.0_dp .or. alpha >= 1.0_dp) then
         allocate(estimate(0), rep_dist(0), resamples(0, 0), center(0))
         status = ftsa_invalid_input
         return
      end if
      if (present(seed)) call set_rng_seed(seed)
      allocate(estimate(size(data, 1)), rep_dist(nb), resamples(size(data, 1), nb), center(size(data, 1)))
      call ftsa_mean_coordinate(data, estimate, na_rm=.false.)
      do b = 1, nb
         resamples(:, b) = 0.0_dp
         do i = 1, n
            call random_number(u)
            draw = min(n, 1 + int(u * real(n, dp)))
            resamples(:, b) = resamples(:, b) + data(:, draw)
         end do
         resamples(:, b) = resamples(:, b) / real(n, dp)
      end do
      use_boot_center = .false.
      if (present(center_boot)) use_boot_center = center_boot
      if (use_boot_center) then
         do i = 1, size(data, 1)
            center(i) = sum(resamples(i, :)) / real(nb, dp)
         end do
      else
         center = estimate
      end if
      do b = 1, nb
         rep_dist(b) = sqrt(sum((resamples(:, b) - center)**2))
      end do
      allocate(sorted(nb))
      sorted = rep_dist
      call sort_in_place(sorted)
      keep = max(1, floor((1.0_dp - alpha) * real(nb, dp)))
      max_dist = sorted(min(nb, keep))
   end subroutine ftsa_bootstrap_mean

   pure real(dp) function round_places(x, digits) result(value)
      real(dp), intent(in) :: x !! Finite or IEEE special value to round.
      integer, intent(in) :: digits !! Number of decimal places retained.
      real(dp) :: scale

      if (.not. ieee_is_finite(x)) then
         value = x
         return
      end if
      scale = 10.0_dp**digits
      value = anint(x * scale) / scale
   end function round_places

   pure function lower_ascii(text) result(lowered)
      character(len=*), intent(in) :: text !! ASCII text to convert to lower case.
      character(len=len(text)) :: lowered
      integer :: code, i

      lowered = text
      do i = 1, len(text)
         code = iachar(text(i:i))
         if (code >= iachar('A') .and. code <= iachar('Z')) lowered(i:i) = achar(code + 32)
      end do
   end function lower_ascii

   pure subroutine sort_in_place(x)
      real(dp), intent(inout) :: x(:) !! Values sorted in ascending order in place.
      real(dp) :: key
      integer :: i, j

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
   end subroutine sort_in_place

   subroutine set_rng_seed(seed)
      integer, intent(in) :: seed !! User seed expanded deterministically to the processor RNG seed vector.
      integer, allocatable :: values(:)
      integer :: i, n

      call random_seed(size=n)
      allocate(values(n))
      do i = 1, n
         values(i) = modulo(seed + 104729 * i, huge(1) - 1)
         if (values(i) <= 0) values(i) = i
      end do
      call random_seed(put=values)
   end subroutine set_rng_seed

end module ftsa_statistics
