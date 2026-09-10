module bivkld_kernel
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_positive_inf, ieee_quiet_nan, ieee_value
   use bivkld_kinds, only : dp
   implicit none
   private

   real(dp), parameter :: pi = acos(-1.0_dp)

   integer, parameter, public :: BIVKLD_SUCCESS = 0
   integer, parameter, public :: BIVKLD_INVALID_INPUT = 1
   integer, parameter, public :: BIVKLD_INVALID_BANDWIDTH = 2
   integer, parameter, public :: BIVKLD_NUMERICAL_FAILURE = 3
   integer, parameter, public :: BIVKLD_INVALID_OPTION = 4

   type, public :: biv_sample
      real(dp), allocatable :: values(:, :)
   end type biv_sample

   type, public :: bivkld_estimate
      real(dp) :: estimate = 0.0_dp
      real(dp) :: bandwidth_x(2, 2) = 0.0_dp
      real(dp) :: bandwidth_y(2, 2) = 0.0_dp
      real(dp), allocatable :: axis1(:)
      real(dp), allocatable :: axis2(:)
      real(dp), allocatable :: density_x(:)
      real(dp), allocatable :: density_y(:)
      real(dp) :: cell_area = 0.0_dp
      real(dp) :: range(2, 2) = 0.0_dp
      character(len=8) :: standardize = 'none'
      character(len=8) :: bandwidth_selector = 'scv'
   end type bivkld_estimate

   public :: biv_kld
   public :: biv_kld_matrix
   public :: hns_bandwidth
   public :: hscv_bandwidth
   public :: kde_density

contains

   pure subroutine biv_kld(x, y, estimate, Hx, Hy, bandwidth, grid_size, limits, standardize, details, info)
      real(dp), intent(in) :: x(:, :) !! First bivariate sample, with observations in rows and exactly two columns.
      real(dp), intent(in) :: y(:, :) !! Second bivariate sample, with observations in rows and exactly two columns.
      real(dp), intent(out) :: estimate !! Estimated directed divergence D(x || y), or NaN when validation fails.
      real(dp), intent(in), optional :: Hx(2, 2) !! Optional symmetric positive-definite bandwidth matrix for x.
      real(dp), intent(in), optional :: Hy(2, 2) !! Optional symmetric positive-definite bandwidth matrix for y.
      character(len=*), intent(in), optional :: bandwidth !! Selector name: 'scv' (default) or 'normal'.
      integer, intent(in), optional :: grid_size(2) !! Grid counts by coordinate; each must be at least 10.
      real(dp), intent(in), optional :: limits(2, 2) !! Integration limits, one [minimum, maximum] row per coordinate.
      character(len=*), intent(in), optional :: standardize !! Standardization mode: 'none', 'pooled', or 'separate'.
      type(bivkld_estimate), intent(out), optional :: details !! Optional bandwidth, grid, density, and settings output.
      integer, intent(out), optional :: info !! Status code; zero succeeds and positive values report errors.

      real(dp), allocatable :: xs(:, :), ys(:, :), axis1(:), axis2(:), px(:), py(:)
      real(dp) :: bandwidth_x(2, 2), bandwidth_y(2, 2), cell_area, range_used(2, 2)
      integer :: ngrid(2), status
      character(len=8) :: selector, std_mode

      estimate = ieee_value(0.0_dp, ieee_quiet_nan)
      call assign_info(info, BIVKLD_SUCCESS)
      if (.not. valid_sample(x) .or. .not. valid_sample(y)) then
         call assign_info(info, BIVKLD_INVALID_INPUT)
         return
      end if

      call resolve_options(bandwidth, standardize, grid_size, selector, std_mode, ngrid, status)
      if (status /= BIVKLD_SUCCESS) then
         call assign_info(info, status)
         return
      end if

      call standardize_pair(x, y, std_mode, xs, ys, status)
      if (status /= BIVKLD_SUCCESS) then
         call assign_info(info, status)
         return
      end if

      if (present(Hx)) then
         bandwidth_x = 0.5_dp*(Hx + transpose(Hx))
         if (.not. is_spd_2x2(Hx)) then
            call assign_info(info, BIVKLD_INVALID_BANDWIDTH)
            return
         end if
      else
         call select_bandwidth(xs, selector, bandwidth_x, status)
         if (status /= BIVKLD_SUCCESS) then
            call assign_info(info, status)
            return
         end if
      end if

      if (present(Hy)) then
         bandwidth_y = 0.5_dp*(Hy + transpose(Hy))
         if (.not. is_spd_2x2(Hy)) then
            call assign_info(info, BIVKLD_INVALID_BANDWIDTH)
            return
         end if
      else
         call select_bandwidth(ys, selector, bandwidth_y, status)
         if (status /= BIVKLD_SUCCESS) then
            call assign_info(info, status)
            return
         end if
      end if

      call pair_grid(xs, ys, bandwidth_x, bandwidth_y, ngrid, limits, axis1, axis2, cell_area, range_used, status)
      if (status /= BIVKLD_SUCCESS) then
         call assign_info(info, status)
         return
      end if
      allocate(px(size(axis1)*size(axis2)), py(size(axis1)*size(axis2)))
      call kde_on_grid(xs, bandwidth_x, axis1, axis2, px, status)
      if (status /= BIVKLD_SUCCESS) then
         call assign_info(info, status)
         return
      end if
      call kde_on_grid(ys, bandwidth_y, axis1, axis2, py, status)
      if (status /= BIVKLD_SUCCESS) then
         call assign_info(info, status)
         return
      end if

      estimate = integrated_kld(px, py, cell_area)
      if (present(details)) then
         details%estimate = estimate
         details%bandwidth_x = bandwidth_x
         details%bandwidth_y = bandwidth_y
         details%axis1 = axis1
         details%axis2 = axis2
         details%density_x = px
         details%density_y = py
         details%cell_area = cell_area
         details%range = range_used
         details%standardize = std_mode
         details%bandwidth_selector = selector
      end if
   end subroutine biv_kld

   pure subroutine biv_kld_matrix(samples, divergence, H, bandwidth, grid_size, limits, standardize, &
                             selected_H, axis1_out, axis2_out, info)
      type(biv_sample), intent(in) :: samples(:) !! At least two samples; each values component is n-by-2 with n >= 3.
      real(dp), allocatable, intent(out) :: divergence(:, :) !! Directed matrix with element (i,j) equal to D(i || j).
      real(dp), intent(in), optional :: H(:, :, :) !! Optional (2,2,ngroup) bandwidth array, supplied for all groups.
      character(len=*), intent(in), optional :: bandwidth !! Selector name: 'scv' (default) or 'normal'.
      integer, intent(in), optional :: grid_size(2) !! Common grid counts by coordinate; each must be at least 10.
      real(dp), intent(in), optional :: limits(2, 2) !! Common limits, one increasing range per coordinate row.
      character(len=*), intent(in), optional :: standardize !! Standardization mode: 'none', 'pooled', or 'separate'.
      real(dp), allocatable, intent(out), optional :: selected_H(:, :, :) !! Selected bandwidths with shape (2,2,ngroup).
      real(dp), allocatable, intent(out), optional :: axis1_out(:) !! Optional first-coordinate integration grid.
      real(dp), allocatable, intent(out), optional :: axis2_out(:) !! Optional second-coordinate integration grid.
      integer, intent(out), optional :: info !! Status code; zero succeeds and positive values report errors.

      type(biv_sample), allocatable :: work(:)
      real(dp), allocatable :: bandwidths(:, :, :), density(:, :), axis1(:), axis2(:)
      real(dp) :: cell_area, range_used(2, 2)
      integer :: ngrid(2), status, ng, i, j
      character(len=8) :: selector, std_mode

      allocate(divergence(0, 0))
      call assign_info(info, BIVKLD_SUCCESS)
      ng = size(samples)
      if (ng < 2) then
         call assign_info(info, BIVKLD_INVALID_INPUT)
         return
      end if
      do i = 1, ng
         if (.not. allocated(samples(i)%values)) then
            call assign_info(info, BIVKLD_INVALID_INPUT)
            return
         end if
         if (.not. valid_sample(samples(i)%values)) then
            call assign_info(info, BIVKLD_INVALID_INPUT)
            return
         end if
      end do

      call resolve_options(bandwidth, standardize, grid_size, selector, std_mode, ngrid, status)
      if (status /= BIVKLD_SUCCESS) then
         call assign_info(info, status)
         return
      end if
      call standardize_groups(samples, std_mode, work, status)
      if (status /= BIVKLD_SUCCESS) then
         call assign_info(info, status)
         return
      end if

      allocate(bandwidths(2, 2, ng))
      if (present(H)) then
         if (size(H, 1) /= 2 .or. size(H, 2) /= 2 .or. size(H, 3) /= ng) then
            call assign_info(info, BIVKLD_INVALID_BANDWIDTH)
            return
         end if
         do i = 1, ng
            if (.not. is_spd_2x2(H(:, :, i))) then
               call assign_info(info, BIVKLD_INVALID_BANDWIDTH)
               return
            end if
            bandwidths(:, :, i) = 0.5_dp*(H(:, :, i) + transpose(H(:, :, i)))
         end do
      else
         do i = 1, ng
            call select_bandwidth(work(i)%values, selector, bandwidths(:, :, i), status)
            if (status /= BIVKLD_SUCCESS) then
               call assign_info(info, status)
               return
            end if
         end do
      end if

      call group_grid(work, bandwidths, ngrid, limits, axis1, axis2, cell_area, range_used, status)
      if (status /= BIVKLD_SUCCESS) then
         call assign_info(info, status)
         return
      end if
      allocate(density(size(axis1)*size(axis2), ng))
      do i = 1, ng
         call kde_on_grid(work(i)%values, bandwidths(:, :, i), axis1, axis2, density(:, i), status)
         if (status /= BIVKLD_SUCCESS) then
            call assign_info(info, status)
            return
         end if
      end do

      deallocate(divergence)
      allocate(divergence(ng, ng))
      divergence = 0.0_dp
      do j = 1, ng
         do i = 1, ng
            if (i /= j) divergence(i, j) = integrated_kld(density(:, i), density(:, j), cell_area)
         end do
      end do
      if (present(selected_H)) then
         allocate(selected_H(2, 2, ng))
         selected_H = bandwidths
      end if
      if (present(axis1_out)) axis1_out = axis1
      if (present(axis2_out)) axis2_out = axis2
   end subroutine biv_kld_matrix

   pure subroutine hns_bandwidth(x, H, info)
      real(dp), intent(in) :: x(:, :) !! Bivariate sample, with observations in rows and exactly two columns.
      real(dp), intent(out) :: H(2, 2) !! Normal-scale bandwidth matrix using the bivariate Hns scaling rule.
      integer, intent(out), optional :: info !! Status; zero succeeds and positive values report invalid input.

      real(dp) :: covariance(2, 2), factor
      integer :: status

      H = 0.0_dp
      if (.not. valid_sample(x)) then
         call assign_info(info, BIVKLD_INVALID_INPUT)
         return
      end if
      call sample_covariance(x, covariance, status)
      if (status /= BIVKLD_SUCCESS) then
         call assign_info(info, status)
         return
      end if
      factor = real(size(x, 1), dp)**(-1.0_dp/3.0_dp)
      H = factor*covariance
      call assign_info(info, BIVKLD_SUCCESS)
   end subroutine hns_bandwidth

   pure subroutine hscv_bandwidth(x, H, info)
      real(dp), intent(in) :: x(:, :) !! Bivariate sample, with observations in rows and exactly two columns.
      real(dp), intent(out) :: H(2, 2) !! Full symmetric positive-definite bandwidth minimizing the local SCV criterion.
      integer, intent(out), optional :: info !! Status; zero succeeds and positive values report failure.

      real(dp) :: covariance(2, 2), pilot(2, 2), start_H(2, 2), p0(3), pbest(3), factor
      integer :: status

      H = 0.0_dp
      if (.not. valid_sample(x)) then
         call assign_info(info, BIVKLD_INVALID_INPUT)
         return
      end if
      call sample_covariance(x, covariance, status)
      if (status /= BIVKLD_SUCCESS) then
         call assign_info(info, status)
         return
      end if
      factor = (1.0_dp/(3.0_dp*real(size(x, 1), dp)))**0.25_dp
      pilot = factor*covariance
      call hns_bandwidth(x, start_H, status)
      if (status /= BIVKLD_SUCCESS) then
         call assign_info(info, status)
         return
      end if
      call bandwidth_to_parameters(start_H, p0)
      call nelder_mead_scv(x, pilot, p0, pbest, status)
      call parameters_to_bandwidth(pbest, H)
      if (.not. is_spd_2x2(H)) status = BIVKLD_NUMERICAL_FAILURE
      call assign_info(info, status)
   end subroutine hscv_bandwidth

   pure subroutine kde_density(x, H, points, density, info)
      real(dp), intent(in) :: x(:, :) !! Bivariate sample, with observations in rows and exactly two columns.
      real(dp), intent(in) :: H(2, 2) !! Symmetric positive-definite Gaussian kernel bandwidth matrix.
      real(dp), intent(in) :: points(:, :) !! Evaluation points, one point per row and exactly two columns.
      real(dp), intent(out) :: density(:) !! KDE values; its length must equal the number of evaluation points.
      integer, intent(out), optional :: info !! Status; zero succeeds and positive values report errors.

      real(dp) :: inverse_H(2, 2), det_H, delta1, delta2, quadratic, normalizer
      integer :: i, j

      density = 0.0_dp
      if (.not. valid_sample(x) .or. size(points, 2) /= 2 .or. size(density) /= size(points, 1)) then
         call assign_info(info, BIVKLD_INVALID_INPUT)
         return
      end if
      if (any(.not. ieee_is_finite(points)) .or. .not. is_spd_2x2(H)) then
         call assign_info(info, BIVKLD_INVALID_BANDWIDTH)
         return
      end if
      det_H = determinant_2x2(H)
      call inverse_2x2(H, inverse_H)
      normalizer = 1.0_dp/(2.0_dp*pi*sqrt(det_H)*real(size(x, 1), dp))
      do i = 1, size(points, 1)
         do j = 1, size(x, 1)
            delta1 = points(i, 1) - x(j, 1)
            delta2 = points(i, 2) - x(j, 2)
            quadratic = delta1*(inverse_H(1, 1)*delta1 + inverse_H(1, 2)*delta2) &
                        + delta2*(inverse_H(2, 1)*delta1 + inverse_H(2, 2)*delta2)
            density(i) = density(i) + exp(-0.5_dp*max(quadratic, 0.0_dp))
         end do
      end do
      density = normalizer*density
      call assign_info(info, BIVKLD_SUCCESS)
   end subroutine kde_density

   pure subroutine select_bandwidth(x, selector, H, status)
      real(dp), intent(in) :: x(:, :) !! Bivariate sample used to choose a bandwidth matrix.
      character(len=*), intent(in) :: selector !! Bandwidth selector name, either 'scv' or 'normal'.
      real(dp), intent(out) :: H(2, 2) !! Selected symmetric positive-definite bandwidth matrix.
      integer, intent(out) :: status !! Status code from the selected bandwidth algorithm.

      select case (trim(selector))
      case ('normal')
         call hns_bandwidth(x, H, status)
      case ('scv')
         call hscv_bandwidth(x, H, status)
      case default
         H = 0.0_dp
         status = BIVKLD_INVALID_OPTION
      end select
   end subroutine select_bandwidth

   pure subroutine resolve_options(bandwidth, standardize, grid_size, selector, std_mode, ngrid, status)
      character(len=*), intent(in), optional :: bandwidth !! Optional bandwidth selector requested by the caller.
      character(len=*), intent(in), optional :: standardize !! Optional standardization mode requested by the caller.
      integer, intent(in), optional :: grid_size(2) !! Optional grid dimensions; each must contain at least 10 points.
      character(len=8), intent(out) :: selector !! Resolved bandwidth selector, 'scv' or 'normal'.
      character(len=8), intent(out) :: std_mode !! Resolved standardization mode, 'none', 'pooled', or 'separate'.
      integer, intent(out) :: ngrid(2) !! Resolved grid dimensions for the two coordinates.
      integer, intent(out) :: status !! Zero for valid options, otherwise BIVKLD_INVALID_OPTION or BIVKLD_INVALID_INPUT.

      selector = 'scv'
      if (present(bandwidth)) selector = trim(adjustl(bandwidth))
      std_mode = 'none'
      if (present(standardize)) std_mode = trim(adjustl(standardize))
      ngrid = [100, 100]
      if (present(grid_size)) ngrid = grid_size
      status = BIVKLD_SUCCESS
      if (trim(selector) /= 'scv' .and. trim(selector) /= 'normal') status = BIVKLD_INVALID_OPTION
      if (trim(std_mode) /= 'none' .and. trim(std_mode) /= 'pooled' .and. trim(std_mode) /= 'separate') then
         status = BIVKLD_INVALID_OPTION
      end if
      if (any(ngrid < 10)) status = BIVKLD_INVALID_INPUT
   end subroutine resolve_options

   pure subroutine standardize_pair(x, y, mode, xs, ys, status)
      real(dp), intent(in) :: x(:, :) !! First valid bivariate sample.
      real(dp), intent(in) :: y(:, :) !! Second valid bivariate sample.
      character(len=*), intent(in) :: mode !! Standardization mode: 'none', 'pooled', or 'separate'.
      real(dp), allocatable, intent(out) :: xs(:, :) !! Transformed copy of x.
      real(dp), allocatable, intent(out) :: ys(:, :) !! Transformed copy of y.
      integer, intent(out) :: status !! Zero on success, otherwise invalid-input status for zero coordinate variance.

      real(dp) :: center(2), spread(2)

      xs = x
      ys = y
      status = BIVKLD_SUCCESS
      select case (trim(mode))
      case ('none')
         return
      case ('pooled')
         call pooled_pair_moments(x, y, center, spread)
         if (any(spread <= 0.0_dp) .or. any(.not. ieee_is_finite(spread))) then
            status = BIVKLD_INVALID_INPUT
            return
         end if
         call apply_standardization(xs, center, spread)
         call apply_standardization(ys, center, spread)
      case ('separate')
         call sample_moments(xs, center, spread)
         if (any(spread <= 0.0_dp) .or. any(.not. ieee_is_finite(spread))) then
            status = BIVKLD_INVALID_INPUT
            return
         end if
         call apply_standardization(xs, center, spread)
         call sample_moments(ys, center, spread)
         if (any(spread <= 0.0_dp) .or. any(.not. ieee_is_finite(spread))) then
            status = BIVKLD_INVALID_INPUT
            return
         end if
         call apply_standardization(ys, center, spread)
      case default
         status = BIVKLD_INVALID_OPTION
      end select
   end subroutine standardize_pair

   pure subroutine standardize_groups(samples, mode, work, status)
      type(biv_sample), intent(in) :: samples(:) !! Valid input samples to copy and optionally standardize.
      character(len=*), intent(in) :: mode !! Standardization mode: 'none', 'pooled', or 'separate'.
      type(biv_sample), allocatable, intent(out) :: work(:) !! Transformed copies of all input samples.
      integer, intent(out) :: status !! Zero on success, otherwise invalid-input status for zero coordinate variance.

      real(dp) :: center(2), spread(2)
      integer :: i

      allocate(work(size(samples)))
      do i = 1, size(samples)
         work(i)%values = samples(i)%values
      end do
      status = BIVKLD_SUCCESS
      select case (trim(mode))
      case ('none')
         return
      case ('pooled')
         call pooled_group_moments(work, center, spread)
         if (any(spread <= 0.0_dp) .or. any(.not. ieee_is_finite(spread))) then
            status = BIVKLD_INVALID_INPUT
            return
         end if
         do i = 1, size(work)
            call apply_standardization(work(i)%values, center, spread)
         end do
      case ('separate')
         do i = 1, size(work)
            call sample_moments(work(i)%values, center, spread)
            if (any(spread <= 0.0_dp) .or. any(.not. ieee_is_finite(spread))) then
               status = BIVKLD_INVALID_INPUT
               return
            end if
            call apply_standardization(work(i)%values, center, spread)
         end do
      case default
         status = BIVKLD_INVALID_OPTION
      end select
   end subroutine standardize_groups

   pure subroutine sample_moments(x, center, spread)
      real(dp), intent(in) :: x(:, :) !! Bivariate sample for which column moments are required.
      real(dp), intent(out) :: center(2) !! Arithmetic mean of each coordinate.
      real(dp), intent(out) :: spread(2) !! Unbiased sample standard deviation of each coordinate.

      integer :: j

      do j = 1, 2
         center(j) = sum(x(:, j))/real(size(x, 1), dp)
         spread(j) = sqrt(sum((x(:, j) - center(j))**2)/real(size(x, 1) - 1, dp))
      end do
   end subroutine sample_moments

   pure subroutine pooled_pair_moments(x, y, center, spread)
      real(dp), intent(in) :: x(:, :) !! First bivariate sample contributing to pooled moments.
      real(dp), intent(in) :: y(:, :) !! Second bivariate sample contributing to pooled moments.
      real(dp), intent(out) :: center(2) !! Pooled coordinate means across both samples.
      real(dp), intent(out) :: spread(2) !! Pooled unbiased coordinate standard deviations across all observations.

      integer :: j, ntotal

      ntotal = size(x, 1) + size(y, 1)
      do j = 1, 2
         center(j) = (sum(x(:, j)) + sum(y(:, j)))/real(ntotal, dp)
         spread(j) = sqrt((sum((x(:, j) - center(j))**2) + sum((y(:, j) - center(j))**2)) &
                          /real(ntotal - 1, dp))
      end do
   end subroutine pooled_pair_moments

   pure subroutine pooled_group_moments(samples, center, spread)
      type(biv_sample), intent(in) :: samples(:) !! Collection of bivariate samples contributing to pooled moments.
      real(dp), intent(out) :: center(2) !! Pooled coordinate means across every observation in every group.
      real(dp), intent(out) :: spread(2) !! Pooled unbiased coordinate standard deviations across every observation.

      real(dp) :: sums(2), ss(2)
      integer :: i, j, ntotal

      sums = 0.0_dp
      ntotal = 0
      do i = 1, size(samples)
         ntotal = ntotal + size(samples(i)%values, 1)
         do j = 1, 2
            sums(j) = sums(j) + sum(samples(i)%values(:, j))
         end do
      end do
      center = sums/real(ntotal, dp)
      ss = 0.0_dp
      do i = 1, size(samples)
         do j = 1, 2
            ss(j) = ss(j) + sum((samples(i)%values(:, j) - center(j))**2)
         end do
      end do
      spread = sqrt(ss/real(ntotal - 1, dp))
   end subroutine pooled_group_moments

   pure subroutine apply_standardization(x, center, spread)
      real(dp), intent(inout) :: x(:, :) !! Bivariate sample transformed in place by coordinate centering and scaling.
      real(dp), intent(in) :: center(2) !! Coordinate centers subtracted from each observation.
      real(dp), intent(in) :: spread(2) !! Positive coordinate scales dividing each centered observation.

      integer :: j

      do j = 1, 2
         x(:, j) = (x(:, j) - center(j))/spread(j)
      end do
   end subroutine apply_standardization

   pure subroutine sample_covariance(x, covariance, status)
      real(dp), intent(in) :: x(:, :) !! Valid bivariate sample whose unbiased covariance matrix is required.
      real(dp), intent(out) :: covariance(2, 2) !! Unbiased two-by-two sample covariance matrix.
      integer, intent(out) :: status !! Zero for positive-definite covariance, otherwise invalid-bandwidth status.

      real(dp) :: center(2), spread_unused(2), z1, z2
      integer :: i

      call sample_moments(x, center, spread_unused)
      covariance = 0.0_dp
      do i = 1, size(x, 1)
         z1 = x(i, 1) - center(1)
         z2 = x(i, 2) - center(2)
         covariance(1, 1) = covariance(1, 1) + z1*z1
         covariance(1, 2) = covariance(1, 2) + z1*z2
         covariance(2, 2) = covariance(2, 2) + z2*z2
      end do
      covariance = covariance/real(size(x, 1) - 1, dp)
      covariance(2, 1) = covariance(1, 2)
      if (is_spd_2x2(covariance)) then
         status = BIVKLD_SUCCESS
      else
         status = BIVKLD_INVALID_BANDWIDTH
      end if
   end subroutine sample_covariance

   pure subroutine pair_grid(x, y, Hx, Hy, ngrid, limits, axis1, axis2, cell_area, range_used, status)
      real(dp), intent(in) :: x(:, :) !! First transformed bivariate sample used to determine default endpoints.
      real(dp), intent(in) :: y(:, :) !! Second transformed bivariate sample used to determine default endpoints.
      real(dp), intent(in) :: Hx(2, 2) !! Bandwidth matrix for x, used in default padding.
      real(dp), intent(in) :: Hy(2, 2) !! Bandwidth matrix for y, used in default padding.
      integer, intent(in) :: ngrid(2) !! Number of equally spaced integration points in each coordinate.
      real(dp), intent(in), optional :: limits(2, 2) !! Optional explicit coordinate ranges, one row per coordinate.
      real(dp), allocatable, intent(out) :: axis1(:) !! Equally spaced first-coordinate grid.
      real(dp), allocatable, intent(out) :: axis2(:) !! Equally spaced second-coordinate grid.
      real(dp), intent(out) :: cell_area !! Product of adjacent coordinate spacings used by rectangular integration.
      real(dp), intent(out) :: range_used(2, 2) !! Actual coordinate limits used to construct the grid.
      integer, intent(out) :: status !! Zero on success, otherwise invalid-input status for malformed limits.

      real(dp) :: lo, hi, span, kernel_sd, padding
      integer :: j

      if (present(limits)) then
         if (.not. valid_limits(limits)) then
            status = BIVKLD_INVALID_INPUT
            return
         end if
         range_used = limits
      else
         do j = 1, 2
            lo = min(minval(x(:, j)), minval(y(:, j)))
            hi = max(maxval(x(:, j)), maxval(y(:, j)))
            span = hi - lo
            kernel_sd = max(sqrt(Hx(j, j)), sqrt(Hy(j, j)))
            padding = max(0.15_dp*span, 3.0_dp*kernel_sd, sqrt(epsilon(1.0_dp)))
            range_used(j, 1) = lo - padding
            range_used(j, 2) = hi + padding
         end do
      end if
      call make_axes(range_used, ngrid, axis1, axis2, cell_area)
      status = BIVKLD_SUCCESS
   end subroutine pair_grid

   pure subroutine group_grid(samples, bandwidths, ngrid, limits, axis1, axis2, cell_area, range_used, status)
      type(biv_sample), intent(in) :: samples(:) !! Transformed groups used to determine default endpoints.
      real(dp), intent(in) :: bandwidths(:, :, :) !! Positive-definite group bandwidth matrices, shape (2,2,ngroups).
      integer, intent(in) :: ngrid(2) !! Number of equally spaced integration points in each coordinate.
      real(dp), intent(in), optional :: limits(2, 2) !! Optional explicit coordinate ranges, one row per coordinate.
      real(dp), allocatable, intent(out) :: axis1(:) !! Equally spaced first-coordinate grid.
      real(dp), allocatable, intent(out) :: axis2(:) !! Equally spaced second-coordinate grid.
      real(dp), intent(out) :: cell_area !! Product of adjacent coordinate spacings used by rectangular integration.
      real(dp), intent(out) :: range_used(2, 2) !! Actual common coordinate limits used by all groups.
      integer, intent(out) :: status !! Zero on success, otherwise invalid-input status for malformed limits.

      real(dp) :: lo, hi, span, kernel_sd, padding
      integer :: i, j

      if (present(limits)) then
         if (.not. valid_limits(limits)) then
            status = BIVKLD_INVALID_INPUT
            return
         end if
         range_used = limits
      else
         do j = 1, 2
            lo = huge(1.0_dp)
            hi = -huge(1.0_dp)
            kernel_sd = 0.0_dp
            do i = 1, size(samples)
               lo = min(lo, minval(samples(i)%values(:, j)))
               hi = max(hi, maxval(samples(i)%values(:, j)))
               kernel_sd = max(kernel_sd, sqrt(bandwidths(j, j, i)))
            end do
            span = hi - lo
            padding = max(0.15_dp*span, 3.0_dp*kernel_sd, sqrt(epsilon(1.0_dp)))
            range_used(j, 1) = lo - padding
            range_used(j, 2) = hi + padding
         end do
      end if
      call make_axes(range_used, ngrid, axis1, axis2, cell_area)
      status = BIVKLD_SUCCESS
   end subroutine group_grid

   pure subroutine make_axes(range_used, ngrid, axis1, axis2, cell_area)
      real(dp), intent(in) :: range_used(2, 2) !! Increasing coordinate ranges used for equally spaced axes.
      integer, intent(in) :: ngrid(2) !! Number of points requested on each coordinate axis.
      real(dp), allocatable, intent(out) :: axis1(:) !! First coordinate axis, including both endpoints.
      real(dp), allocatable, intent(out) :: axis2(:) !! Second coordinate axis, including both endpoints.
      real(dp), intent(out) :: cell_area !! Rectangular integration cell area from adjacent grid spacings.

      integer :: i

      allocate(axis1(ngrid(1)), axis2(ngrid(2)))
      do i = 1, ngrid(1)
         axis1(i) = range_used(1, 1) + real(i - 1, dp)*(range_used(1, 2) - range_used(1, 1)) &
                    /real(ngrid(1) - 1, dp)
      end do
      do i = 1, ngrid(2)
         axis2(i) = range_used(2, 1) + real(i - 1, dp)*(range_used(2, 2) - range_used(2, 1)) &
                    /real(ngrid(2) - 1, dp)
      end do
      cell_area = (axis1(2) - axis1(1))*(axis2(2) - axis2(1))
   end subroutine make_axes

   pure subroutine kde_on_grid(x, H, axis1, axis2, density, status)
      real(dp), intent(in) :: x(:, :) !! Valid bivariate sample for KDE evaluation.
      real(dp), intent(in) :: H(2, 2) !! Symmetric positive-definite Gaussian bandwidth matrix.
      real(dp), intent(in) :: axis1(:) !! First-coordinate evaluation axis.
      real(dp), intent(in) :: axis2(:) !! Second-coordinate evaluation axis.
      real(dp), intent(out) :: density(:) !! Flattened density values with coordinate 1 varying fastest.
      integer, intent(out) :: status !! Zero on success, otherwise a KDE validation status code.

      real(dp), allocatable :: points(:, :)
      integer :: i, j, k

      if (size(density) /= size(axis1)*size(axis2)) then
         status = BIVKLD_INVALID_INPUT
         return
      end if
      allocate(points(size(axis1)*size(axis2), 2))
      k = 0
      do j = 1, size(axis2)
         do i = 1, size(axis1)
            k = k + 1
            points(k, 1) = axis1(i)
            points(k, 2) = axis2(j)
         end do
      end do
      call kde_density(x, H, points, density, status)
   end subroutine kde_on_grid

   pure function integrated_kld(p, q, cell_area) result(value)
      real(dp), intent(in) :: p(:) !! First density evaluated on the common grid.
      real(dp), intent(in) :: q(:) !! Second density evaluated at the same grid points as p.
      real(dp), intent(in) :: cell_area !! Area represented by each rectangular grid point in the numerical sum.
      real(dp) :: value

      integer :: i

      value = 0.0_dp
      do i = 1, size(p)
         if (p(i) > 0.0_dp) then
            if (q(i) <= 0.0_dp) then
               value = ieee_value(0.0_dp, ieee_positive_inf)
               return
            end if
            value = value + p(i)*(log(p(i)) - log(q(i)))
         end if
      end do
      value = value*cell_area
      if (ieee_is_finite(value) .and. value < 0.0_dp &
          .and. abs(value) < 100.0_dp*epsilon(1.0_dp)) value = 0.0_dp
   end function integrated_kld

   pure function valid_sample(x) result(ok)
      real(dp), intent(in) :: x(:, :) !! Candidate bivariate sample matrix.
      logical :: ok

      ok = size(x, 2) == 2 .and. size(x, 1) >= 3
      if (ok) ok = all(ieee_is_finite(x))
   end function valid_sample

   pure function valid_limits(limits) result(ok)
      real(dp), intent(in) :: limits(2, 2) !! Candidate coordinate limits, one [minimum, maximum] row per coordinate.
      logical :: ok

      ok = all(ieee_is_finite(limits)) .and. all(limits(:, 1) < limits(:, 2))
   end function valid_limits

   pure function determinant_2x2(a) result(det_a)
      real(dp), intent(in) :: a(2, 2) !! Two-by-two matrix whose determinant is required.
      real(dp) :: det_a

      det_a = a(1, 1)*a(2, 2) - a(1, 2)*a(2, 1)
   end function determinant_2x2

   pure subroutine inverse_2x2(a, inverse_a)
      real(dp), intent(in) :: a(2, 2) !! Nonsingular two-by-two matrix to invert.
      real(dp), intent(out) :: inverse_a(2, 2) !! Matrix inverse of a.

      real(dp) :: det_a

      det_a = determinant_2x2(a)
      inverse_a(1, 1) = a(2, 2)/det_a
      inverse_a(1, 2) = -a(1, 2)/det_a
      inverse_a(2, 1) = -a(2, 1)/det_a
      inverse_a(2, 2) = a(1, 1)/det_a
   end subroutine inverse_2x2

   pure function is_spd_2x2(a) result(ok)
      real(dp), intent(in) :: a(2, 2) !! Candidate symmetric positive-definite two-by-two matrix.
      logical :: ok

      real(dp) :: tolerance

      tolerance = sqrt(epsilon(1.0_dp))
      ok = all(ieee_is_finite(a))
      if (.not. ok) return
      ok = abs(a(1, 2) - a(2, 1)) <= tolerance &
           .and. a(1, 1) > 0.0_dp .and. a(2, 2) > 0.0_dp &
           .and. determinant_2x2(a) > 0.0_dp
   end function is_spd_2x2

   pure subroutine bandwidth_to_parameters(H, parameters)
      real(dp), intent(in) :: H(2, 2) !! Symmetric positive-definite bandwidth matrix to parameterize.
      real(dp), intent(out) :: parameters(3) !! Log marginal scales and Fisher-transformed correlation parameter.

      real(dp) :: sd1, sd2, rho

      sd1 = sqrt(H(1, 1))
      sd2 = sqrt(H(2, 2))
      rho = H(1, 2)/(sd1*sd2)
      rho = max(-0.999999_dp, min(0.999999_dp, rho))
      parameters = [log(sd1), atanh(rho), log(sd2)]
   end subroutine bandwidth_to_parameters

   pure subroutine parameters_to_bandwidth(parameters, H)
      real(dp), intent(in) :: parameters(3) !! Log marginal scales and Fisher-transformed correlation parameter.
      real(dp), intent(out) :: H(2, 2) !! Symmetric positive-definite bandwidth matrix represented by parameters.

      real(dp) :: sd1, sd2, rho

      sd1 = exp(max(-40.0_dp, min(40.0_dp, parameters(1))))
      sd2 = exp(max(-40.0_dp, min(40.0_dp, parameters(3))))
      rho = tanh(parameters(2))
      H(1, 1) = sd1*sd1
      H(2, 2) = sd2*sd2
      H(1, 2) = rho*sd1*sd2
      H(2, 1) = H(1, 2)
   end subroutine parameters_to_bandwidth

   pure function scv_objective(x, pilot, parameters) result(value)
      real(dp), intent(in) :: x(:, :) !! Bivariate sample used by the smoothed cross-validation objective.
      real(dp), intent(in) :: pilot(2, 2) !! Positive-definite pilot bandwidth entering the SCV bias estimate.
      real(dp), intent(in) :: parameters(3) !! Candidate log-scale/correlation bandwidth parameters.
      real(dp) :: value

      real(dp) :: H(2, 2), sigma1(2, 2), sigma2(2, 2), sigma3(2, 2), det_H, bias2

      call parameters_to_bandwidth(parameters, H)
      det_H = determinant_2x2(H)
      if (.not. ieee_is_finite(det_H) .or. det_H <= tiny(1.0_dp)) then
         value = huge(1.0_dp)
         return
      end if
      sigma1 = 2.0_dp*H + 2.0_dp*pilot
      sigma2 = H + 2.0_dp*pilot
      sigma3 = 2.0_dp*pilot
      bias2 = gaussian_pair_functional(x, sigma1) - 2.0_dp*gaussian_pair_functional(x, sigma2) &
              + gaussian_pair_functional(x, sigma3)
      bias2 = max(bias2, 0.0_dp)
      value = 1.0_dp/(4.0_dp*pi*real(size(x, 1), dp)*sqrt(det_H)) + bias2
      if (.not. ieee_is_finite(value)) value = huge(1.0_dp)
   end function scv_objective

   pure function gaussian_pair_functional(x, sigma) result(value)
      real(dp), intent(in) :: x(:, :) !! Bivariate sample whose all-pairs Gaussian functional is required.
      real(dp), intent(in) :: sigma(2, 2) !! Positive-definite covariance used for Gaussian pair differences.
      real(dp) :: value

      real(dp) :: inverse_sigma(2, 2), det_sigma, delta1, delta2, quadratic, normalizer
      integer :: i, j, n

      n = size(x, 1)
      det_sigma = determinant_2x2(sigma)
      if (det_sigma <= 0.0_dp) then
         value = huge(1.0_dp)
         return
      end if
      call inverse_2x2(sigma, inverse_sigma)
      normalizer = 1.0_dp/(2.0_dp*pi*sqrt(det_sigma)*real(n*n, dp))
      value = 0.0_dp
      do j = 1, n
         do i = 1, n
            delta1 = x(i, 1) - x(j, 1)
            delta2 = x(i, 2) - x(j, 2)
            quadratic = delta1*(inverse_sigma(1, 1)*delta1 + inverse_sigma(1, 2)*delta2) &
                        + delta2*(inverse_sigma(2, 1)*delta1 + inverse_sigma(2, 2)*delta2)
            value = value + exp(-0.5_dp*max(quadratic, 0.0_dp))
         end do
      end do
      value = normalizer*value
   end function gaussian_pair_functional

   pure subroutine nelder_mead_scv(x, pilot, p0, pbest, status)
      real(dp), intent(in) :: x(:, :) !! Bivariate sample used to evaluate the SCV objective.
      real(dp), intent(in) :: pilot(2, 2) !! Positive-definite pilot bandwidth held fixed during optimization.
      real(dp), intent(in) :: p0(3) !! Initial log-scale/correlation bandwidth parameters.
      real(dp), intent(out) :: pbest(3) !! Best parameter vector found by the Nelder-Mead search.
      integer, intent(out) :: status !! Zero on convergence; numerical-failure status otherwise.

      real(dp) :: simplex(3, 4), f(4), centroid(3), xr(3), xe(3), xc(3), fr, fe, fc
      real(dp), parameter :: alpha = 1.0_dp, gamma = 2.0_dp, rho = 0.5_dp, sigma = 0.5_dp
      real(dp), parameter :: step = 0.20_dp, tolerance = 1.0e-8_dp
      integer, parameter :: max_iter = 400
      integer :: i, iter

      simplex(:, 1) = p0
      do i = 1, 3
         simplex(:, i + 1) = p0
         simplex(i, i + 1) = simplex(i, i + 1) + step
      end do
      do i = 1, 4
         f(i) = scv_objective(x, pilot, simplex(:, i))
      end do
      status = BIVKLD_NUMERICAL_FAILURE
      do iter = 1, max_iter
         call sort_simplex(simplex, f)
         if (maxval(abs(f(2:4) - f(1))) <= tolerance*(1.0_dp + abs(f(1)))) then
            status = BIVKLD_SUCCESS
            exit
         end if
         centroid = sum(simplex(:, 1:3), dim=2)/3.0_dp
         xr = centroid + alpha*(centroid - simplex(:, 4))
         fr = scv_objective(x, pilot, xr)
         if (fr < f(1)) then
            xe = centroid + gamma*(xr - centroid)
            fe = scv_objective(x, pilot, xe)
            if (fe < fr) then
               simplex(:, 4) = xe
               f(4) = fe
            else
               simplex(:, 4) = xr
               f(4) = fr
            end if
         else if (fr < f(3)) then
            simplex(:, 4) = xr
            f(4) = fr
         else
            if (fr < f(4)) then
               xc = centroid + rho*(xr - centroid)
               fc = scv_objective(x, pilot, xc)
               if (fc <= fr) then
                  simplex(:, 4) = xc
                  f(4) = fc
                  cycle
               end if
            else
               xc = centroid - rho*(centroid - simplex(:, 4))
               fc = scv_objective(x, pilot, xc)
               if (fc < f(4)) then
                  simplex(:, 4) = xc
                  f(4) = fc
                  cycle
               end if
            end if
            do i = 2, 4
               simplex(:, i) = simplex(:, 1) + sigma*(simplex(:, i) - simplex(:, 1))
               f(i) = scv_objective(x, pilot, simplex(:, i))
            end do
         end if
      end do
      call sort_simplex(simplex, f)
      pbest = simplex(:, 1)
      if (all(ieee_is_finite(pbest)) .and. ieee_is_finite(f(1))) then
         if (status == BIVKLD_NUMERICAL_FAILURE) status = BIVKLD_SUCCESS
      end if
   end subroutine nelder_mead_scv

   pure subroutine sort_simplex(simplex, f)
      real(dp), intent(inout) :: simplex(3, 4) !! Nelder-Mead simplex columns, reordered so objective values ascend.
      real(dp), intent(inout) :: f(4) !! Objective values corresponding to simplex columns, sorted in ascending order.

      real(dp) :: temp_p(3), temp_f
      integer :: i, j

      do i = 1, 3
         do j = i + 1, 4
            if (f(j) < f(i)) then
               temp_f = f(i)
               f(i) = f(j)
               f(j) = temp_f
               temp_p = simplex(:, i)
               simplex(:, i) = simplex(:, j)
               simplex(:, j) = temp_p
            end if
         end do
      end do
   end subroutine sort_simplex

   pure subroutine assign_info(info, value)
      integer, intent(out), optional :: info !! Optional caller variable receiving the supplied status code.
      integer, intent(in) :: value !! Status code to assign to info when the optional argument is present.

      if (present(info)) info = value
   end subroutine assign_info

end module bivkld_kernel
