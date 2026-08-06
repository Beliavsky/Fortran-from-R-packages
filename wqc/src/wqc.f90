! SPDX-License-Identifier: GPL-3.0-only
module wqc
   use wqc_kinds, only : dp
   use wqc_random, only : normal_rng
   use wqc_statistics, only : mean_value, quantile_type7_sorted, sample_sd, sort_real
   use waveslim, only : mra, mra_result
   use qcsis_mod, only : qc, qc_result
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   implicit none
   private

   type, public :: wqc_pair_result
      real(dp), allocatable :: quantiles(:)
      real(dp), allocatable :: estimated_qc(:, :)
      real(dp), allocatable :: ci_lower(:, :)
      real(dp), allocatable :: ci_upper(:, :)
      character(len=16) :: wavelet = 'la8'
      integer :: levels = 0
      integer :: simulations = 0
   end type wqc_pair_result

   type, public :: wqc_series_result
      character(len=:), allocatable :: name
      type(wqc_pair_result) :: analysis
   end type wqc_series_result

   type, public :: wqc_multi_result
      type(wqc_series_result), allocatable :: series(:)
   end type wqc_multi_result

   public :: dp
   public :: quantile_correlation
   public :: quantile_correlation_analysis
   public :: apply_quantile_correlation

contains

   function quantile_correlation(x, y, quantiles, stat, errmsg) result(rho)
      real(dp), intent(in) :: x(:)
      real(dp), intent(in) :: y(:)
      real(dp), intent(in) :: quantiles(:)
      integer, optional, intent(out) :: stat
      character(len=:), allocatable, optional, intent(out) :: errmsg
      real(dp), allocatable :: rho(:)

      type(qc_result) :: fit
      integer :: status
      character(len=:), allocatable :: message

      allocate(character(len=1) :: message)
      message = ''
      fit = qc(x, y, quantiles, stat=status, errmsg=message)
      if (status /= 0) then
         allocate(rho(0))
      else
         rho = fit%rho
      end if
      call set_status(status, message, stat, errmsg)
   end function quantile_correlation


   function quantile_correlation_analysis(x, y, quantiles, wf, j_levels, n_sim, seed, stat, errmsg) result(fit)
      real(dp), intent(in) :: x(:)
      real(dp), intent(in) :: y(:)
      real(dp), intent(in) :: quantiles(:)
      character(len=*), optional, intent(in) :: wf
      integer, optional, intent(in) :: j_levels
      integer, optional, intent(in) :: n_sim
      integer, optional, intent(in) :: seed
      integer, optional, intent(out) :: stat
      character(len=:), allocatable, optional, intent(out) :: errmsg
      type(wqc_pair_result) :: fit

      type(normal_rng) :: rng
      integer :: status
      character(len=:), allocatable :: message

      call rng%seed(seed)
      call analyze_with_rng(x, y, quantiles, wf, j_levels, n_sim, rng, fit, status, message)
      call set_status(status, message, stat, errmsg)
   end function quantile_correlation_analysis


   function apply_quantile_correlation(data, quantiles, wf, j_levels, n_sim, seed, series_names, stat, errmsg) result(fit)
      real(dp), intent(in) :: data(:, :)
      real(dp), intent(in) :: quantiles(:)
      character(len=*), optional, intent(in) :: wf
      integer, optional, intent(in) :: j_levels
      integer, optional, intent(in) :: n_sim
      integer, optional, intent(in) :: seed
      character(len=*), optional, intent(in) :: series_names(:)
      integer, optional, intent(out) :: stat
      character(len=:), allocatable, optional, intent(out) :: errmsg
      type(wqc_multi_result) :: fit

      type(normal_rng) :: rng
      integer :: i, n_targets, status
      character(len=:), allocatable :: message
      character(len=32) :: generated_name

      call initialize_status(stat, errmsg)
      if (size(data, 1) < 2) then
         allocate(fit%series(0))
         call set_status(1, 'data must contain at least two observations', stat, errmsg)
         return
      else if (size(data, 2) < 2) then
         allocate(fit%series(0))
         call set_status(2, 'data must contain a reference column and at least one target column', stat, errmsg)
         return
      else if (.not. all(ieee_is_finite(data))) then
         allocate(fit%series(0))
         call set_status(3, 'data must contain only finite values', stat, errmsg)
         return
      end if

      n_targets = size(data, 2) - 1
      if (present(series_names)) then
         if (size(series_names) /= n_targets) then
            allocate(fit%series(0))
            call set_status(4, 'series_names must have one entry per target column', stat, errmsg)
            return
         end if
      end if

      allocate(fit%series(n_targets))
      call rng%seed(seed)
      do i = 1, n_targets
         if (present(series_names)) then
            fit%series(i)%name = trim(series_names(i))
         else
            write(generated_name, '(a,i0)') 'series_', i
            fit%series(i)%name = trim(generated_name)
         end if
         call analyze_with_rng(data(:, 1), data(:, i + 1), quantiles, wf, j_levels, n_sim, rng, &
            fit%series(i)%analysis, status, message)
         if (status /= 0) then
            call set_status(status, 'target column '//trim(integer_string(i + 1))//': '//message, stat, errmsg)
            return
         end if
      end do
   end function apply_quantile_correlation


   subroutine analyze_with_rng(x, y, quantiles, wf, j_levels, n_sim, rng, fit, stat, errmsg)
      real(dp), intent(in) :: x(:)
      real(dp), intent(in) :: y(:)
      real(dp), intent(in) :: quantiles(:)
      character(len=*), optional, intent(in) :: wf
      integer, optional, intent(in) :: j_levels
      integer, optional, intent(in) :: n_sim
      type(normal_rng), intent(inout) :: rng
      type(wqc_pair_result), intent(out) :: fit
      integer, intent(out) :: stat
      character(len=:), allocatable, intent(out) :: errmsg

      type(mra_result) :: decomp_x, decomp_y
      real(dp), allocatable :: sim_x(:), sim_y(:), sim_qc(:, :), work(:), rho(:)
      real(dp) :: mean_x, mean_y, sd_x, sd_y
      character(len=16) :: wavelet
      integer :: i, j, levels, simulations, status
      character(len=:), allocatable :: message

      stat = 0
      errmsg = ''
      call validate_pair_inputs(x, y, quantiles, status, message)
      if (status /= 0) then
         stat = status
         errmsg = message
         call empty_pair_result(fit)
         return
      end if

      wavelet = 'la8'
      if (present(wf)) wavelet = trim(wf)
      levels = 8
      if (present(j_levels)) levels = j_levels
      simulations = 1000
      if (present(n_sim)) simulations = n_sim
      if (levels < 1) then
         stat = 6
         errmsg = 'j_levels must be positive'
         call empty_pair_result(fit)
         return
      else if (simulations < 1) then
         stat = 7
         errmsg = 'n_sim must be positive'
         call empty_pair_result(fit)
         return
      end if

      decomp_x = mra(x, wf=wavelet, j_levels=levels, method='modwt', boundary='periodic')
      if (.not. decomp_x%status%ok()) then
         stat = 8
         errmsg = 'waveslim mra failed for x: '//decomp_x%status%message
         call empty_pair_result(fit)
         return
      end if
      decomp_y = mra(y, wf=wavelet, j_levels=levels, method='modwt', boundary='periodic')
      if (.not. decomp_y%status%ok()) then
         stat = 9
         errmsg = 'waveslim mra failed for y: '//decomp_y%status%message
         call empty_pair_result(fit)
         return
      end if

      fit%quantiles = quantiles
      allocate(fit%estimated_qc(levels, size(quantiles)))
      allocate(fit%ci_lower(levels, size(quantiles)))
      allocate(fit%ci_upper(levels, size(quantiles)))
      fit%wavelet = wavelet
      fit%levels = levels
      fit%simulations = simulations
      allocate(sim_x(size(x)), sim_y(size(y)), sim_qc(simulations, size(quantiles)), work(simulations))

      do j = 1, levels
         rho = quantile_correlation(decomp_x%detail(j)%values, decomp_y%detail(j)%values, quantiles, status, message)
         if (status /= 0) then
            stat = 10
            errmsg = 'quantile correlation failed at level '//trim(integer_string(j))//': '//message
            call empty_pair_result(fit)
            return
         end if
         fit%estimated_qc(j, :) = rho

         mean_x = mean_value(decomp_x%detail(j)%values)
         mean_y = mean_value(decomp_y%detail(j)%values)
         sd_x = sample_sd(decomp_x%detail(j)%values)
         sd_y = sample_sd(decomp_y%detail(j)%values)
         if (sd_x <= 0.0_dp .or. sd_y <= 0.0_dp) then
            stat = 11
            errmsg = 'wavelet detail coefficients must have positive sample standard deviations'
            call empty_pair_result(fit)
            return
         end if

         do i = 1, simulations
            call rng%fill_normal(sim_x, mean_x, sd_x)
            call rng%fill_normal(sim_y, mean_y, sd_y)
            rho = quantile_correlation(sim_x, sim_y, quantiles, status, message)
            if (status /= 0) then
               stat = 12
               errmsg = 'simulated quantile correlation failed: '//message
               call empty_pair_result(fit)
               return
            end if
            sim_qc(i, :) = rho
         end do

         do i = 1, size(quantiles)
            work = sim_qc(:, i)
            call sort_real(work)
            fit%ci_lower(j, i) = quantile_type7_sorted(work, 0.025_dp)
            fit%ci_upper(j, i) = quantile_type7_sorted(work, 0.975_dp)
         end do
      end do
   end subroutine analyze_with_rng


   subroutine validate_pair_inputs(x, y, quantiles, stat, errmsg)
      real(dp), intent(in) :: x(:), y(:), quantiles(:)
      integer, intent(out) :: stat
      character(len=:), allocatable, intent(out) :: errmsg

      stat = 0
      errmsg = ''
      if (size(x) /= size(y)) then
         stat = 1
         errmsg = 'x and y must have the same length'
      else if (size(x) < 2) then
         stat = 2
         errmsg = 'x and y must contain at least two observations'
      else if (.not. all(ieee_is_finite(x)) .or. .not. all(ieee_is_finite(y))) then
         stat = 3
         errmsg = 'x and y must contain only finite values'
      else if (size(quantiles) == 0) then
         stat = 4
         errmsg = 'quantiles must contain at least one probability'
      else if (.not. all(ieee_is_finite(quantiles)) .or. any(quantiles <= 0.0_dp) .or. any(quantiles >= 1.0_dp)) then
         stat = 5
         errmsg = 'quantiles must be finite and strictly between zero and one'
      end if
   end subroutine validate_pair_inputs


   subroutine empty_pair_result(fit)
      type(wqc_pair_result), intent(out) :: fit

      allocate(fit%quantiles(0), fit%estimated_qc(0, 0), fit%ci_lower(0, 0), fit%ci_upper(0, 0))
      fit%levels = 0
      fit%simulations = 0
   end subroutine empty_pair_result


   subroutine initialize_status(stat, errmsg)
      integer, optional, intent(out) :: stat
      character(len=:), allocatable, optional, intent(out) :: errmsg

      if (present(stat)) stat = 0
      if (present(errmsg)) errmsg = ''
   end subroutine initialize_status


   subroutine set_status(code, message, stat, errmsg)
      integer, intent(in) :: code
      character(len=*), intent(in) :: message
      integer, optional, intent(out) :: stat
      character(len=:), allocatable, optional, intent(out) :: errmsg

      if (present(stat)) stat = code
      if (present(errmsg)) errmsg = trim(message)
   end subroutine set_status


   function integer_string(value) result(text)
      integer, intent(in) :: value
      character(len=:), allocatable :: text

      character(len=32) :: buffer

      write(buffer, '(i0)') value
      text = trim(buffer)
   end function integer_string

end module wqc
