! SPDX-License-Identifier: GPL-2.0-or-later
! This file is part of a modern Fortran translation of garchx.
! Copyright (C) 2026 translation contributors.
! Original garchx package copyright (C) Genaro Sucarrat.
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 2 of the License, or
! (at your option) any later version.
module garchx_model
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   use garchx_kinds, only : dp
   use garchx_math, only : random_normal_vector, empirical_quantile, mean_value
   use garchx_linalg, only : invert_matrix
   use garchx_optimize, only : bounded_nelder_mead, numerical_hessian
   implicit none
   private

   type, public :: garchx_spec
      integer, allocatable :: arch_lags(:)
      integer, allocatable :: garch_lags(:)
      integer, allocatable :: asym_lags(:)
      integer :: xreg_count = 0
   end type garchx_spec

   type, public :: garchx_fit
      type(garchx_spec) :: spec
      real(dp), allocatable :: par(:)
      real(dp), allocatable :: y(:)
      real(dp), allocatable :: xreg(:, :)
      real(dp), allocatable :: sigma2(:)
      real(dp), allocatable :: residuals(:)
      real(dp), allocatable :: hessian(:, :)
      real(dp), allocatable :: vcov(:, :)
      real(dp) :: objective = huge(1.0_dp)
      real(dp) :: loglik = -huge(1.0_dp)
      real(dp) :: backcast = 0.0_dp
      integer :: objective_mode = 1
      integer :: status = 1
      integer :: iterations = 0
      character(len=16) :: vcov_type = 'ordinary'
   end type garchx_fit

   public :: make_garchx_spec, garchx_parameter_count, garchx_max_lag
   public :: garchx_filter, garchx_filter_derivatives, garchx_objective_value
   public :: garchx_simulate, fit_garchx, garchx_covariance
   public :: garchx_forecast, garchx_quantile_path, refit_garchx
   public :: garchx_asymptotic_covariance, garchx_loglikelihood
contains
   subroutine make_garchx_spec(spec, arch_lags, garch_lags, asym_lags, xreg_count)
      type(garchx_spec), intent(out) :: spec
      integer, intent(in), optional :: arch_lags(:), garch_lags(:), asym_lags(:)
      integer, intent(in), optional :: xreg_count
      if (present(arch_lags)) then
         allocate(spec%arch_lags(size(arch_lags)))
         spec%arch_lags = arch_lags
      else
         allocate(spec%arch_lags(0))
      end if
      if (present(garch_lags)) then
         allocate(spec%garch_lags(size(garch_lags)))
         spec%garch_lags = garch_lags
      else
         allocate(spec%garch_lags(0))
      end if
      if (present(asym_lags)) then
         allocate(spec%asym_lags(size(asym_lags)))
         spec%asym_lags = asym_lags
      else
         allocate(spec%asym_lags(0))
      end if
      if (present(xreg_count)) spec%xreg_count = xreg_count
   end subroutine make_garchx_spec

   pure integer function garchx_parameter_count(spec) result(npar)
      type(garchx_spec), intent(in) :: spec
      npar = 1 + size(spec%arch_lags) + size(spec%garch_lags) + &
             size(spec%asym_lags) + spec%xreg_count
   end function garchx_parameter_count

   pure integer function garchx_max_lag(spec) result(max_lag)
      type(garchx_spec), intent(in) :: spec
      max_lag = 0
      if (size(spec%arch_lags) > 0) max_lag = max(max_lag, maxval(spec%arch_lags))
      if (size(spec%garch_lags) > 0) max_lag = max(max_lag, maxval(spec%garch_lags))
      if (size(spec%asym_lags) > 0) max_lag = max(max_lag, maxval(spec%asym_lags))
   end function garchx_max_lag

   pure integer function garchx_garch_order(spec) result(order_value)
      type(garchx_spec), intent(in) :: spec
      if (size(spec%garch_lags) > 0) then
         order_value = maxval(spec%garch_lags)
      else
         order_value = 0
      end if
   end function garchx_garch_order

   pure logical function valid_spec(spec) result(ok)
      type(garchx_spec), intent(in) :: spec
      ok = spec%xreg_count >= 0
      if (size(spec%arch_lags) > 0) ok = ok .and. all(spec%arch_lags > 0)
      if (size(spec%garch_lags) > 0) ok = ok .and. all(spec%garch_lags > 0)
      if (size(spec%asym_lags) > 0) ok = ok .and. all(spec%asym_lags > 0)
   end function valid_spec

   subroutine parameter_offsets(spec, arch_start, garch_start, asym_start, xreg_start)
      type(garchx_spec), intent(in) :: spec
      integer, intent(out) :: arch_start, garch_start, asym_start, xreg_start
      arch_start = 2
      garch_start = arch_start + size(spec%arch_lags)
      asym_start = garch_start + size(spec%garch_lags)
      xreg_start = asym_start + size(spec%asym_lags)
   end subroutine parameter_offsets

   subroutine validate_inputs(y, spec, pars, xreg, status)
      real(dp), intent(in) :: y(:), pars(:)
      type(garchx_spec), intent(in) :: spec
      real(dp), intent(in), optional :: xreg(:, :)
      integer, intent(out) :: status
      status = 0
      if (.not. valid_spec(spec)) then
         status = 1
      else if (size(pars) /= garchx_parameter_count(spec)) then
         status = 2
      else if (size(y) <= garchx_max_lag(spec)) then
         status = 3
      else if (.not. all(ieee_is_finite(y)) .or. .not. all(ieee_is_finite(pars))) then
         status = 4
      else if (spec%xreg_count > 0) then
         if (.not. present(xreg)) then
            status = 5
         else if (size(xreg, 1) /= size(y) .or. size(xreg, 2) /= spec%xreg_count) then
            status = 6
         else if (.not. all(ieee_is_finite(xreg))) then
            status = 7
         end if
      end if
   end subroutine validate_inputs

   subroutine build_direct_component(y, spec, pars, direct, xreg)
      real(dp), intent(in) :: y(:), pars(:)
      type(garchx_spec), intent(in) :: spec
      real(dp), intent(out) :: direct(:)
      real(dp), intent(in), optional :: xreg(:, :)
      integer :: n, i, j, lag_value
      integer :: arch_start, garch_start, asym_start, xreg_start
      real(dp), allocatable :: y2(:), asym_y2(:)
      real(dp) :: y2mean, asym_mean

      n = size(y)
      allocate(y2(n), asym_y2(n))
      y2 = y*y
      asym_y2 = merge(y2, 0.0_dp, y < 0.0_dp)
      y2mean = mean_value(y2)
      asym_mean = mean_value(asym_y2)
      call parameter_offsets(spec, arch_start, garch_start, asym_start, xreg_start)
      direct = pars(1)
      do j = 1, size(spec%arch_lags)
         lag_value = spec%arch_lags(j)
         do i = 1, n
            if (i > lag_value) then
               direct(i) = direct(i) + pars(arch_start+j-1)*y2(i-lag_value)
            else
               direct(i) = direct(i) + pars(arch_start+j-1)*y2mean
            end if
         end do
      end do
      do j = 1, size(spec%asym_lags)
         lag_value = spec%asym_lags(j)
         do i = 1, n
            if (i > lag_value) then
               direct(i) = direct(i) + pars(asym_start+j-1)*asym_y2(i-lag_value)
            else
               direct(i) = direct(i) + pars(asym_start+j-1)*asym_mean
            end if
         end do
      end do
      if (spec%xreg_count > 0) then
         direct = direct + matmul(xreg, pars(xreg_start:xreg_start+spec%xreg_count-1))
      end if
   end subroutine build_direct_component

   subroutine garchx_filter(y, spec, pars, sigma2, status, xreg, backcast)
      real(dp), intent(in) :: y(:), pars(:)
      type(garchx_spec), intent(in) :: spec
      real(dp), allocatable, intent(out) :: sigma2(:)
      integer, intent(out) :: status
      real(dp), intent(in), optional :: xreg(:, :)
      real(dp), intent(in), optional :: backcast
      integer :: n, i, j, q, lag_value
      integer :: arch_start, garch_start, asym_start, xreg_start
      real(dp) :: initial_variance
      real(dp), allocatable :: direct(:)

      call validate_inputs(y, spec, pars, xreg, status)
      if (status /= 0) then
         allocate(sigma2(0))
         return
      end if
      n = size(y)
      allocate(sigma2(n), direct(n))
      call build_direct_component(y, spec, pars, direct, xreg)
      if (present(backcast)) then
         initial_variance = backcast
      else
         initial_variance = mean_value(y*y)
      end if
      if (initial_variance < 0.0_dp .or. .not. ieee_is_finite(initial_variance)) then
         status = 8
         sigma2 = 0.0_dp
         return
      end if
      q = garchx_garch_order(spec)
      call parameter_offsets(spec, arch_start, garch_start, asym_start, xreg_start)
      if (q == 0) then
         sigma2 = direct
      else
         sigma2 = initial_variance
         do i = q+1, n
            sigma2(i) = direct(i)
            do j = 1, size(spec%garch_lags)
               lag_value = spec%garch_lags(j)
               sigma2(i) = sigma2(i) + pars(garch_start+j-1)*sigma2(i-lag_value)
            end do
         end do
      end if
      if (.not. all(ieee_is_finite(sigma2))) status = 9
   end subroutine garchx_filter

   subroutine garchx_filter_derivatives(y, spec, pars, sigma2, deriv, status, xreg, backcast)
      real(dp), intent(in) :: y(:), pars(:)
      type(garchx_spec), intent(in) :: spec
      real(dp), allocatable, intent(out) :: sigma2(:), deriv(:, :)
      integer, intent(out) :: status
      real(dp), intent(in), optional :: xreg(:, :)
      real(dp), intent(in), optional :: backcast
      integer :: n, k, i, j, q, lag_value
      integer :: arch_start, garch_start, asym_start, xreg_start
      real(dp), allocatable :: direct(:), y2(:), asym_y2(:)
      real(dp) :: y2mean, asym_mean, initial_variance

      call validate_inputs(y, spec, pars, xreg, status)
      if (status /= 0) then
         allocate(sigma2(0), deriv(0, 0))
         return
      end if
      n = size(y)
      k = size(pars)
      allocate(sigma2(n), deriv(n, k), direct(n), y2(n), asym_y2(n))
      y2 = y*y
      asym_y2 = merge(y2, 0.0_dp, y < 0.0_dp)
      y2mean = mean_value(y2)
      asym_mean = mean_value(asym_y2)
      call parameter_offsets(spec, arch_start, garch_start, asym_start, xreg_start)
      call build_direct_component(y, spec, pars, direct, xreg)
      deriv = 0.0_dp
      deriv(:, 1) = 1.0_dp
      do j = 1, size(spec%arch_lags)
         lag_value = spec%arch_lags(j)
         do i = 1, n
            if (i > lag_value) then
               deriv(i, arch_start+j-1) = y2(i-lag_value)
            else
               deriv(i, arch_start+j-1) = y2mean
            end if
         end do
      end do
      do j = 1, size(spec%asym_lags)
         lag_value = spec%asym_lags(j)
         do i = 1, n
            if (i > lag_value) then
               deriv(i, asym_start+j-1) = asym_y2(i-lag_value)
            else
               deriv(i, asym_start+j-1) = asym_mean
            end if
         end do
      end do
      if (spec%xreg_count > 0) then
         deriv(:, xreg_start:xreg_start+spec%xreg_count-1) = xreg
      end if
      if (present(backcast)) then
         initial_variance = backcast
      else
         initial_variance = y2mean
      end if
      q = garchx_garch_order(spec)
      if (q == 0) then
         sigma2 = direct
      else
         sigma2 = initial_variance
         deriv(1:q, :) = 0.0_dp
         do i = q+1, n
            sigma2(i) = direct(i)
            do j = 1, size(spec%garch_lags)
               lag_value = spec%garch_lags(j)
               sigma2(i) = sigma2(i) + pars(garch_start+j-1)*sigma2(i-lag_value)
               deriv(i, :) = deriv(i, :) + pars(garch_start+j-1)*deriv(i-lag_value, :)
               deriv(i, garch_start+j-1) = deriv(i, garch_start+j-1) + sigma2(i-lag_value)
            end do
         end do
      end if
      if (.not. all(ieee_is_finite(sigma2)) .or. .not. all(ieee_is_finite(deriv))) status = 9
   end subroutine garchx_filter_derivatives

   subroutine garchx_objective_value(y, spec, pars, value, status, xreg, backcast, &
                                     objective_mode, sigma2_min, penalty)
      real(dp), intent(in) :: y(:), pars(:)
      type(garchx_spec), intent(in) :: spec
      real(dp), intent(out) :: value
      integer, intent(out) :: status
      real(dp), intent(in), optional :: xreg(:, :), backcast, sigma2_min, penalty
      integer, intent(in), optional :: objective_mode
      integer :: first, mode, xreg_start
      integer :: arch_start, garch_start, asym_start
      real(dp) :: floor_value, penalty_value
      real(dp), allocatable :: sigma2(:), terms(:)
      logical, allocatable :: include(:)

      mode = 1
      if (present(objective_mode)) mode = objective_mode
      floor_value = epsilon(1.0_dp)
      if (present(sigma2_min)) floor_value = sigma2_min
      penalty_value = 1.0e30_dp
      if (present(penalty)) penalty_value = penalty
      if (mode /= 0 .and. mode /= 1) then
         value = penalty_value
         status = 10
         return
      end if
      call parameter_offsets(spec, arch_start, garch_start, asym_start, xreg_start)
      if (spec%xreg_count > 0) then
         if (pars(1)+sum(pars(xreg_start:xreg_start+spec%xreg_count-1)) <= 0.0_dp) then
            value = penalty_value
            status = 11
            return
         end if
      end if
      call garchx_filter(y, spec, pars, sigma2, status, xreg, backcast)
      if (status /= 0) then
         value = penalty_value
         return
      end if
      first = garchx_max_lag(spec)+1
      if (any(sigma2(first:) <= 0.0_dp) .or. .not. all(ieee_is_finite(sigma2(first:)))) then
         value = penalty_value
         status = 12
         return
      end if
      allocate(terms(size(y)-first+1), include(size(y)-first+1))
      terms = y(first:)**2/max(sigma2(first:), floor_value) + log(max(sigma2(first:), floor_value))
      if (mode == 0) then
         include = abs(y(first:)) > tiny(1.0_dp)
         if (count(include) == 0) then
            value = penalty_value
            status = 13
            return
         end if
         value = sum(terms, mask=include)/real(size(terms), dp)
      else
         value = sum(terms)/real(size(terms), dp)
      end if
      status = 0
   end subroutine garchx_objective_value

   subroutine garchx_simulate(n, spec, pars, y, sigma2, innovations, status, xreg, &
                              supplied_innovations, back_innovations, back_sigma2, back_xreg)
      integer, intent(in) :: n
      type(garchx_spec), intent(in) :: spec
      real(dp), intent(in) :: pars(:)
      real(dp), allocatable, intent(out) :: y(:), sigma2(:), innovations(:)
      integer, intent(out) :: status
      real(dp), intent(in), optional :: xreg(:, :), supplied_innovations(:)
      real(dp), intent(in), optional :: back_innovations(:), back_sigma2(:), back_xreg(:, :)
      integer :: m, total, i, t, j, lag_value
      integer :: arch_start, garch_start, asym_start, xreg_start
      real(dp) :: denom, default_variance, z2mean
      real(dp), allocatable :: zext(:), z2ext(:), iext(:), sext(:), xext(:, :)

      status = 0
      if (.not. valid_spec(spec) .or. size(pars) /= garchx_parameter_count(spec) .or. n < 1) then
         status = 1
         allocate(y(0), sigma2(0), innovations(0))
         return
      end if
      if (spec%xreg_count > 0) then
         if (.not. present(xreg)) then
            status = 2
            allocate(y(0), sigma2(0), innovations(0))
            return
         else if (size(xreg, 1) /= n .or. size(xreg, 2) /= spec%xreg_count) then
            status = 3
            allocate(y(0), sigma2(0), innovations(0))
            return
         end if
      end if
      if (present(supplied_innovations)) then
         if (size(supplied_innovations) /= n) then
            status = 4
            allocate(y(0), sigma2(0), innovations(0))
            return
         end if
      end if
      m = garchx_max_lag(spec)
      total = m+n
      allocate(y(n), sigma2(n), innovations(n), zext(total), z2ext(total), iext(total), sext(total))
      allocate(xext(total, max(1, spec%xreg_count)))
      xext = 0.0_dp
      if (present(supplied_innovations)) then
         innovations = supplied_innovations
      else
         call random_normal_vector(innovations)
      end if
      if (m > 0) then
         if (present(back_innovations)) then
            if (size(back_innovations) /= m) then
               status = 5
               return
            end if
            zext(1:m) = back_innovations
         else
            zext(1:m) = 0.0_dp
         end if
      end if
      zext(m+1:) = innovations
      z2mean = mean_value(innovations*innovations)
      if (m > 0) then
         if (present(back_innovations)) then
            z2ext(1:m) = back_innovations*back_innovations
            iext(1:m) = merge(1.0_dp, 0.0_dp, back_innovations < 0.0_dp)
         else
            z2ext(1:m) = z2mean
            iext(1:m) = 0.0_dp
         end if
      end if
      z2ext(m+1:) = innovations*innovations
      iext(m+1:) = merge(1.0_dp, 0.0_dp, innovations < 0.0_dp)
      call parameter_offsets(spec, arch_start, garch_start, asym_start, xreg_start)
      denom = 1.0_dp
      if (size(spec%arch_lags) > 0) denom = denom-sum(pars(arch_start:garch_start-1))
      if (size(spec%garch_lags) > 0) denom = denom-sum(pars(garch_start:asym_start-1))
      if (abs(denom) <= epsilon(1.0_dp)) then
         status = 6
         return
      end if
      default_variance = pars(1)/denom
      if (.not. ieee_is_finite(default_variance)) then
         status = 7
         return
      end if
      if (m > 0) then
         if (present(back_sigma2)) then
            if (size(back_sigma2) /= m) then
               status = 8
               return
            end if
            sext(1:m) = back_sigma2
         else
            sext(1:m) = default_variance
         end if
      end if
      if (spec%xreg_count > 0) then
         if (m > 0) then
            if (present(back_xreg)) then
               if (size(back_xreg, 1) /= m .or. size(back_xreg, 2) /= spec%xreg_count) then
                  status = 9
                  return
               end if
               xext(1:m, :) = back_xreg
            else
               do j = 1, spec%xreg_count
                  xext(1:m, j) = mean_value(xreg(:, j))
               end do
            end if
         end if
         xext(m+1:, :) = xreg
      end if
      do t = 1, n
         i = m+t
         sext(i) = pars(1)
         do j = 1, size(spec%arch_lags)
            lag_value = spec%arch_lags(j)
            sext(i) = sext(i) + pars(arch_start+j-1)*z2ext(i-lag_value)*sext(i-lag_value)
         end do
         do j = 1, size(spec%garch_lags)
            lag_value = spec%garch_lags(j)
            sext(i) = sext(i) + pars(garch_start+j-1)*sext(i-lag_value)
         end do
         do j = 1, size(spec%asym_lags)
            lag_value = spec%asym_lags(j)
            sext(i) = sext(i) + pars(asym_start+j-1)*iext(i-lag_value)* &
                      z2ext(i-lag_value)*sext(i-lag_value)
         end do
         if (spec%xreg_count > 0) then
            sext(i) = sext(i) + dot_product(xext(i, :), &
                      pars(xreg_start:xreg_start+spec%xreg_count-1))
         end if
         if (sext(i) <= 0.0_dp .or. .not. ieee_is_finite(sext(i))) then
            status = 10
            return
         end if
      end do
      sigma2 = sext(m+1:)
      y = sqrt(sigma2)*innovations
   end subroutine garchx_simulate

   subroutine fit_garchx(y, spec, fit, xreg, initial, lower, upper, backcast, &
                         vcov_type, objective_mode, max_iter, rel_tol, bandwidth, kernel_weights)
      real(dp), intent(in) :: y(:)
      type(garchx_spec), intent(in) :: spec
      type(garchx_fit), intent(out) :: fit
      real(dp), intent(in), optional :: xreg(:, :), initial(:), lower(:), upper(:), backcast
      character(len=*), intent(in), optional :: vcov_type
      integer, intent(in), optional :: objective_mode, max_iter
      real(dp), intent(in), optional :: rel_tol, bandwidth, kernel_weights(:)
      integer :: npar, arch_start, garch_start, asym_start, xreg_start, i, local_status
      real(dp) :: penalty_value, tol_value
      real(dp), allocatable :: x0(:), lo(:), hi(:), xbest(:)

      fit%spec = spec
      allocate(fit%y(size(y)))
      fit%y = y
      if (spec%xreg_count > 0 .and. present(xreg)) then
         allocate(fit%xreg(size(xreg, 1), size(xreg, 2)))
         fit%xreg = xreg
      else
         allocate(fit%xreg(size(y), 0))
      end if
      fit%backcast = mean_value(y*y)
      if (present(backcast)) fit%backcast = backcast
      if (present(objective_mode)) fit%objective_mode = objective_mode
      if (present(vcov_type)) fit%vcov_type = adjustl(vcov_type)
      npar = garchx_parameter_count(spec)
      allocate(x0(npar), lo(npar), hi(npar))
      call parameter_offsets(spec, arch_start, garch_start, asym_start, xreg_start)
      x0 = 0.01_dp
      x0(1) = 0.1_dp
      do i = 1, size(spec%arch_lags)
         x0(arch_start+i-1) = 0.1_dp/max(1.0_dp, real(size(spec%arch_lags)*spec%arch_lags(i), dp))
      end do
      do i = 1, size(spec%garch_lags)
         x0(garch_start+i-1) = 0.7_dp/max(1.0_dp, real(size(spec%garch_lags)*spec%garch_lags(i), dp))
      end do
      do i = 1, size(spec%asym_lags)
         x0(asym_start+i-1) = 0.02_dp/real(spec%asym_lags(i), dp)
      end do
      if (present(initial)) then
         if (size(initial) /= npar) then
            fit%status = 20
            return
         end if
         x0 = initial
      end if
      lo = 0.0_dp
      hi = huge(1.0_dp)**0.25_dp
      if (present(lower)) then
         if (size(lower) /= npar) then
            fit%status = 21
            return
         end if
         lo = lower
      end if
      if (present(upper)) then
         if (size(upper) /= npar) then
            fit%status = 22
            return
         end if
         hi = upper
      end if
      call garchx_objective_value(y, spec, x0, penalty_value, local_status, xreg, fit%backcast, &
                                  fit%objective_mode, penalty=1.0e30_dp)
      if (local_status /= 0 .or. .not. ieee_is_finite(penalty_value)) penalty_value = 1.0e20_dp
      tol_value = 1.0e-8_dp
      if (present(rel_tol)) tol_value = rel_tol
      call bounded_nelder_mead(objective_local, x0, lo, hi, xbest, fit%objective, &
                               fit%status, fit%iterations, max_iter, tol_value)
      allocate(fit%par(npar))
      fit%par = xbest
      call garchx_filter(y, spec, fit%par, fit%sigma2, local_status, xreg, fit%backcast)
      if (local_status /= 0) then
         fit%status = 23
         return
      end if
      allocate(fit%residuals(size(y)))
      fit%residuals = y/sqrt(max(fit%sigma2, epsilon(1.0_dp)))
      fit%loglik = garchx_loglikelihood(y, fit%sigma2, fit%objective_mode, garchx_max_lag(spec))
      call numerical_hessian(objective_local, fit%par, fit%hessian)
      call garchx_covariance(y, spec, fit%par, fit%hessian, trim(fit%vcov_type), fit%vcov, &
                            local_status, xreg, fit%backcast, fit%objective_mode, bandwidth, kernel_weights)
      if (local_status /= 0 .and. fit%status == 0) fit%status = 24
   contains
      function objective_local(pars_local) result(value)
         real(dp), intent(in) :: pars_local(:)
         real(dp) :: value
         integer :: objective_status
         call garchx_objective_value(y, spec, pars_local, value, objective_status, xreg, &
                                     fit%backcast, fit%objective_mode, penalty=penalty_value)
         if (objective_status /= 0) value = penalty_value
      end function objective_local
   end subroutine fit_garchx

   real(dp) function garchx_loglikelihood(y, sigma2, objective_mode, max_lag) result(value)
      real(dp), intent(in) :: y(:), sigma2(:)
      integer, intent(in) :: objective_mode, max_lag
      integer :: first
      real(dp), parameter :: log_two_pi = log(2.0_dp*acos(-1.0_dp))
      logical, allocatable :: include(:)
      real(dp), allocatable :: terms(:)
      first = max_lag+1
      allocate(terms(size(y)-first+1))
      terms = -0.5_dp*(log_two_pi+log(sigma2(first:))+y(first:)**2/sigma2(first:))
      if (objective_mode == 0) then
         allocate(include(size(terms)))
         include = abs(y(first:)) > tiny(1.0_dp)
         value = sum(terms, mask=include)
      else
         value = sum(terms)
      end if
   end function garchx_loglikelihood

   subroutine garchx_covariance(y, spec, pars, hessian, vcov_type, vcov, status, xreg, &
                                backcast, objective_mode, bandwidth, kernel_weights)
      real(dp), intent(in) :: y(:), pars(:), hessian(:, :)
      type(garchx_spec), intent(in) :: spec
      character(len=*), intent(in) :: vcov_type
      real(dp), allocatable, intent(out) :: vcov(:, :)
      integer, intent(out) :: status
      real(dp), intent(in), optional :: xreg(:, :), backcast, bandwidth, kernel_weights(:)
      integer, intent(in), optional :: objective_mode
      integer :: first, n_eff, k, i, j, lag_count, mode, inv_status
      real(dp) :: kappa, gamma_n, weight
      real(dp), allocatable :: sigma2(:), deriv(:, :), hinv(:, :), sample_resid(:)
      real(dp), allocatable :: meat(:, :), score(:, :), weighted_deriv(:, :)
      logical, allocatable :: include(:)

      mode = 1
      if (present(objective_mode)) mode = objective_mode
      call garchx_filter_derivatives(y, spec, pars, sigma2, deriv, status, xreg, backcast)
      if (status /= 0) then
         allocate(vcov(0, 0))
         return
      end if
      call invert_matrix(hessian, hinv, inv_status)
      if (inv_status /= 0) then
         status = 30
         allocate(vcov(size(pars), size(pars)))
         vcov = 0.0_dp
         return
      end if
      first = garchx_max_lag(spec)+1
      n_eff = size(y)-first+1
      k = size(pars)
      allocate(vcov(k, k), sample_resid(n_eff), include(n_eff))
      sample_resid = y(first:)/sqrt(sigma2(first:))
      include = .true.
      if (mode == 0) include = abs(y(first:)) > tiny(1.0_dp)
      select case (trim(adjustl(vcov_type)))
      case ('ordinary')
         if (count(include) == 0) then
            status = 31
            vcov = 0.0_dp
            return
         end if
         kappa = sum(sample_resid**4, mask=include)/real(count(include), dp)
         vcov = (kappa-1.0_dp)*hinv/real(n_eff, dp)
      case ('robust')
         allocate(weighted_deriv(n_eff, k), meat(k, k))
         do i = 1, n_eff
            weighted_deriv(i, :) = y(first+i-1)**2/(sigma2(first+i-1)**2)*deriv(first+i-1, :)
         end do
         meat = matmul(transpose(weighted_deriv), weighted_deriv)/real(n_eff, dp)-hessian
         vcov = matmul(hinv, matmul(meat, hinv))/real(n_eff, dp)
      case ('hac')
         allocate(score(n_eff, k), meat(k, k))
         do i = 1, n_eff
            score(i, :) = (1.0_dp/sigma2(first+i-1)- &
                           y(first+i-1)**2/sigma2(first+i-1)**2)*deriv(first+i-1, :)
         end do
         meat = matmul(transpose(score), score)/real(n_eff, dp)
         if (present(bandwidth)) then
            if (bandwidth < 0.0_dp) then
               status = 32
               vcov = 0.0_dp
               return
            end if
            gamma_n = bandwidth
         else
            gamma_n = 4.0_dp*(real(size(y), dp)/100.0_dp)**(2.0_dp/9.0_dp)
         end if
         lag_count = floor(gamma_n)
         do j = 1, lag_count
            if (present(kernel_weights) .and. j <= size(kernel_weights)) then
               weight = kernel_weights(j)
            else
               weight = 1.0_dp-real(j, dp)/gamma_n
            end if
            meat = meat + weight*( &
               matmul(transpose(score(j+1:n_eff, :)), score(1:n_eff-j, :))/real(n_eff-j, dp) + &
               matmul(transpose(score(1:n_eff-j, :)), score(j+1:n_eff, :))/real(n_eff-j, dp))
         end do
         vcov = matmul(hinv, matmul(meat, hinv))/real(n_eff, dp)
      case default
         status = 33
         vcov = 0.0_dp
         return
      end select
      status = 0
   end subroutine garchx_covariance

   subroutine garchx_forecast(fit, n_ahead, forecast, status, future_xreg, n_sim, paths)
      type(garchx_fit), intent(in) :: fit
      integer, intent(in) :: n_ahead
      real(dp), allocatable, intent(out) :: forecast(:)
      integer, intent(out) :: status
      real(dp), intent(in), optional :: future_xreg(:, :)
      integer, intent(in), optional :: n_sim
      real(dp), allocatable, intent(out), optional :: paths(:, :)
      integer :: simulations, s, t, m, draw_index, first
      real(dp) :: u
      real(dp), allocatable :: z(:), ysim(:), ssim(:), zout(:), all_paths(:, :)
      real(dp), allocatable :: back_z(:), back_s(:), back_x(:, :)

      status = 0
      if (n_ahead < 1) then
         status = 40
         allocate(forecast(0))
         return
      end if
      if (fit%spec%xreg_count > 0) then
         if (.not. present(future_xreg)) then
            status = 41
            allocate(forecast(0))
            return
         else if (size(future_xreg, 1) /= n_ahead .or. &
                  size(future_xreg, 2) /= fit%spec%xreg_count) then
            status = 42
            allocate(forecast(0))
            return
         end if
      end if
      simulations = merge(1, 5000, n_ahead == 1)
      if (present(n_sim)) simulations = n_sim
      if (simulations < 1) then
         status = 43
         allocate(forecast(0))
         return
      end if
      allocate(all_paths(n_ahead, simulations), z(n_ahead))
      m = garchx_max_lag(fit%spec)
      first = m+1
      allocate(back_z(m), back_s(m))
      if (m > 0) then
         back_z = fit%residuals(size(fit%residuals)-m+1:)
         back_s = fit%sigma2(size(fit%sigma2)-m+1:)
         if (fit%spec%xreg_count > 0) then
            allocate(back_x(m, fit%spec%xreg_count))
            back_x = fit%xreg(size(fit%xreg, 1)-m+1:, :)
         end if
      end if
      do s = 1, simulations
         do t = 1, n_ahead
            call random_number(u)
            draw_index = first + int(u*real(size(fit%residuals)-first+1, dp))
            draw_index = min(size(fit%residuals), max(first, draw_index))
            z(t) = fit%residuals(draw_index)
         end do
         if (fit%spec%xreg_count > 0) then
            call garchx_simulate(n_ahead, fit%spec, fit%par, ysim, ssim, zout, status, &
                                 xreg=future_xreg, supplied_innovations=z, &
                                 back_innovations=back_z, back_sigma2=back_s, back_xreg=back_x)
         else
            call garchx_simulate(n_ahead, fit%spec, fit%par, ysim, ssim, zout, status, &
                                 supplied_innovations=z, back_innovations=back_z, &
                                 back_sigma2=back_s)
         end if
         if (status /= 0) then
            allocate(forecast(0))
            return
         end if
         all_paths(:, s) = ssim
      end do
      allocate(forecast(n_ahead))
      do t = 1, n_ahead
         forecast(t) = sum(all_paths(t, :))/real(simulations, dp)
      end do
      if (present(paths)) then
         allocate(paths(n_ahead, simulations))
         paths = all_paths
      end if
   end subroutine garchx_forecast

   subroutine garchx_quantile_path(fit, probabilities, quantiles)
      type(garchx_fit), intent(in) :: fit
      real(dp), intent(in) :: probabilities(:)
      real(dp), allocatable, intent(out) :: quantiles(:, :)
      integer :: j, first
      first = garchx_max_lag(fit%spec)+1
      allocate(quantiles(size(fit%y)-first+1, size(probabilities)))
      do j = 1, size(probabilities)
         quantiles(:, j) = sqrt(fit%sigma2(first:))* &
                            empirical_quantile(fit%residuals(first:), probabilities(j))
      end do
   end subroutine garchx_quantile_path

   subroutine refit_garchx(old_fit, new_y, new_fit, status, new_xreg, reestimate, &
                           max_iter, rel_tol)
      type(garchx_fit), intent(in) :: old_fit
      real(dp), intent(in) :: new_y(:)
      type(garchx_fit), intent(out) :: new_fit
      integer, intent(out) :: status
      real(dp), intent(in), optional :: new_xreg(:, :)
      logical, intent(in), optional :: reestimate
      integer, intent(in), optional :: max_iter
      real(dp), intent(in), optional :: rel_tol
      logical :: do_reestimate
      integer :: filter_status

      do_reestimate = .false.
      if (present(reestimate)) do_reestimate = reestimate
      if (do_reestimate) then
         call fit_garchx(new_y, old_fit%spec, new_fit, new_xreg, old_fit%par, &
                         backcast=mean_value(new_y*new_y), vcov_type=trim(old_fit%vcov_type), &
                         objective_mode=old_fit%objective_mode, max_iter=max_iter, rel_tol=rel_tol)
         status = new_fit%status
      else
         new_fit%spec = old_fit%spec
         allocate(new_fit%par(size(old_fit%par)), new_fit%y(size(new_y)))
         new_fit%par = old_fit%par
         new_fit%y = new_y
         new_fit%backcast = mean_value(new_y*new_y)
         new_fit%objective_mode = old_fit%objective_mode
         new_fit%vcov_type = old_fit%vcov_type
         if (old_fit%spec%xreg_count > 0 .and. present(new_xreg)) then
            allocate(new_fit%xreg(size(new_xreg, 1), size(new_xreg, 2)))
            new_fit%xreg = new_xreg
         else
            allocate(new_fit%xreg(size(new_y), 0))
         end if
         call garchx_filter(new_y, old_fit%spec, old_fit%par, new_fit%sigma2, &
                            filter_status, new_xreg, new_fit%backcast)
         if (filter_status /= 0) then
            status = filter_status
            new_fit%status = filter_status
            return
         end if
         allocate(new_fit%residuals(size(new_y)))
         new_fit%residuals = new_y/sqrt(new_fit%sigma2)
         new_fit%objective = 0.0_dp
         call garchx_objective_value(new_y, old_fit%spec, old_fit%par, new_fit%objective, &
                                     filter_status, new_xreg, new_fit%backcast, old_fit%objective_mode)
         new_fit%loglik = garchx_loglikelihood(new_y, new_fit%sigma2, old_fit%objective_mode, &
                                               garchx_max_lag(old_fit%spec))
         allocate(new_fit%hessian(size(old_fit%hessian, 1), size(old_fit%hessian, 2)))
         allocate(new_fit%vcov(size(old_fit%vcov, 1), size(old_fit%vcov, 2)))
         new_fit%hessian = old_fit%hessian
         new_fit%vcov = old_fit%vcov
         new_fit%status = 0
         status = 0
      end if
   end subroutine refit_garchx

   subroutine garchx_asymptotic_covariance(pars, spec, n, avar, status, xreg, &
                                            innovations, e_eta4, objective_mode, backcast)
      real(dp), intent(in) :: pars(:)
      type(garchx_spec), intent(in) :: spec
      integer, intent(in) :: n
      real(dp), allocatable, intent(out) :: avar(:, :)
      integer, intent(out) :: status
      real(dp), intent(in), optional :: xreg(:, :), innovations(:), e_eta4, backcast
      integer, intent(in), optional :: objective_mode
      real(dp), allocatable :: y(:), sigma2(:), z(:), hessian(:, :), hinv(:, :)
      real(dp) :: fourth, back_value
      integer :: inv_status, mode

      mode = 1
      if (present(objective_mode)) mode = objective_mode
      call garchx_simulate(n, spec, pars, y, sigma2, z, status, xreg, innovations)
      if (status /= 0) then
         allocate(avar(0, 0))
         return
      end if
      back_value = mean_value(y*y)
      if (present(backcast)) back_value = backcast
      call numerical_hessian(objective_local, pars, hessian)
      call invert_matrix(hessian, hinv, inv_status)
      if (inv_status /= 0) then
         status = 50
         allocate(avar(size(pars), size(pars)))
         avar = 0.0_dp
         return
      end if
      if (present(e_eta4)) then
         fourth = e_eta4
      else
         if (mode == 0) then
            fourth = sum(z**4, mask=abs(z) > tiny(1.0_dp))/real(count(abs(z) > tiny(1.0_dp)), dp)
         else
            fourth = sum(z**4)/real(n, dp)
         end if
      end if
      allocate(avar(size(pars), size(pars)))
      avar = (fourth-1.0_dp)*hinv
      status = 0
   contains
      function objective_local(pars_local) result(value)
         real(dp), intent(in) :: pars_local(:)
         real(dp) :: value
         integer :: local_status
         call garchx_objective_value(y, spec, pars_local, value, local_status, xreg, &
                                     back_value, mode, penalty=1.0e30_dp)
         if (local_status /= 0) value = 1.0e30_dp
      end function objective_local
   end subroutine garchx_asymptotic_covariance
end module garchx_model
