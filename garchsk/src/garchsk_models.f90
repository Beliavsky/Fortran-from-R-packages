! SPDX-License-Identifier: GPL-2.0-or-later
! Derived from GARCHSK 0.1.0, Copyright (C) 2021 Kei Nakagawa.
module garchsk_models
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   use garchsk_kinds, only : dp
   use garchsk_stats, only : mean_value, sample_variance, sample_skewness, sample_kurtosis, all_finite
   use garchsk_types, only : moment_path, forecast_result
   implicit none
   private
   public :: garchsk_construct, gjrsk_construct
   public :: garchsk_negative_log_likelihood, gjrsk_negative_log_likelihood
   public :: garchsk_constraints, gjrsk_constraints
   public :: garchsk_parameters_valid, gjrsk_parameters_valid
   public :: garchsk_forecast, gjrsk_forecast

contains

   function garchsk_construct(params, data) result(path)
      real(dp), intent(in) :: params(:), data(:)
      type(moment_path) :: path
      real(dp) :: u, eta
      integer :: n, t

      n = size(data)
      allocate(path%mu(n), path%h(n), path%skewness(n), path%kurtosis(n))
      if (size(params) /= 10 .or. n < 2 .or. .not. all_finite(data)) then
         path%message = 'garchsk_construct: invalid dimensions or non-finite data'
         return
      end if
      path%mu = mean_value(data)
      path%h = sample_variance(data)
      path%skewness = sample_skewness(data)
      path%kurtosis = sample_kurtosis(data)
      if (path%h(1) <= tiny(1.0_dp)) then
         path%message = 'garchsk_construct: data variance must be positive'
         return
      end if
      do t = 2, n
         path%mu(t) = params(1) * data(t-1)
         u = data(t-1) - path%mu(t-1)
         eta = u / sqrt(path%h(t-1))
         path%h(t) = params(2) + params(3) * u**2 + params(4) * path%h(t-1)
         path%skewness(t) = params(5) + params(6) * eta**3 + params(7) * path%skewness(t-1)
         path%kurtosis(t) = params(8) + params(9) * eta**4 + params(10) * path%kurtosis(t-1)
         if (path%h(t) <= tiny(1.0_dp) .or. .not. ieee_is_finite(path%h(t))) then
            path%message = 'garchsk_construct: non-positive conditional variance'
            return
         end if
      end do
      path%success = .true.
      path%message = 'ok'
   end function garchsk_construct

   function gjrsk_construct(params, data) result(path)
      real(dp), intent(in) :: params(:), data(:)
      type(moment_path) :: path
      real(dp) :: u, eta, indicator
      integer :: n, t

      n = size(data)
      allocate(path%mu(n), path%h(n), path%skewness(n), path%kurtosis(n))
      if (size(params) /= 13 .or. n < 2 .or. .not. all_finite(data)) then
         path%message = 'gjrsk_construct: invalid dimensions or non-finite data'
         return
      end if
      path%mu = mean_value(data)
      path%h = sample_variance(data)
      path%skewness = sample_skewness(data)
      path%kurtosis = sample_kurtosis(data)
      if (path%h(1) <= tiny(1.0_dp)) then
         path%message = 'gjrsk_construct: data variance must be positive'
         return
      end if
      do t = 2, n
         path%mu(t) = params(1) * data(t-1)
         u = data(t-1) - path%mu(t-1)
         eta = u / sqrt(path%h(t-1))
         indicator = merge(1.0_dp, 0.0_dp, u < 0.0_dp)
         path%h(t) = params(2) + params(3) * u**2 + params(4) * u**2 * indicator + params(5) * path%h(t-1)
         path%skewness(t) = params(6) + params(7) * eta**3 + params(8) * eta**3 * indicator + &
            params(9) * path%skewness(t-1)
         path%kurtosis(t) = params(10) + params(11) * eta**4 + params(12) * eta**4 * indicator + &
            params(13) * path%kurtosis(t-1)
         if (path%h(t) <= tiny(1.0_dp) .or. .not. ieee_is_finite(path%h(t))) then
            path%message = 'gjrsk_construct: non-positive conditional variance'
            return
         end if
      end do
      path%success = .true.
      path%message = 'ok'
   end function gjrsk_construct

   function garchsk_negative_log_likelihood(params, data, include_constant) result(value)
      real(dp), intent(in) :: params(:), data(:)
      logical, intent(in), optional :: include_constant
      real(dp) :: value
      type(moment_path) :: path
      real(dp) :: z, expansion, normalizer, term, constant
      logical :: add_constant
      integer :: t

      add_constant = .false.
      if (present(include_constant)) add_constant = include_constant
      value = 1.0e12_dp + constraint_penalty(params, 1)
      if (.not. garchsk_parameters_valid(params)) return
      path = garchsk_construct(params, data)
      if (.not. path%success) return
      constant = merge(log(2.0_dp * acos(-1.0_dp)), 0.0_dp, add_constant)
      value = 0.0_dp
      do t = 2, size(data)
         z = (data(t) - path%mu(t)) / sqrt(path%h(t))
         expansion = 1.0_dp + path%skewness(t) * (z**3 - 3.0_dp*z) / 6.0_dp + &
            (path%kurtosis(t) - 3.0_dp) * (z**4 - 6.0_dp*z**2 + 3.0_dp) / 24.0_dp
         normalizer = 1.0_dp + path%skewness(t)**2 / 6.0_dp + (path%kurtosis(t) - 3.0_dp)**2 / 24.0_dp
         if (abs(expansion) <= tiny(1.0_dp) .or. normalizer <= 0.0_dp) then
            value = 1.0e12_dp
            return
         end if
         term = 0.5_dp * (log(path%h(t)) + z*z + constant) - 2.0_dp*log(abs(expansion)) + log(normalizer)
         if (.not. ieee_is_finite(term)) then
            value = 1.0e12_dp
            return
         end if
         value = value + term
      end do
   end function garchsk_negative_log_likelihood

   function gjrsk_negative_log_likelihood(params, data, include_constant) result(value)
      real(dp), intent(in) :: params(:), data(:)
      logical, intent(in), optional :: include_constant
      real(dp) :: value
      type(moment_path) :: path
      real(dp) :: z, expansion, normalizer, term, constant
      logical :: add_constant
      integer :: t

      add_constant = .false.
      if (present(include_constant)) add_constant = include_constant
      value = 1.0e12_dp + constraint_penalty(params, 2)
      if (.not. gjrsk_parameters_valid(params)) return
      path = gjrsk_construct(params, data)
      if (.not. path%success) return
      constant = merge(log(2.0_dp * acos(-1.0_dp)), 0.0_dp, add_constant)
      value = 0.0_dp
      do t = 2, size(data)
         z = (data(t) - path%mu(t)) / sqrt(path%h(t))
         expansion = 1.0_dp + path%skewness(t) * (z**3 - 3.0_dp*z) / 6.0_dp + &
            (path%kurtosis(t) - 3.0_dp) * (z**4 - 6.0_dp*z**2 + 3.0_dp) / 24.0_dp
         normalizer = 1.0_dp + path%skewness(t)**2 / 6.0_dp + (path%kurtosis(t) - 3.0_dp)**2 / 24.0_dp
         if (abs(expansion) <= tiny(1.0_dp) .or. normalizer <= 0.0_dp) then
            value = 1.0e12_dp
            return
         end if
         term = 0.5_dp * (log(path%h(t)) + z*z + constant) - 2.0_dp*log(abs(expansion)) + log(normalizer)
         if (.not. ieee_is_finite(term)) then
            value = 1.0e12_dp
            return
         end if
         value = value + term
      end do
   end function gjrsk_negative_log_likelihood

   pure function garchsk_constraints(params) result(values)
      real(dp), intent(in) :: params(:)
      real(dp), allocatable :: values(:)
      allocate(values(12))
      if (size(params) /= 10) then
         values = huge(1.0_dp)
         return
      end if
      values = [params(1), params(2), params(3), params(4), params(6), params(7), &
         params(8), params(9), params(10), sum(params(3:4)), sum(params(6:7)), sum(params(9:10))]
   end function garchsk_constraints

   pure function gjrsk_constraints(params) result(values)
      real(dp), intent(in) :: params(:)
      real(dp), allocatable :: values(:)
      allocate(values(15))
      if (size(params) /= 13) then
         values = huge(1.0_dp)
         return
      end if
      values = [params(1), params(2), params(3), params(4), params(5), params(7), params(8), params(9), &
         params(10), params(11), params(12), params(13), sum(params(3:5)), sum(params(7:9)), sum(params(11:13))]
   end function gjrsk_constraints

   pure logical function garchsk_parameters_valid(params) result(ok)
      real(dp), intent(in) :: params(:)
      real(dp), parameter :: margin = 1.0e-10_dp
      ok = .false.
      if (size(params) /= 10) return
      if (.not. all_finite(params)) return
      ok = abs(params(1)) < 1.0_dp - margin .and. params(2) > margin .and. &
         all(params(3:4) >= 0.0_dp) .and. sum(params(3:4)) < 1.0_dp - margin .and. &
         all(abs(params(6:7)) < 1.0_dp) .and. abs(sum(params(6:7))) < 1.0_dp - margin .and. &
         params(8) > margin .and. all(params(9:10) >= 0.0_dp) .and. &
         sum(params(9:10)) < 1.0_dp - margin
   end function garchsk_parameters_valid

   pure logical function gjrsk_parameters_valid(params) result(ok)
      real(dp), intent(in) :: params(:)
      real(dp), parameter :: margin = 1.0e-10_dp
      ok = .false.
      if (size(params) /= 13) return
      if (.not. all_finite(params)) return
      ok = abs(params(1)) < 1.0_dp - margin .and. params(2) > margin .and. &
         all(params(3:5) >= 0.0_dp) .and. sum(params(3:5)) < 1.0_dp - margin .and. &
         all(abs(params(7:9)) < 1.0_dp) .and. abs(sum(params(7:9))) < 1.0_dp - margin .and. &
         params(10) > margin .and. all(params(11:13) >= 0.0_dp) .and. &
         sum(params(11:13)) < 1.0_dp - margin
   end function gjrsk_parameters_valid

   function garchsk_forecast(params, data, max_forecast) result(fcst)
      real(dp), intent(in) :: params(:), data(:)
      integer, intent(in), optional :: max_forecast
      type(forecast_result) :: fcst
      type(moment_path) :: path
      real(dp) :: u, eta
      integer :: horizon, i, n

      horizon = 20
      if (present(max_forecast)) horizon = max_forecast
      if (horizon < 1) then
         fcst%message = 'garchsk_forecast: horizon must be positive'
         return
      end if
      allocate(fcst%mu(horizon), fcst%h(horizon), fcst%skewness(horizon), fcst%kurtosis(horizon))
      path = garchsk_construct(params, data)
      if (.not. path%success) then
         fcst%message = path%message
         return
      end if
      n = size(data)
      u = data(n) - path%mu(n)
      eta = u / sqrt(path%h(n))
      fcst%mu(1) = params(1) * data(n)
      fcst%h(1) = params(2) + params(3) * u**2 + params(4) * path%h(n)
      fcst%skewness(1) = params(5) + params(6) * eta**3 + params(7) * path%skewness(n)
      fcst%kurtosis(1) = params(8) + params(9) * eta**4 + params(10) * path%kurtosis(n)
      do i = 2, horizon
         fcst%mu(i) = params(1) * fcst%mu(i-1)
         fcst%h(i) = params(2) + sum(params(3:4)) * fcst%h(i-1)
         fcst%skewness(i) = params(5) + sum(params(6:7)) * fcst%skewness(i-1)
         fcst%kurtosis(i) = params(8) + sum(params(9:10)) * fcst%kurtosis(i-1)
      end do
      fcst%success = all(fcst%h > 0.0_dp)
      if (fcst%success) then
         fcst%message = 'ok'
      else
         fcst%message = 'garchsk_forecast: invalid variance forecast'
      end if
   end function garchsk_forecast

   function gjrsk_forecast(params, data, max_forecast) result(fcst)
      real(dp), intent(in) :: params(:), data(:)
      integer, intent(in), optional :: max_forecast
      type(forecast_result) :: fcst
      type(moment_path) :: path
      real(dp) :: u, eta, indicator
      integer :: horizon, i, n

      horizon = 20
      if (present(max_forecast)) horizon = max_forecast
      if (horizon < 1) then
         fcst%message = 'gjrsk_forecast: horizon must be positive'
         return
      end if
      allocate(fcst%mu(horizon), fcst%h(horizon), fcst%skewness(horizon), fcst%kurtosis(horizon))
      path = gjrsk_construct(params, data)
      if (.not. path%success) then
         fcst%message = path%message
         return
      end if
      n = size(data)
      u = data(n) - path%mu(n)
      eta = u / sqrt(path%h(n))
      indicator = merge(1.0_dp, 0.0_dp, u < 0.0_dp)
      fcst%mu(1) = params(1) * data(n)
      fcst%h(1) = params(2) + params(3) * u**2 + params(4) * u**2 * indicator + params(5) * path%h(n)
      fcst%skewness(1) = params(6) + params(7) * eta**3 + params(8) * eta**3 * indicator + &
         params(9) * path%skewness(n)
      fcst%kurtosis(1) = params(10) + params(11) * eta**4 + params(12) * eta**4 * indicator + &
         params(13) * path%kurtosis(n)
      do i = 2, horizon
         fcst%mu(i) = params(1) * fcst%mu(i-1)
         fcst%h(i) = params(2) + sum(params(3:5)) * fcst%h(i-1)
         fcst%skewness(i) = params(6) + sum(params(7:9)) * fcst%skewness(i-1)
         fcst%kurtosis(i) = params(10) + sum(params(11:13)) * fcst%kurtosis(i-1)
      end do
      fcst%success = all(fcst%h > 0.0_dp)
      if (fcst%success) then
         fcst%message = 'ok'
      else
         fcst%message = 'gjrsk_forecast: invalid variance forecast'
      end if
   end function gjrsk_forecast

   pure real(dp) function constraint_penalty(params, model) result(value)
      real(dp), intent(in) :: params(:)
      integer, intent(in) :: model
      real(dp) :: violation
      value = 0.0_dp
      if (model == 1) then
         if (size(params) /= 10) then
            value = 1.0e10_dp
            return
         end if
         violation = max(0.0_dp, abs(params(1)) - 0.999_dp)
         value = value + violation**2
         value = value + max(0.0_dp, 1.0e-8_dp - params(2))**2
         value = value + sum(max(0.0_dp, -params(3:4))**2)
         value = value + max(0.0_dp, sum(params(3:4)) - 0.999_dp)**2
         value = value + sum(max(0.0_dp, abs(params(6:7)) - 0.999_dp)**2)
         value = value + max(0.0_dp, abs(sum(params(6:7))) - 0.999_dp)**2
         value = value + max(0.0_dp, 1.0e-8_dp - params(8))**2
         value = value + sum(max(0.0_dp, -params(9:10))**2)
         value = value + max(0.0_dp, sum(params(9:10)) - 0.999_dp)**2
      else
         if (size(params) /= 13) then
            value = 1.0e10_dp
            return
         end if
         violation = max(0.0_dp, abs(params(1)) - 0.999_dp)
         value = value + violation**2
         value = value + max(0.0_dp, 1.0e-8_dp - params(2))**2
         value = value + sum(max(0.0_dp, -params(3:5))**2)
         value = value + max(0.0_dp, sum(params(3:5)) - 0.999_dp)**2
         value = value + sum(max(0.0_dp, abs(params(7:9)) - 0.999_dp)**2)
         value = value + max(0.0_dp, abs(sum(params(7:9))) - 0.999_dp)**2
         value = value + max(0.0_dp, 1.0e-8_dp - params(10))**2
         value = value + sum(max(0.0_dp, -params(11:13))**2)
         value = value + max(0.0_dp, sum(params(11:13)) - 0.999_dp)**2
      end if
      value = 1.0e10_dp * value
   end function constraint_penalty

end module garchsk_models
