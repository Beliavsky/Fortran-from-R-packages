! SPDX-License-Identifier: MIT
module fixedincome_curves
   use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan, ieee_is_nan
   use fixedincome_kinds, only : dp
   use fixedincome_types
   use fixedincome_terms, only : daycount, todays, toyears, term, offset_date
   use fixedincome_compounding, only : compounding, compound, implied_rate
   use fixedincome_interpolation, only : interpolate_points, fit_interpolation_model
   implicit none
   private

   public :: spotrate, spotratecurve, forwardrate
   public :: spot_compound, spot_discount, curve_compound, curve_discount, forward_compound
   public :: set_interpolation, interpolation, prepare_interpolation, interpolate, fit_interpolation
   public :: forwardrate_from_curve, forwardrate_between, as_spotratecurve, as_spotrate, as_forwardrate
   public :: interpolation_error, maturities
   public :: first, last, closest, curve_at_terms, insert_curve_points

   interface spotrate
      module procedure spotrate_scalar
      module procedure spotrate_vector
   end interface spotrate

   interface forwardrate
      module procedure forwardrate_scalar
      module procedure forwardrate_vector
   end interface forwardrate

   interface as_spotrate
      module procedure as_spotrate_curve
      module procedure as_spotrate_forward
   end interface as_spotrate

   interface as_forwardrate
      module procedure as_forwardrate_curve
      module procedure as_forwardrate_spot
   end interface as_forwardrate

   interface spot_compound
      module procedure spot_compound_term
      module procedure spot_compound_days
   end interface spot_compound

   interface spot_discount
      module procedure spot_discount_term
      module procedure spot_discount_days
   end interface spot_discount

contains

   function spotrate_scalar(rate, compounding_name, daycount_name, calendar, status) result(spot)
      real(dp), intent(in) :: rate
      character(len=*), intent(in) :: compounding_name, daycount_name, calendar
      integer, intent(out), optional :: status
      type(spot_rate_t) :: spot
      integer :: stat_i
      allocate(spot%rate(1))
      spot%rate = rate
      spot%compounding = compounding(compounding_name)
      spot%daycount = daycount(daycount_name, stat_i)
      spot%calendar = trim(calendar)
      if (present(status)) then
         if (spot%compounding == 0 .or. stat_i /= FI_OK) then
            status = FI_INVALID_ARGUMENT
         else
            status = FI_OK
         end if
      end if
   end function spotrate_scalar

   function spotrate_vector(rate, compounding_name, daycount_name, calendar, status) result(spot)
      real(dp), intent(in) :: rate(:)
      character(len=*), intent(in) :: compounding_name, daycount_name, calendar
      integer, intent(out), optional :: status
      type(spot_rate_t) :: spot
      integer :: stat_i
      allocate(spot%rate(size(rate)))
      spot%rate = rate
      spot%compounding = compounding(compounding_name)
      spot%daycount = daycount(daycount_name, stat_i)
      spot%calendar = trim(calendar)
      if (present(status)) then
         if (spot%compounding == 0 .or. stat_i /= FI_OK) then
            status = FI_INVALID_ARGUMENT
         else
            status = FI_OK
         end if
      end if
   end function spotrate_vector

   function spotratecurve(rate, terms, compounding_name, daycount_name, calendar, &
                          reference_date, status) result(curve)
      real(dp), intent(in) :: rate(:)
      type(term_t), intent(in) :: terms
      character(len=*), intent(in) :: compounding_name, daycount_name, calendar
      integer, intent(in), optional :: reference_date
      integer, intent(out), optional :: status
      type(spot_rate_curve_t) :: curve
      type(daycount_t) :: dc
      type(term_t) :: days_term
      integer, allocatable :: order(:)
      integer :: i, stat_i

      dc = daycount(daycount_name, stat_i)
      if (size(rate) /= terms%size() .or. stat_i /= FI_OK .or. compounding(compounding_name) == 0) then
         allocate(curve%rate(0), curve%term_days(0))
         if (present(status)) status = FI_INVALID_ARGUMENT
         return
      end if
      days_term = todays(dc, terms, stat_i)
      allocate(order(size(rate)))
      order = [(i, i=1,size(rate))]
      call sort_indices(days_term%value, order)
      allocate(curve%rate(size(rate)), curve%term_days(size(rate)))
      curve%rate = rate(order)
      curve%term_days = days_term%value(order)
      curve%compounding = compounding(compounding_name)
      curve%daycount = dc
      curve%calendar = trim(calendar)
      if (present(reference_date)) curve%reference_date = reference_date
      if (present(status)) then
         if (any(curve%term_days < 0.0_dp) .or. any(curve%term_days(2:) <= curve%term_days(:size(rate)-1))) then
            status = FI_INVALID_ARGUMENT
         else
            status = FI_OK
         end if
      end if
   end function spotratecurve

   function forwardrate_scalar(rate, interval, compounding_name, daycount_name, calendar, &
                               reference_date, status) result(forward)
      real(dp), intent(in) :: rate
      type(term_t), intent(in) :: interval
      character(len=*), intent(in) :: compounding_name, daycount_name, calendar
      integer, intent(in), optional :: reference_date
      integer, intent(out), optional :: status
      type(forward_rate_t) :: forward
      real(dp) :: values(1)
      values = rate
      forward = forwardrate_vector(values, interval, compounding_name, daycount_name, calendar, &
                                   reference_date, status)
   end function forwardrate_scalar

   function forwardrate_vector(rate, interval, compounding_name, daycount_name, calendar, &
                               reference_date, status) result(forward)
      real(dp), intent(in) :: rate(:)
      type(term_t), intent(in) :: interval
      character(len=*), intent(in) :: compounding_name, daycount_name, calendar
      integer, intent(in), optional :: reference_date
      integer, intent(out), optional :: status
      type(forward_rate_t) :: forward
      type(term_t) :: days_term
      integer :: stat_i
      forward%daycount = daycount(daycount_name, stat_i)
      if (size(rate) /= interval%size() .or. stat_i /= FI_OK .or. compounding(compounding_name) == 0) then
         allocate(forward%rate(0), forward%interval_days(0))
         if (present(status)) status = FI_INVALID_ARGUMENT
         return
      end if
      days_term = todays(forward%daycount, interval, stat_i)
      allocate(forward%rate(size(rate)), forward%interval_days(size(rate)))
      forward%rate = rate
      forward%interval_days = days_term%value
      forward%compounding = compounding(compounding_name)
      forward%calendar = trim(calendar)
      if (present(reference_date)) forward%reference_date = reference_date
      if (present(status)) then
         if (any(forward%interval_days < 0.0_dp)) then
            status = FI_INVALID_ARGUMENT
         else
            status = FI_OK
         end if
      end if
   end function forwardrate_vector

   function spot_compound_term(spot, t, status) result(factor)
      type(spot_rate_t), intent(in) :: spot
      type(term_t), intent(in) :: t
      integer, intent(out), optional :: status
      real(dp), allocatable :: factor(:)
      type(term_t) :: years
      integer :: stat_i
      years = toyears(spot%daycount, t, stat_i)
      if (stat_i /= FI_OK) then
         allocate(factor(0))
         if (present(status)) status = stat_i
         return
      end if
      factor = broadcast_compound(spot%compounding, years%value, spot%rate, stat_i)
      if (present(status)) status = stat_i
   end function spot_compound_term

   function spot_compound_days(spot, days, status) result(factor)
      type(spot_rate_t), intent(in) :: spot
      real(dp), intent(in) :: days(:)
      integer, intent(out), optional :: status
      real(dp), allocatable :: factor(:)
      type(term_t) :: t
      t = term(days, 'days')
      factor = spot_compound_term(spot, t, status)
   end function spot_compound_days

   function spot_discount_term(spot, t, status) result(value)
      type(spot_rate_t), intent(in) :: spot
      type(term_t), intent(in) :: t
      integer, intent(out), optional :: status
      real(dp), allocatable :: value(:)
      real(dp), allocatable :: factor(:)
      integer :: stat_i
      factor = spot_compound_term(spot, t, stat_i)
      allocate(value(size(factor)))
      if (stat_i /= FI_OK .or. any(abs(factor) <= tiny(1.0_dp))) then
         value = ieee_value(0.0_dp, ieee_quiet_nan)
         if (present(status)) status = FI_INVALID_ARGUMENT
      else
         value = 1.0_dp / factor
         if (present(status)) status = FI_OK
      end if
   end function spot_discount_term

   function spot_discount_days(spot, days, status) result(value)
      type(spot_rate_t), intent(in) :: spot
      real(dp), intent(in) :: days(:)
      integer, intent(out), optional :: status
      real(dp), allocatable :: value(:)
      type(term_t) :: t
      t = term(days, 'days')
      value = spot_discount_term(spot, t, status)
   end function spot_discount_days

   function curve_compound(curve, status) result(factor)
      type(spot_rate_curve_t), intent(in) :: curve
      integer, intent(out), optional :: status
      real(dp), allocatable :: factor(:)
      real(dp), allocatable :: years(:)
      integer :: stat_i
      years = curve%term_days / real(curve%daycount%days_in_base, dp)
      factor = compound(curve%compounding, years, curve%rate, stat_i)
      if (present(status)) status = stat_i
   end function curve_compound

   function curve_discount(curve, status) result(value)
      type(spot_rate_curve_t), intent(in) :: curve
      integer, intent(out), optional :: status
      real(dp), allocatable :: value(:)
      real(dp), allocatable :: factor(:)
      integer :: stat_i
      factor = curve_compound(curve, stat_i)
      allocate(value(size(factor)))
      if (stat_i /= FI_OK .or. any(abs(factor) <= tiny(1.0_dp))) then
         value = ieee_value(0.0_dp, ieee_quiet_nan)
         if (present(status)) status = FI_INVALID_ARGUMENT
      else
         value = 1.0_dp / factor
         if (present(status)) status = FI_OK
      end if
   end function curve_discount

   function forward_compound(forward, status) result(factor)
      type(forward_rate_t), intent(in) :: forward
      integer, intent(out), optional :: status
      real(dp), allocatable :: factor(:), years(:)
      integer :: stat_i
      years = forward%interval_days / real(forward%daycount%days_in_base, dp)
      factor = compound(forward%compounding, years, forward%rate, stat_i)
      if (present(status)) status = stat_i
   end function forward_compound

   pure function interpolation(curve) result(model)
      type(spot_rate_curve_t), intent(in) :: curve
      type(interpolation_t) :: model
      model = curve%interpolation
   end function interpolation

   function prepare_interpolation(model, curve, status) result(prepared)
      type(interpolation_t), intent(in) :: model
      type(spot_rate_curve_t), intent(in) :: curve
      integer, intent(out), optional :: status
      type(interpolation_t) :: prepared
      prepared = model
      if (curve%size() < 2 .and. model%method < INTERP_NELSON_SIEGEL) then
         if (present(status)) status = FI_INVALID_ARGUMENT
      else
         if (present(status)) status = FI_OK
      end if
   end function prepare_interpolation

   function as_spotrate_curve(curve) result(spot)
      type(spot_rate_curve_t), intent(in) :: curve
      type(spot_rate_t) :: spot
      allocate(spot%rate(curve%size()))
      spot%rate = curve%rate
      spot%compounding = curve%compounding
      spot%daycount = curve%daycount
      spot%calendar = curve%calendar
   end function as_spotrate_curve

   function as_spotrate_forward(forward) result(spot)
      type(forward_rate_t), intent(in) :: forward
      type(spot_rate_t) :: spot
      allocate(spot%rate(forward%size()))
      spot%rate = forward%rate
      spot%compounding = forward%compounding
      spot%daycount = forward%daycount
      spot%calendar = forward%calendar
   end function as_spotrate_forward

   function as_forwardrate_curve(curve, status) result(forward)
      type(spot_rate_curve_t), intent(in) :: curve
      integer, intent(out), optional :: status
      type(forward_rate_t) :: forward
      forward = forwardrate_from_curve(curve, status)
   end function as_forwardrate_curve

   function as_forwardrate_spot(spot, interval, reference_date, status) result(forward)
      type(spot_rate_t), intent(in) :: spot
      type(term_t), intent(in) :: interval
      integer, intent(in), optional :: reference_date
      integer, intent(out), optional :: status
      type(forward_rate_t) :: forward
      type(term_t) :: days
      integer :: stat_i
      days = todays(spot%daycount, interval, stat_i)
      if (spot%size() /= interval%size()) then
         allocate(forward%rate(0), forward%interval_days(0))
         if (present(status)) status = FI_SIZE_MISMATCH
         return
      end if
      allocate(forward%rate(spot%size()), forward%interval_days(spot%size()))
      forward%rate = spot%rate
      forward%interval_days = days%value
      forward%compounding = spot%compounding
      forward%daycount = spot%daycount
      forward%calendar = spot%calendar
      if (present(reference_date)) forward%reference_date = reference_date
      if (present(status)) status = FI_OK
   end function as_forwardrate_spot

   subroutine set_interpolation(curve, model, status)
      type(spot_rate_curve_t), intent(inout) :: curve
      type(interpolation_t), intent(in) :: model
      integer, intent(out), optional :: status
      if (curve%size() < 2 .and. model%method < INTERP_NELSON_SIEGEL) then
         if (present(status)) status = FI_INVALID_ARGUMENT
         return
      end if
      curve%interpolation = model
      if (present(status)) status = FI_OK
   end subroutine set_interpolation

   function interpolate(curve, query_days, status) result(rate)
      type(spot_rate_curve_t), intent(in) :: curve
      real(dp), intent(in) :: query_days(:)
      integer, intent(out), optional :: status
      real(dp), allocatable :: rate(:)
      real(dp), allocatable :: factors(:), log_factors(:), log_query(:), query_years(:)
      integer :: stat_i, i
      select case (curve%interpolation%method)
      case (INTERP_FLAT_FORWARD)
         factors = curve_compound(curve, stat_i)
         if (stat_i /= FI_OK .or. any(factors <= 0.0_dp)) then
            allocate(rate(size(query_days)))
            rate = ieee_value(0.0_dp, ieee_quiet_nan)
            if (present(status)) status = FI_INVALID_ARGUMENT
            return
         end if
         log_factors = log(factors)
         log_query = interpolate_points(INTERPOLATION_LINEAR, curve%term_days, log_factors, query_days, status=stat_i)
         allocate(rate(size(query_days)))
         allocate(query_years(size(query_days)))
         query_years = query_days / real(curve%daycount%days_in_base, dp)
         do i = 1, size(query_days)
            if (ieee_is_nan(log_query(i)) .or. abs(query_years(i)) <= tiny(1.0_dp)) then
               rate(i) = ieee_value(0.0_dp, ieee_quiet_nan)
            else
               rate(i) = implied_rate(curve%compounding, query_years(i), exp(log_query(i)))
            end if
         end do
      case (INTERPOLATION_LINEAR, INTERP_LOG_LINEAR, INTERP_NATURAL_SPLINE, &
            INTERP_HERMITE_SPLINE, INTERP_MONOTONE_SPLINE)
         rate = interpolate_points(curve%interpolation%method, curve%term_days, curve%rate, &
                                   query_days, status=stat_i)
      case (INTERP_NELSON_SIEGEL, INTERP_NELSON_SIEGEL_SVENSSON)
         allocate(query_years(size(query_days)))
         query_years = query_days / real(curve%daycount%days_in_base, dp)
         rate = interpolate_points(curve%interpolation%method, [0.0_dp, 1.0_dp], [0.0_dp, 0.0_dp], &
                                   query_years, curve%interpolation%parameters, stat_i)
      case default
         allocate(rate(size(query_days)))
         rate = ieee_value(0.0_dp, ieee_quiet_nan)
         stat_i = FI_NOT_CONFIGURED
      end select
      if (present(status)) status = stat_i
   end function interpolate

   function fit_interpolation(initial_model, curve, max_iterations, tolerance) result(fit)
      type(interpolation_t), intent(in) :: initial_model
      type(spot_rate_curve_t), intent(in) :: curve
      integer, intent(in), optional :: max_iterations
      real(dp), intent(in), optional :: tolerance
      type(fit_result_t) :: fit
      real(dp), allocatable :: terms_years(:)
      terms_years = curve%term_days / real(curve%daycount%days_in_base, dp)
      fit = fit_interpolation_model(initial_model, terms_years, curve%rate, max_iterations, tolerance)
   end function fit_interpolation

   function forwardrate_from_curve(curve, status) result(forward)
      type(spot_rate_curve_t), intent(in) :: curve
      integer, intent(out), optional :: status
      type(forward_rate_t) :: forward
      real(dp), allocatable :: factor(:), ratio(:), years(:)
      integer :: i, stat_i
      if (curve%size() == 0) then
         allocate(forward%rate(0), forward%interval_days(0))
         if (present(status)) status = FI_INVALID_ARGUMENT
         return
      end if
      factor = curve_compound(curve, stat_i)
      allocate(forward%rate(curve%size()), forward%interval_days(curve%size()))
      forward%rate(1) = curve%rate(1)
      forward%interval_days(1) = curve%term_days(1)
      if (curve%size() > 1) then
         allocate(ratio(curve%size()-1), years(curve%size()-1))
         ratio = factor(2:) / factor(:curve%size()-1)
         forward%interval_days(2:) = curve%term_days(2:) - curve%term_days(:curve%size()-1)
         years = forward%interval_days(2:) / real(curve%daycount%days_in_base, dp)
         do i = 1, size(ratio)
            forward%rate(i+1) = implied_rate(curve%compounding, years(i), ratio(i), stat_i)
         end do
      end if
      forward%compounding = curve%compounding
      forward%daycount = curve%daycount
      forward%calendar = curve%calendar
      forward%reference_date = curve%reference_date
      if (present(status)) status = FI_OK
   end function forwardrate_from_curve

   function forwardrate_between(curve, t1_days, t2_days, status) result(forward)
      type(spot_rate_curve_t), intent(in) :: curve
      real(dp), intent(in) :: t1_days, t2_days
      integer, intent(out), optional :: status
      type(forward_rate_t) :: forward
      real(dp), allocatable :: rates(:), factors(:)
      real(dp) :: f1, f2, years, rate_value
      integer :: stat_i
      if (t2_days <= t1_days) then
         allocate(forward%rate(0), forward%interval_days(0))
         if (present(status)) status = FI_INVALID_ARGUMENT
         return
      end if
      if (curve%interpolation%method /= INTERP_NONE) then
         rates = interpolate(curve, [t1_days, t2_days], stat_i)
         factors = compound(curve%compounding, [t1_days, t2_days] / &
            real(curve%daycount%days_in_base, dp), rates, stat_i)
         f1 = factors(1)
         f2 = factors(2)
      else
         f1 = factor_at_exact_term(curve, t1_days, stat_i)
         if (stat_i /= FI_OK) then
            allocate(forward%rate(0), forward%interval_days(0))
            if (present(status)) status = stat_i
            return
         end if
         f2 = factor_at_exact_term(curve, t2_days, stat_i)
      end if
      years = (t2_days - t1_days) / real(curve%daycount%days_in_base, dp)
      rate_value = implied_rate(curve%compounding, years, f2/f1, stat_i)
      allocate(forward%rate(1), forward%interval_days(1))
      forward%rate = rate_value
      forward%interval_days = t2_days - t1_days
      forward%compounding = curve%compounding
      forward%daycount = curve%daycount
      forward%calendar = curve%calendar
      forward%reference_date = curve%reference_date
      if (present(status)) status = stat_i
   end function forwardrate_between

   function as_spotratecurve(forward, reference_date, status) result(curve)
      type(forward_rate_t), intent(in) :: forward
      integer, intent(in), optional :: reference_date
      integer, intent(out), optional :: status
      type(spot_rate_curve_t) :: curve
      real(dp), allocatable :: factors(:), cumulative_factor(:), cumulative_days(:), years(:)
      integer :: i, stat_i
      factors = forward_compound(forward, stat_i)
      allocate(cumulative_factor(forward%size()), cumulative_days(forward%size()), years(forward%size()))
      if (forward%size() > 0) then
         cumulative_factor(1) = factors(1)
         cumulative_days(1) = forward%interval_days(1)
         do i = 2, forward%size()
            cumulative_factor(i) = cumulative_factor(i-1) * factors(i)
            cumulative_days(i) = cumulative_days(i-1) + forward%interval_days(i)
         end do
      end if
      years = cumulative_days / real(forward%daycount%days_in_base, dp)
      allocate(curve%rate(forward%size()), curve%term_days(forward%size()))
      do i = 1, forward%size()
         curve%rate(i) = implied_rate(forward%compounding, years(i), cumulative_factor(i), stat_i)
      end do
      curve%term_days = cumulative_days
      curve%compounding = forward%compounding
      curve%daycount = forward%daycount
      curve%calendar = forward%calendar
      curve%reference_date = forward%reference_date
      if (present(reference_date)) curve%reference_date = reference_date
      if (present(status)) status = FI_OK
   end function as_spotratecurve

   function interpolation_error(curve, status) result(rmse)
      type(spot_rate_curve_t), intent(in) :: curve
      integer, intent(out), optional :: status
      real(dp) :: rmse
      real(dp), allocatable :: fitted(:)
      integer :: stat_i
      fitted = interpolate(curve, curve%term_days, stat_i)
      if (stat_i /= FI_OK .and. stat_i /= FI_OUT_OF_RANGE) then
         rmse = ieee_value(0.0_dp, ieee_quiet_nan)
      else
         rmse = sqrt(sum((fitted - curve%rate)**2) / real(curve%size(), dp))
         stat_i = FI_OK
      end if
      if (present(status)) status = stat_i
   end function interpolation_error

   function maturities(curve, status) result(dates)
      type(spot_rate_curve_t), intent(in) :: curve
      integer, intent(out), optional :: status
      integer, allocatable :: dates(:)
      integer :: i, stat_i
      allocate(dates(curve%size()))
      do i = 1, curve%size()
         dates(i) = offset_date(curve%reference_date, nint(curve%term_days(i)), curve%calendar, stat_i)
         if (stat_i /= FI_OK) then
            if (present(status)) status = FI_UNSUPPORTED_CALENDAR
            return
         end if
      end do
      if (present(status)) status = FI_OK
   end function maturities

   function first(curve, t, status) result(subset)
      type(spot_rate_curve_t), intent(in) :: curve
      type(term_t), intent(in) :: t
      integer, intent(out), optional :: status
      type(spot_rate_curve_t) :: subset
      type(term_t) :: days
      integer :: stat_i, n
      days = todays(curve%daycount, t, stat_i)
      if (days%size() /= 1) then
         allocate(subset%rate(0), subset%term_days(0))
         if (present(status)) status = FI_INVALID_ARGUMENT
         return
      end if
      n = count(curve%term_days <= days%value(1))
      subset = copy_curve_slice(curve, 1, n)
      if (present(status)) status = FI_OK
   end function first

   function last(curve, t, status) result(subset)
      type(spot_rate_curve_t), intent(in) :: curve
      type(term_t), intent(in) :: t
      integer, intent(out), optional :: status
      type(spot_rate_curve_t) :: subset
      type(term_t) :: days
      real(dp) :: cutoff
      integer :: stat_i, start
      days = todays(curve%daycount, t, stat_i)
      if (days%size() /= 1 .or. curve%size() == 0) then
         allocate(subset%rate(0), subset%term_days(0))
         if (present(status)) status = FI_INVALID_ARGUMENT
         return
      end if
      cutoff = curve%term_days(curve%size()) - days%value(1)
      start = 1
      do while (start <= curve%size() .and. curve%term_days(start) < cutoff)
         start = start + 1
      end do
      subset = copy_curve_slice(curve, start, curve%size())
      if (present(status)) status = FI_OK
   end function last

   function closest(curve, t, status) result(subset)
      type(spot_rate_curve_t), intent(in) :: curve
      type(term_t), intent(in) :: t
      integer, intent(out), optional :: status
      type(spot_rate_curve_t) :: subset
      type(term_t) :: days
      integer :: stat_i, index_min
      days = todays(curve%daycount, t, stat_i)
      if (days%size() /= 1 .or. curve%size() == 0) then
         allocate(subset%rate(0), subset%term_days(0))
         if (present(status)) status = FI_INVALID_ARGUMENT
         return
      end if
      index_min = minloc(abs(curve%term_days - days%value(1)), dim=1)
      subset = copy_curve_slice(curve, index_min, index_min)
      if (present(status)) status = FI_OK
   end function closest

   function curve_at_terms(curve, query_days, status) result(subset)
      type(spot_rate_curve_t), intent(in) :: curve
      real(dp), intent(in) :: query_days(:)
      integer, intent(out), optional :: status
      type(spot_rate_curve_t) :: subset
      real(dp), allocatable :: rates(:)
      integer :: i, idx, stat_i
      allocate(subset%rate(size(query_days)), subset%term_days(size(query_days)))
      subset%term_days = query_days
      if (curve%interpolation%method /= INTERP_NONE) then
         rates = interpolate(curve, query_days, stat_i)
         subset%rate = rates
      else
         do i = 1, size(query_days)
            idx = exact_index(curve%term_days, query_days(i))
            if (idx == 0) then
               subset%rate(i) = ieee_value(0.0_dp, ieee_quiet_nan)
            else
               subset%rate(i) = curve%rate(idx)
            end if
         end do
         stat_i = FI_OK
      end if
      subset%compounding = curve%compounding
      subset%daycount = curve%daycount
      subset%calendar = curve%calendar
      subset%reference_date = curve%reference_date
      if (curve%interpolation%propagate) subset%interpolation = curve%interpolation
      if (present(status)) status = stat_i
   end function curve_at_terms

   subroutine insert_curve_points(curve, term_days, rates, status)
      type(spot_rate_curve_t), intent(inout) :: curve
      real(dp), intent(in) :: term_days(:), rates(:)
      integer, intent(out), optional :: status
      real(dp), allocatable :: all_terms(:), all_rates(:), unique_terms(:), unique_rates(:)
      integer, allocatable :: order(:)
      integer :: i, j, n, idx
      if (size(term_days) /= size(rates)) then
         if (present(status)) status = FI_SIZE_MISMATCH
         return
      end if
      n = curve%size() + size(rates)
      allocate(all_terms(n), all_rates(n), order(n))
      all_terms = [curve%term_days, term_days]
      all_rates = [curve%rate, rates]
      order = [(i, i=1,n)]
      call sort_indices(all_terms, order)
      allocate(unique_terms(n), unique_rates(n))
      j = 0
      do i = 1, n
         idx = order(i)
         if (j == 0) then
            j = 1
            unique_terms(j) = all_terms(idx)
            unique_rates(j) = all_rates(idx)
         else if (abs(all_terms(idx)-unique_terms(j)) > 32.0_dp*epsilon(1.0_dp)* &
            max(1.0_dp,abs(all_terms(idx)))) then
            j = j + 1
            unique_terms(j) = all_terms(idx)
            unique_rates(j) = all_rates(idx)
         else
            unique_rates(j) = all_rates(idx)
         end if
      end do
      curve%term_days = unique_terms(:j)
      curve%rate = unique_rates(:j)
      if (present(status)) status = FI_OK
   end subroutine insert_curve_points

   function broadcast_compound(method, times, rates, status) result(factors)
      integer, intent(in) :: method
      real(dp), intent(in) :: times(:), rates(:)
      integer, intent(out) :: status
      real(dp), allocatable :: factors(:), use_times(:), use_rates(:)
      if (size(times) == size(rates)) then
         factors = compound(method, times, rates, status)
      else if (size(times) == 1) then
         allocate(use_times(size(rates)))
         use_times = times(1)
         factors = compound(method, use_times, rates, status)
      else if (size(rates) == 1) then
         allocate(use_rates(size(times)))
         use_rates = rates(1)
         factors = compound(method, times, use_rates, status)
      else
         allocate(factors(0))
         status = FI_SIZE_MISMATCH
      end if
   end function broadcast_compound

   real(dp) function factor_at_exact_term(curve, term_days, status) result(factor)
      type(spot_rate_curve_t), intent(in) :: curve
      real(dp), intent(in) :: term_days
      integer, intent(out) :: status
      integer :: idx
      idx = exact_index(curve%term_days, term_days)
      if (idx == 0) then
         factor = ieee_value(0.0_dp, ieee_quiet_nan)
         status = FI_OUT_OF_RANGE
      else
         factor = compound(curve%compounding, term_days / real(curve%daycount%days_in_base, dp), &
                           curve%rate(idx), status)
      end if
   end function factor_at_exact_term

   pure integer function exact_index(x, value) result(idx)
      real(dp), intent(in) :: x(:), value
      integer :: i
      idx = 0
      do i = 1, size(x)
         if (abs(x(i)-value) <= 32.0_dp*epsilon(1.0_dp)*max(1.0_dp,abs(value))) then
            idx = i
            return
         end if
      end do
   end function exact_index

   function copy_curve_slice(curve, first_index, last_index) result(subset)
      type(spot_rate_curve_t), intent(in) :: curve
      integer, intent(in) :: first_index, last_index
      type(spot_rate_curve_t) :: subset
      integer :: n
      n = max(0, last_index-first_index+1)
      allocate(subset%rate(n), subset%term_days(n))
      if (n > 0) then
         subset%rate = curve%rate(first_index:last_index)
         subset%term_days = curve%term_days(first_index:last_index)
      end if
      subset%compounding = curve%compounding
      subset%daycount = curve%daycount
      subset%calendar = curve%calendar
      subset%reference_date = curve%reference_date
      if (curve%interpolation%propagate) subset%interpolation = curve%interpolation
   end function copy_curve_slice

   subroutine sort_indices(values, order)
      real(dp), intent(in) :: values(:)
      integer, intent(inout) :: order(:)
      integer :: i, j, key
      do i = 2, size(order)
         key = order(i)
         j = i - 1
         do while (j >= 1)
            if (values(order(j)) <= values(key)) exit
            order(j+1) = order(j)
            j = j - 1
         end do
         order(j+1) = key
      end do
   end subroutine sort_indices

end module fixedincome_curves
