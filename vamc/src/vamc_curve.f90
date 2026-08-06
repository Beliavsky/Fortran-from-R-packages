module vamc_curve
  use vamc_kinds, only: dp
  use vamc_status, only: status_type, vamc_invalid_argument, vamc_dimension_error, vamc_numerical_error
  use vamc_dates, only: date_type, add_days, add_years, frac_year, gen_schedule, roll_date, &
                        bdc_following, dcc_from_string, bdc_from_string, calendar_from_string, &
                        operator(==), operator(<), operator(>)
  implicit none
  private

  type, public :: yield_curve_type
    type(date_type) :: curve_date
    type(date_type) :: settle_date
    type(date_type), allocatable :: observation_dates(:)
    real(dp), allocatable :: discount_factors(:)
    real(dp), allocatable :: zero_rates(:)
    real(dp), allocatable :: forward_rates(:)
    real(dp), allocatable :: day_counts(:)
    integer :: day_count_convention = 0
  contains
    procedure :: discount => curve_discount
    procedure :: zero_rate => curve_zero_rate
  end type yield_curve_type

  public :: build_curve, build_curve_named, log_linear_discount, present_value_swap
contains
  real(dp) function log_linear_discount(date_in, observation_dates, discount_factors, curve_date, day_count_convention, status)
    type(date_type), intent(in) :: date_in, observation_dates(:), curve_date
    real(dp), intent(in) :: discount_factors(:)
    integer, intent(in) :: day_count_convention
    type(status_type), intent(inout), optional :: status
    integer :: n, i, lo, hi
    real(dp) :: t1, t2, tin, temp
    if (present(status)) call status%clear()
    n = size(observation_dates)
    if (n < 2 .or. size(discount_factors) /= n) then
      log_linear_discount = 0.0_dp
      if (present(status)) call status%set(vamc_dimension_error, 'At least two conformable curve observations are required.')
      return
    end if
    if (date_in < observation_dates(1)) then
      log_linear_discount = 0.0_dp
      if (present(status)) call status%set(vamc_invalid_argument, 'Date is earlier than the curve date.')
      return
    end if
    do i = 1, n
      if (date_in == observation_dates(i)) then
        log_linear_discount = discount_factors(i)
        return
      end if
    end do
    lo = n - 1
    hi = n
    do i = 1, n - 1
      if (date_in > observation_dates(i) .and. date_in < observation_dates(i+1)) then
        lo = i
        hi = i + 1
        exit
      end if
    end do
    if (discount_factors(lo) <= 0.0_dp .or. discount_factors(hi) <= 0.0_dp) then
      log_linear_discount = 0.0_dp
      if (present(status)) call status%set(vamc_numerical_error, 'Discount factors must be positive.')
      return
    end if
    t1 = frac_year(curve_date, observation_dates(lo), day_count_convention)
    t2 = frac_year(curve_date, observation_dates(hi), day_count_convention)
    tin = frac_year(curve_date, date_in, day_count_convention)
    if (abs(t2 - t1) <= epsilon(1.0_dp)) then
      log_linear_discount = discount_factors(lo)
      return
    end if
    temp = log(discount_factors(hi)) * (tin - t1) + log(discount_factors(lo)) * (t2 - tin)
    log_linear_discount = exp(temp / (t2 - t1))
  end function log_linear_discount

  real(dp) function present_value_swap(rate, tenor_years, fixed_frequency, fixed_dcc, floating_frequency, floating_dcc, &
                                        calendar, bdc, curve_date, yield_curve_dcc, payment_dates, discount_factors, &
                                        settle_date, holidays, source_compatible_holidays, status)
    real(dp), intent(in) :: rate
    integer, intent(in) :: tenor_years, fixed_frequency, fixed_dcc, floating_frequency, floating_dcc
    integer, intent(in) :: calendar, bdc, yield_curve_dcc
    type(date_type), intent(in) :: curve_date, payment_dates(:), settle_date
    real(dp), intent(in) :: discount_factors(:)
    type(date_type), intent(in), optional :: holidays(:)
    logical, intent(in), optional :: source_compatible_holidays
    type(status_type), intent(inout), optional :: status
    type(date_type), allocatable :: fixed_schedule(:), floating_schedule(:)
    type(status_type) :: local_status
    real(dp) :: fixed_pv, floating_pv, df, dt, forward
    integer :: i
    if (present(status)) call status%clear()
    call gen_schedule(settle_date, fixed_frequency, tenor_years, calendar, bdc, fixed_schedule, holidays, &
                      source_compatible_holidays, local_status)
    if (.not. local_status%ok() .or. size(fixed_schedule) < 2) then
      present_value_swap = 0.0_dp
      if (present(status)) call status%set(vamc_invalid_argument, 'Invalid fixed-leg schedule.')
      return
    end if
    fixed_pv = 0.0_dp
    do i = 2, size(fixed_schedule)
      df = log_linear_discount(fixed_schedule(i), payment_dates, discount_factors, curve_date, yield_curve_dcc, local_status)
      if (.not. local_status%ok()) then
        present_value_swap = 0.0_dp
        if (present(status)) call status%set(local_status%code, local_status%message)
        return
      end if
      fixed_pv = fixed_pv + frac_year(fixed_schedule(i-1), fixed_schedule(i), fixed_dcc) * rate * df
    end do
    call gen_schedule(settle_date, floating_frequency, tenor_years, calendar, bdc, floating_schedule, holidays, &
                      source_compatible_holidays, local_status)
    if (.not. local_status%ok() .or. size(floating_schedule) < 2) then
      present_value_swap = 0.0_dp
      if (present(status)) call status%set(vamc_invalid_argument, 'Invalid floating-leg schedule.')
      return
    end if
    floating_pv = 0.0_dp
    do i = 2, size(floating_schedule)
      df = log_linear_discount(floating_schedule(i), payment_dates, discount_factors, curve_date, yield_curve_dcc, local_status)
      dt = frac_year(floating_schedule(i-1), floating_schedule(i), floating_dcc)
      if (dt <= 0.0_dp) cycle
      forward = (log_linear_discount(floating_schedule(i-1), payment_dates, discount_factors, curve_date, &
                                     yield_curve_dcc, local_status) / df - 1.0_dp) / dt
      floating_pv = floating_pv + dt * forward * df
    end do
    present_value_swap = fixed_pv - floating_pv
  end function present_value_swap

  subroutine build_curve(swap_rates, tenors, fixed_frequency, fixed_dcc, floating_frequency, floating_dcc, calendar, bdc, &
                         curve_date, settlement_days, yield_curve_dcc, curve, holidays, source_compatible_holidays, status)
    real(dp), intent(in) :: swap_rates(:)
    integer, intent(in) :: tenors(:), fixed_frequency, fixed_dcc, floating_frequency, floating_dcc, calendar, bdc
    type(date_type), intent(in) :: curve_date
    integer, intent(in) :: settlement_days, yield_curve_dcc
    type(yield_curve_type), intent(out) :: curve
    type(date_type), intent(in), optional :: holidays(:)
    logical, intent(in), optional :: source_compatible_holidays
    type(status_type), intent(inout), optional :: status
    integer :: n, i, iter
    real(dp) :: r0, r1, f0, f1, dx, f, denom
    type(status_type) :: local_status
    if (present(status)) call status%clear()
    n = size(swap_rates)
    if (n < 1 .or. size(tenors) /= n) then
      if (present(status)) call status%set(vamc_dimension_error, &
        'Number of swap rates differs from number of tenors.')
      return
    end if
    if (any(tenors <= 0) .or. any(swap_rates < 0.0_dp) .or. settlement_days < 0) then
      if (present(status)) call status%set(vamc_invalid_argument, 'Invalid curve inputs.')
      return
    end if
    curve%curve_date = curve_date
    curve%settle_date = curve_date
    do i = 1, settlement_days
      curve%settle_date = roll_date(add_days(curve%settle_date, 1), bdc_following, calendar, holidays, &
                                    source_compatible_holidays)
    end do
    curve%day_count_convention = yield_curve_dcc
    allocate(curve%observation_dates(n+1), curve%discount_factors(n+1), curve%zero_rates(n+1), &
             curve%forward_rates(n+1), curve%day_counts(n+1))
    curve%observation_dates(1) = curve_date
    curve%discount_factors = 1.0_dp
    do i = 1, n
      curve%observation_dates(i+1) = add_years(curve%settle_date, tenors(i))
    end do
    do i = 1, n
      r0 = exp(-swap_rates(i) * frac_year(curve_date, curve%observation_dates(i+1), yield_curve_dcc))
      curve%discount_factors(i+1) = r0
      f0 = present_value_swap(swap_rates(i), tenors(i), fixed_frequency, fixed_dcc, floating_frequency, floating_dcc, &
                              calendar, bdc, curve_date, yield_curve_dcc, curve%observation_dates(1:i+1), &
                              curve%discount_factors(1:i+1), curve%settle_date, holidays, source_compatible_holidays, local_status)
      if (.not. local_status%ok()) then
        if (present(status)) call status%set(local_status%code, local_status%message)
        return
      end if
      r1 = r0 + 0.001_dp
      curve%discount_factors(i+1) = r1
      f1 = present_value_swap(swap_rates(i), tenors(i), fixed_frequency, fixed_dcc, floating_frequency, floating_dcc, &
                              calendar, bdc, curve_date, yield_curve_dcc, curve%observation_dates(1:i+1), &
                              curve%discount_factors(1:i+1), curve%settle_date, holidays, source_compatible_holidays, local_status)
      do iter = 1, 1000
        if (abs(r1-r0) < 1.0e-10_dp) exit
        denom = f1 - f0
        if (abs(denom) < 1.0e-16_dp) then
          dx = 0.5_dp * (r0+r1)
        else
          dx = r0 - f0 * (r1-r0) / denom
        end if
        if (dx <= 0.0_dp .or. .not. (dx < huge(1.0_dp))) dx = max(1.0e-12_dp, 0.5_dp*(r0+r1))
        curve%discount_factors(i+1) = dx
        f = present_value_swap(swap_rates(i), tenors(i), fixed_frequency, fixed_dcc, floating_frequency, floating_dcc, &
                               calendar, bdc, curve_date, yield_curve_dcc, curve%observation_dates(1:i+1), &
                               curve%discount_factors(1:i+1), curve%settle_date, holidays, source_compatible_holidays, local_status)
        r0 = r1
        f0 = f1
        r1 = dx
        f1 = f
      end do
      curve%discount_factors(i+1) = r1
    end do
    do i = 1, n+1
      curve%day_counts(i) = frac_year(curve_date, curve%observation_dates(i), yield_curve_dcc)
      if (curve%day_counts(i) <= 0.0_dp) then
        curve%zero_rates(i) = 0.0_dp
      else
        curve%zero_rates(i) = -log(curve%discount_factors(i)) / curve%day_counts(i)
      end if
    end do
    do i = 1, n
      curve%forward_rates(i) = ( &
        log_linear_discount(curve%observation_dates(i), curve%observation_dates, curve%discount_factors, curve_date, &
                            yield_curve_dcc, local_status) / &
        log_linear_discount(curve%observation_dates(i+1), curve%observation_dates, curve%discount_factors, curve_date, &
                            yield_curve_dcc, local_status) - 1.0_dp) / &
        frac_year(curve%observation_dates(i), curve%observation_dates(i+1), yield_curve_dcc)
    end do
    curve%forward_rates(n+1) = curve%forward_rates(n)
  end subroutine build_curve

  subroutine build_curve_named(swap_rates, tenors, fixed_frequency, fixed_dcc_name, floating_frequency, floating_dcc_name, &
                               calendar_name, bdc_name, curve_date, settlement_days, yield_curve_dcc_name, curve, holidays, &
                               source_compatible_holidays, status)
    real(dp), intent(in) :: swap_rates(:)
    integer, intent(in) :: tenors(:), fixed_frequency, floating_frequency, settlement_days
    character(len=*), intent(in) :: fixed_dcc_name, floating_dcc_name, calendar_name, bdc_name, yield_curve_dcc_name
    type(date_type), intent(in) :: curve_date
    type(yield_curve_type), intent(out) :: curve
    type(date_type), intent(in), optional :: holidays(:)
    logical, intent(in), optional :: source_compatible_holidays
    type(status_type), intent(inout), optional :: status
    integer :: fd, fld, cal, bd, yd
    fd = dcc_from_string(fixed_dcc_name)
    fld = dcc_from_string(floating_dcc_name)
    cal = calendar_from_string(calendar_name)
    bd = bdc_from_string(bdc_name)
    yd = dcc_from_string(yield_curve_dcc_name)
    if (fd == 0 .or. fld == 0 .or. cal == 0 .or. bd == 0 .or. yd == 0) then
      if (present(status)) call status%set(vamc_invalid_argument, 'Invalid day-count, calendar, or business-day convention.')
      return
    end if
    call build_curve(swap_rates, tenors, fixed_frequency, fd, floating_frequency, fld, cal, bd, curve_date, &
                     settlement_days, yd, curve, holidays, source_compatible_holidays, status)
  end subroutine build_curve_named

  real(dp) function curve_discount(self, date, status)
    class(yield_curve_type), intent(in) :: self
    type(date_type), intent(in) :: date
    type(status_type), intent(inout), optional :: status
    curve_discount = log_linear_discount(date, self%observation_dates, self%discount_factors, self%curve_date, &
                                         self%day_count_convention, status)
  end function curve_discount

  real(dp) function curve_zero_rate(self, date, status)
    class(yield_curve_type), intent(in) :: self
    type(date_type), intent(in) :: date
    type(status_type), intent(inout), optional :: status
    real(dp) :: t, df
    t = frac_year(self%curve_date, date, self%day_count_convention)
    if (t <= 0.0_dp) then
      curve_zero_rate = 0.0_dp
    else
      df = self%discount(date, status)
      curve_zero_rate = -log(df) / t
    end if
  end function curve_zero_rate
end module vamc_curve
