! SPDX-License-Identifier: GPL-2.0-or-later
module ragtop_cashflows
   use ragtop_kinds, only : dp
   use ragtop_types, only : dividend_schedule, cashflow_schedule, market_spec
   use ragtop_term_structures, only : discount_factor
   use ragtop_math, only : cubic_spline_interp
   implicit none
   private
   public :: time_adjusted_dividend_one, time_adjusted_dividends, shift_for_dividends, adjust_for_dividends
   public :: value_from_prior_coupons, accelerated_coupon_value, coupon_value_at_exercise
   public :: cashflows_between

contains

   pure real(dp) function time_adjusted_dividend_one(div_time, fixed, proportional, t_final, r, h, spot, spot0) result(v)
      real(dp), intent(in) :: div_time, fixed, proportional, t_final, r, h, spot, spot0
      real(dp) :: dt
      dt = div_time-t_final
      v = exp(-(r+h)*dt)*(fixed+proportional*spot/max(spot0,sqrt(tiny(1.0_dp))))
   end function time_adjusted_dividend_one

   subroutine time_adjusted_dividends(dividends, t_final, r, h, spot, spot0, values)
      type(dividend_schedule), intent(in) :: dividends
      real(dp), intent(in) :: t_final, r, h, spot(:), spot0
      real(dp), intent(out) :: values(:)
      integer :: i
      values = 0.0_dp
      if (.not. allocated(dividends%time)) return
      do i = 1, size(dividends%time)
         values = values + exp(-(r+h)*(dividends%time(i)-t_final))* &
            (dividends%fixed(i)+dividends%proportional(i)*spot/max(spot0,sqrt(tiny(1.0_dp))))
      end do
   end subroutine time_adjusted_dividends

   subroutine shift_for_dividends(values_before, stock, dividend_amount, values_after)
      real(dp), intent(in) :: values_before(:), stock(:), dividend_amount(:)
      real(dp), intent(out) :: values_after(:)
      integer :: i, n
      real(dp) :: shifted
      n = size(stock)
      do i = 1, n
         shifted = max(0.0_dp,stock(i)-dividend_amount(i))
         values_after(i) = cubic_spline_interp(stock,values_before,shifted)
      end do
   end subroutine shift_for_dividends

   subroutine adjust_for_dividends(values, t, dt, r, h, stock, spot0, dividends, adjusted)
      real(dp), intent(in) :: values(:,:), t, dt, r, h(:), stock(:), spot0
      type(dividend_schedule), intent(in) :: dividends
      real(dp), intent(out) :: adjusted(:,:)
      type(dividend_schedule) :: relevant
      real(dp), allocatable :: div_sum(:)
      integer :: i, j, nrel
      logical, allocatable :: mask(:)
      adjusted = values
      if (.not. allocated(dividends%time)) return
      allocate(mask(size(dividends%time)))
      mask = dividends%time > t .and. dividends%time <= t+dt+100.0_dp*epsilon(1.0_dp)
      nrel = count(mask)
      if (nrel == 0) return
      allocate(relevant%time(nrel), relevant%fixed(nrel), relevant%proportional(nrel))
      j = 0
      do i = 1, size(mask)
         if (mask(i)) then
            j = j+1
            relevant%time(j) = dividends%time(i)
            relevant%fixed(j) = dividends%fixed(i)
            relevant%proportional(j) = dividends%proportional(i)
         end if
      end do
      allocate(div_sum(size(stock)))
      do j = 1, size(values,2)
         div_sum = 0.0_dp
         do i = 1, nrel
            div_sum = div_sum + exp(-(r+h)*(relevant%time(i)-t))* &
               (relevant%fixed(i)+relevant%proportional(i)*stock/max(spot0,sqrt(tiny(1.0_dp))))
         end do
         call shift_for_dividends(values(:,j),stock,div_sum,adjusted(:,j))
      end do
   end subroutine adjust_for_dividends

   pure real(dp) function value_from_prior_coupons(t, coupons, market, model_t) result(v)
      real(dp), intent(in) :: t
      type(cashflow_schedule), intent(in) :: coupons
      type(market_spec), intent(in) :: market
      real(dp), intent(in), optional :: model_t
      real(dp) :: mt
      integer :: i
      mt = 0.0_dp
      if (present(model_t)) mt = model_t
      v = 0.0_dp
      if (.not. allocated(coupons%time)) return
      do i = 1, size(coupons%time)
         if (coupons%time(i) <= t .and. coupons%time(i) > mt) then
            v = v + coupons%amount(i)* &
                discount_factor(market,coupons%time(i),t)
         end if
      end do
   end function value_from_prior_coupons

   pure real(dp) function accelerated_coupon_value(t, coupons, market, acceleration_t) result(v)
      real(dp), intent(in) :: t
      type(cashflow_schedule), intent(in) :: coupons
      type(market_spec), intent(in) :: market
      real(dp), intent(in), optional :: acceleration_t
      real(dp) :: at
      integer :: i
      at = huge(1.0_dp)
      if (present(acceleration_t)) at = acceleration_t
      v = 0.0_dp
      if (.not. allocated(coupons%time)) return
      do i = 1, size(coupons%time)
         if (coupons%time(i) > t .and. coupons%time(i) <= at) then
            v = v + coupons%amount(i)*discount_factor(market,coupons%time(i),t)
         end if
      end do
   end function accelerated_coupon_value

   pure real(dp) function coupon_value_at_exercise(t, coupons, market, model_t, accelerate_future, acceleration_t) result(v)
      real(dp), intent(in) :: t
      type(cashflow_schedule), intent(in) :: coupons
      type(market_spec), intent(in) :: market
      real(dp), intent(in), optional :: model_t, acceleration_t
      logical, intent(in), optional :: accelerate_future
      logical :: accel
      accel = .false.
      if (present(accelerate_future)) accel = accelerate_future
      v = value_from_prior_coupons(t,coupons,market,model_t)
      if (accel) v = v + accelerated_coupon_value(t,coupons,market,acceleration_t)
   end function coupon_value_at_exercise

   pure real(dp) function cashflows_between(coupons, t0, t1, market, value_time) result(v)
      type(cashflow_schedule), intent(in) :: coupons
      real(dp), intent(in) :: t0, t1
      type(market_spec), intent(in) :: market
      real(dp), intent(in), optional :: value_time
      real(dp) :: vt
      integer :: i
      vt = t1
      if (present(value_time)) vt = value_time
      v = 0.0_dp
      if (.not. allocated(coupons%time)) return
      do i = 1, size(coupons%time)
         if (coupons%time(i) > t0 .and. coupons%time(i) <= t1+100.0_dp*epsilon(1.0_dp)) then
            v = v + coupons%amount(i)*discount_factor(market,coupons%time(i),vt)
         end if
      end do
   end function cashflows_between

end module ragtop_cashflows
