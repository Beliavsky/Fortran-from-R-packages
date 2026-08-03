! SPDX-License-Identifier: GPL-2.0-or-later
module ragtop_term_structures
   use ragtop_kinds, only : dp
   use ragtop_constants, only : ragtop_ok, ragtop_invalid_argument
   use ragtop_types, only : discount_curve, volatility_curve, market_spec
   implicit none
   private
   public :: initialize_discount_curve, initialize_volatility_curve
   public :: discount_factor, cumulative_variance, survival_probability
   public :: default_intensity, effective_rate, effective_volatility, effective_hazard

contains

   subroutine initialize_discount_curve(curve, time, rate, status)
      type(discount_curve), intent(out) :: curve
      real(dp), intent(in) :: time(:), rate(:)
      integer, intent(out), optional :: status
      integer :: n, i, st
      st = ragtop_ok
      n = size(time)
      if (n < 1 .or. size(rate) /= n .or. any(time < 0.0_dp)) then
         st = ragtop_invalid_argument
         if (present(status)) status = st
         return
      end if
      do i = 2, n
         if (time(i) <= time(i-1)) st = ragtop_invalid_argument
      end do
      if (st /= ragtop_ok) then
         if (present(status)) status = st
         return
      end if
      allocate(curve%time(n), curve%rate(n), curve%df(n), curve%forward_rate(n))
      curve%time = time
      curve%rate = rate
      curve%df = exp(-time*rate)
      if (n == 1) then
         curve%forward_rate(1) = rate(1)
      else
         do i = 1, n-1
            curve%forward_rate(i) = -log(curve%df(i+1)/curve%df(i))/(time(i+1)-time(i))
         end do
         curve%forward_rate(n) = curve%forward_rate(n-1)
      end if
      if (present(status)) status = st
   end subroutine initialize_discount_curve

   subroutine initialize_volatility_curve(curve, time, volatility, status)
      type(volatility_curve), intent(out) :: curve
      real(dp), intent(in) :: time(:), volatility(:)
      integer, intent(out), optional :: status
      real(dp), allocatable :: aug_time(:), aug_var(:)
      real(dp) :: dv, dt
      integer :: n, i, st
      st = ragtop_ok
      n = size(time)
      if (n < 1 .or. size(volatility) /= n .or. any(time <= 0.0_dp) .or. any(volatility < 0.0_dp)) then
         st = ragtop_invalid_argument
         if (present(status)) status = st
         return
      end if
      do i = 2, n
         if (time(i) <= time(i-1)) st = ragtop_invalid_argument
      end do
      if (st /= ragtop_ok) then
         if (present(status)) status = st
         return
      end if
      allocate(curve%time(n), curve%volatility(n), curve%cumulative_variance(n), curve%forward_volatility(n))
      curve%time = time
      curve%volatility = volatility
      curve%cumulative_variance = time*volatility*volatility
      allocate(aug_time(0:n), aug_var(0:n))
      aug_time(0) = 0.0_dp
      aug_var(0) = 0.0_dp
      aug_time(1:n) = time
      aug_var(1:n) = curve%cumulative_variance
      do i = 1, n
         dv = aug_var(i)-aug_var(i-1)
         dt = aug_time(i)-aug_time(i-1)
         if (dv < -100.0_dp*epsilon(1.0_dp)) then
            st = ragtop_invalid_argument
            curve%forward_volatility = 0.0_dp
            if (present(status)) status = st
            return
         end if
         curve%forward_volatility(i) = sqrt(max(0.0_dp,dv/dt))
      end do
      if (present(status)) status = st
   end subroutine initialize_volatility_curve

   pure real(dp) function curve_df0(curve, x) result(df)
      type(discount_curve), intent(in) :: curve
      real(dp), intent(in) :: x
      integer :: i, n
      n = size(curve%time)
      if (x <= 0.0_dp) then
         df = 1.0_dp
      else if (x < curve%time(1)) then
         df = exp(-curve%rate(1)*x)
      else
         i = n
         do while (i > 1)
            if (x >= curve%time(i)) exit
            i = i-1
         end do
         df = curve%df(i)*exp(-curve%forward_rate(i)*(x-curve%time(i)))
      end if
   end function curve_df0

   pure real(dp) function discount_factor(market, t_maturity, t_model) result(df)
      type(market_spec), intent(in) :: market
      real(dp), intent(in) :: t_maturity
      real(dp), intent(in), optional :: t_model
      real(dp) :: t
      t = 0.0_dp
      if (present(t_model)) t = t_model
      if (market%use_rate_curve .and. allocated(market%rates%time)) then
         df = curve_df0(market%rates,t_maturity)/curve_df0(market%rates,t)
      else
         df = exp(-market%short_rate*(t_maturity-t))
      end if
   end function discount_factor

   pure real(dp) function curve_cv0(curve, x) result(cv)
      type(volatility_curve), intent(in) :: curve
      real(dp), intent(in) :: x
      integer :: i, n
      real(dp) :: anchor_t, anchor_var
      n = size(curve%time)
      if (x <= 0.0_dp) then
         cv = 0.0_dp
      else if (x < curve%time(1)) then
         cv = curve%forward_volatility(1)**2*x
      else
         i = n
         do while (i > 1)
            if (x >= curve%time(i)) exit
            i = i-1
         end do
         anchor_t = curve%time(i)
         anchor_var = curve%cumulative_variance(i)
         cv = anchor_var + curve%forward_volatility(i)**2*(x-anchor_t)
      end if
   end function curve_cv0

   pure real(dp) function cumulative_variance(market, t_maturity, t_model) result(cv)
      type(market_spec), intent(in) :: market
      real(dp), intent(in) :: t_maturity
      real(dp), intent(in), optional :: t_model
      real(dp) :: t
      t = 0.0_dp
      if (present(t_model)) t = t_model
      if (market%use_vol_curve .and. allocated(market%vols%time)) then
         cv = curve_cv0(market%vols,t_maturity)-curve_cv0(market%vols,t)
      else
         cv = market%volatility**2*(t_maturity-t)
      end if
   end function cumulative_variance

   elemental real(dp) function linked_hazard(base, constant_fraction, power, reference_spot, spot) result(h)
      real(dp), intent(in) :: base, constant_fraction, power, reference_spot, spot
      real(dp) :: safe_spot
      safe_spot = max(spot, sqrt(tiny(1.0_dp)))
      h = base*(constant_fraction+(1.0_dp-constant_fraction)*(reference_spot/safe_spot)**power)
      h = max(0.0_dp,h)
   end function linked_hazard

   elemental real(dp) function default_intensity(market, t, spot) result(h)
      type(market_spec), intent(in) :: market
      real(dp), intent(in) :: t, spot
      if (market%use_hazard_link) then
         h = 0.0_dp*t + linked_hazard(market%hazard%base_intensity, market%hazard%constant_fraction, &
                           market%hazard%power, market%hazard%reference_spot, spot)
      else
         h = max(0.0_dp,market%default_intensity)
      end if
   end function default_intensity

   pure real(dp) function survival_probability(market, t_maturity, t_model, spot) result(p)
      type(market_spec), intent(in) :: market
      real(dp), intent(in) :: t_maturity
      real(dp), intent(in), optional :: t_model, spot
      real(dp) :: t, s, h
      t = 0.0_dp
      if (present(t_model)) t = t_model
      s = market%hazard%reference_spot
      if (present(spot)) s = spot
      h = default_intensity(market,t,s)
      p = exp(-h*(t_maturity-t))
   end function survival_probability

   pure real(dp) function effective_rate(market, maturity) result(r)
      type(market_spec), intent(in) :: market
      real(dp), intent(in) :: maturity
      if (maturity <= 0.0_dp) then
         r = 0.0_dp
      else
         r = -log(discount_factor(market,maturity,0.0_dp))/maturity
      end if
   end function effective_rate

   pure real(dp) function effective_volatility(market, maturity) result(v)
      type(market_spec), intent(in) :: market
      real(dp), intent(in) :: maturity
      if (maturity <= 0.0_dp) then
         v = 0.0_dp
      else
         v = sqrt(max(0.0_dp,cumulative_variance(market,maturity,0.0_dp)/maturity))
      end if
   end function effective_volatility

   pure real(dp) function effective_hazard(market, spot) result(h)
      type(market_spec), intent(in) :: market
      real(dp), intent(in) :: spot
      h = default_intensity(market,0.0_dp,spot)
   end function effective_hazard

end module ragtop_term_structures
