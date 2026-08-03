! SPDX-License-Identifier: GPL-3.0-or-later
module qbc_curves
   use qbc_kinds, only : dp
   use qbc_status, only : qbc_success, qbc_invalid_argument, qbc_size_mismatch
   use qbc_dates, only : qbc_date, discount_time
   use qbc_types, only : qbc_curve, qbc_rate_continuous, qbc_rate_discrete
   implicit none
   private
   public :: curve_rate, interpolate_curve, spot_to_forward, forward_to_spot
   public :: discount_factor, discount_factors, curve_from_market_yields

contains

   pure real(dp) function curve_rate(curve, term) result(rate)
      type(qbc_curve), intent(in) :: curve
      real(dp), intent(in) :: term
      integer :: i, n
      real(dp) :: w
      n = size(curve%terms)
      if (n == 0) then
         rate = 0.0_dp
      else if (n == 1 .or. term <= curve%terms(1)) then
         rate = curve%rates(1)
      else if (term >= curve%terms(n)) then
         rate = curve%rates(n)
      else
         i = 1
         do while (i < n .and. term > curve%terms(i + 1))
            i = i + 1
         end do
         if (curve%approximation == 1) then
            rate = curve%rates(i)
         else
            w = (term - curve%terms(i)) / (curve%terms(i + 1) - curve%terms(i))
            rate = (1.0_dp - w) * curve%rates(i) + w * curve%rates(i + 1)
         end if
      end if
   end function curve_rate

   subroutine interpolate_curve(curve, nodes, values, status)
      type(qbc_curve), intent(in) :: curve
      real(dp), intent(in) :: nodes(:)
      real(dp), allocatable, intent(out) :: values(:)
      integer, intent(out), optional :: status
      integer :: i, st
      st = qbc_success
      if (.not. allocated(curve%terms) .or. .not. allocated(curve%rates) .or. &
          size(curve%terms) /= size(curve%rates)) then
         allocate(values(0)); st = qbc_invalid_argument
      else
         allocate(values(size(nodes)))
         do i = 1, size(nodes)
            values(i) = curve_rate(curve, nodes(i))
         end do
      end if
      if (present(status)) status = st
   end subroutine interpolate_curve

   subroutine spot_to_forward(terms, spot, approximation, forward, out_terms, status)
      real(dp), intent(in) :: terms(:), spot(:)
      integer, intent(in) :: approximation
      real(dp), allocatable, intent(out) :: forward(:), out_terms(:)
      integer, intent(out), optional :: status
      real(dp), allocatable :: t(:), s(:)
      integer :: i, n, st
      st = qbc_success
      if (size(terms) /= size(spot) .or. size(terms) == 0) then
         allocate(forward(0), out_terms(0)); st = qbc_size_mismatch
         if (present(status)) status = st
         return
      end if
      if (abs(terms(1)) > 100.0_dp * epsilon(1.0_dp)) then
         allocate(t(size(terms) + 1), s(size(spot) + 1))
         t = [0.0_dp, terms]
         s = [spot(1), spot]
      else
         allocate(t(size(terms)), s(size(spot)))
         t = terms; s = spot
      end if
      n = size(t)
      allocate(forward(n), out_terms(n))
      out_terms = t
      if (approximation == 1) then
         forward = s
      else if (approximation == 2) then
         forward(1) = s(1)
         do i = 2, n
            if (t(i) <= t(i-1)) then
               st = qbc_invalid_argument
               forward(i) = forward(i-1)
            else
               forward(i) = ((s(i) - s(i-1)) / (t(i) - t(i-1))) * t(i-1) + s(i)
            end if
         end do
      else
         st = qbc_invalid_argument
         forward = 0.0_dp
      end if
      if (present(status)) status = st
   end subroutine spot_to_forward

   subroutine forward_to_spot(terms, forward, approximation, spot, out_terms, status)
      real(dp), intent(in) :: terms(:), forward(:)
      integer, intent(in) :: approximation
      real(dp), allocatable, intent(out) :: spot(:), out_terms(:)
      integer, intent(out), optional :: status
      real(dp), allocatable :: t(:), f(:)
      real(dp) :: area, dt
      integer :: i, n, st
      st = qbc_success
      if (size(terms) /= size(forward) .or. size(terms) == 0) then
         allocate(spot(0), out_terms(0)); st = qbc_size_mismatch
         if (present(status)) status = st
         return
      end if
      if (abs(terms(1)) > 100.0_dp * epsilon(1.0_dp)) then
         allocate(t(size(terms) + 1), f(size(forward) + 1))
         t = [0.0_dp, terms]
         f = [forward(1), forward]
      else
         allocate(t(size(terms)), f(size(forward)))
         t = terms; f = forward
      end if
      n = size(t)
      allocate(spot(n), out_terms(n))
      out_terms = t
      spot(1) = f(1)
      area = 0.0_dp
      do i = 2, n
         dt = t(i) - t(i-1)
         if (dt <= 0.0_dp) then
            st = qbc_invalid_argument
         else if (approximation == 1) then
            area = area + f(i) * dt
         else if (approximation == 2) then
            area = area + 0.5_dp * (f(i-1) + f(i)) * dt
         else
            st = qbc_invalid_argument
         end if
         if (t(i) > 0.0_dp) spot(i) = area / t(i)
      end do
      if (present(status)) status = st
   end subroutine forward_to_spot

   pure real(dp) function discount_factor(rate, term, rate_type, frequency) result(df)
      real(dp), intent(in) :: rate, term
      integer, intent(in) :: rate_type
      integer, intent(in), optional :: frequency
      integer :: freq
      freq = 1
      if (present(frequency)) freq = max(1, frequency)
      if (rate_type == qbc_rate_continuous) then
         df = exp(-rate * term)
      else
         if (1.0_dp + rate / real(freq, dp) <= 0.0_dp) then
            df = huge(1.0_dp)
         else
            df = (1.0_dp + rate / real(freq, dp)) ** (-term * real(freq, dp))
         end if
      end if
   end function discount_factor

   subroutine discount_factors(dates, rates, analysis_date, rate_type, frequency, factors, status)
      type(qbc_date), intent(in) :: dates(:), analysis_date
      real(dp), intent(in) :: rates(:)
      integer, intent(in) :: rate_type
      integer, intent(in), optional :: frequency
      real(dp), allocatable, intent(out) :: factors(:)
      integer, intent(out), optional :: status
      integer :: i, freq, st
      real(dp) :: rate
      st = qbc_success
      freq = 1
      if (present(frequency)) freq = frequency
      if (size(rates) /= 1 .and. size(rates) < size(dates)) then
         allocate(factors(0)); st = qbc_size_mismatch
         if (present(status)) status = st
         return
      end if
      allocate(factors(size(dates)))
      do i = 1, size(dates)
         if (size(rates) == 1) then
            rate = rates(1)
         else
            rate = rates(i)
         end if
         factors(i) = discount_factor(rate, discount_time(analysis_date, dates(i)), rate_type, freq)
      end do
      if (present(status)) status = st
   end subroutine discount_factors

   subroutine curve_from_market_yields(terms, yields, nodes, approximation, rate_type, frequency, &
                                       output_curve, forward_output, status)
      real(dp), intent(in) :: terms(:), yields(:), nodes(:)
      integer, intent(in), optional :: approximation, rate_type, frequency
      type(qbc_curve), intent(out) :: output_curve
      logical, intent(in), optional :: forward_output
      integer, intent(out), optional :: status
      type(qbc_curve) :: base
      real(dp), allocatable :: values(:), fwd(:), out_terms(:)
      integer :: approx, rtype, freq, st
      logical :: fwd_out
      st = qbc_success
      approx = 1; rtype = qbc_rate_continuous; freq = 1; fwd_out = .false.
      if (present(approximation)) approx = approximation
      if (present(rate_type)) rtype = rate_type
      if (present(frequency)) freq = frequency
      if (present(forward_output)) fwd_out = forward_output
      if (size(terms) /= size(yields) .or. size(terms) == 0) then
         allocate(output_curve%terms(0), output_curve%rates(0)); st = qbc_size_mismatch
      else
         allocate(base%terms(size(terms)), base%rates(size(yields)))
         base%terms = terms
         if (rtype == qbc_rate_discrete) then
            base%rates = real(freq, dp) * log(1.0_dp + yields / real(freq, dp))
         else
            base%rates = yields
         end if
         base%approximation = approx
         base%rate_type = qbc_rate_continuous
         base%frequency = freq
         call interpolate_curve(base, nodes, values, st)
         if (fwd_out) then
            call spot_to_forward(nodes, values, approx, fwd, out_terms, st)
            allocate(output_curve%terms(size(out_terms)), output_curve%rates(size(fwd)))
            output_curve%terms = out_terms
            output_curve%rates = fwd
         else
            allocate(output_curve%terms(size(nodes)), output_curve%rates(size(nodes)))
            output_curve%terms = nodes
            output_curve%rates = values
         end if
         output_curve%approximation = approx
         output_curve%rate_type = qbc_rate_continuous
         output_curve%frequency = freq
      end if
      if (present(status)) status = st
   end subroutine curve_from_market_yields

end module qbc_curves
