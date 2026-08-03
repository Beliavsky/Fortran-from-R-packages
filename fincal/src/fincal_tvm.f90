! SPDX-License-Identifier: GPL-2.0-or-later
module fincal_tvm
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_value, ieee_quiet_nan
   use fincal_kinds, only : dp
   use fincal_status, only : fincal_ok, fincal_invalid_input, fincal_no_root
   use fincal_types, only : root_result
   implicit none
   private

   public :: fv_simple, fv_annuity, fv, fv_uneven
   public :: pv_simple, pv_annuity, pv, pv_uneven, pv_perpetuity
   public :: pmt, n_period, discount_rate, npv, irr, irr2
   public :: perpetuity_rate, r_perpetuity, solve_discount_rate, solve_irr

   interface r_perpetuity
      module procedure perpetuity_rate
   end interface r_perpetuity

   integer, parameter :: root_discount_rate = 1
   integer, parameter :: root_irr = 2
contains
   elemental pure function fv_simple(r, n, present_value) result(value)
      real(dp), intent(in) :: r, n, present_value
      real(dp) :: value
      value = -present_value * (1.0_dp + r) ** n
   end function fv_simple

   elemental pure function fv_annuity(r, n, payment, payment_type) result(value)
      real(dp), intent(in) :: r, n, payment
      integer, intent(in), optional :: payment_type
      real(dp) :: value
      integer :: kind

      kind = 0
      if (present(payment_type)) kind = payment_type
      if (abs(r) <= sqrt(epsilon(1.0_dp))) then
         value = -payment * n
      else
         value = -payment / r * ((1.0_dp + r) ** n - 1.0_dp) * (1.0_dp + r) ** kind
      end if
   end function fv_annuity

   elemental pure function fv(r, n, present_value, payment, payment_type) result(value)
      real(dp), intent(in) :: r, n
      real(dp), intent(in), optional :: present_value, payment
      integer, intent(in), optional :: payment_type
      real(dp) :: value, pv_value, pmt_value
      integer :: kind

      pv_value = 0.0_dp
      pmt_value = 0.0_dp
      kind = 0
      if (present(present_value)) pv_value = present_value
      if (present(payment)) pmt_value = payment
      if (present(payment_type)) kind = payment_type
      value = fv_simple(r, n, pv_value) + fv_annuity(r, n, pmt_value, kind)
   end function fv

   pure function fv_uneven(r, cash_flows) result(value)
      real(dp), intent(in) :: r, cash_flows(:)
      real(dp) :: value
      integer :: i, n

      value = 0.0_dp
      n = size(cash_flows)
      do i = 1, n
         value = value + fv_simple(r, real(n - i, dp), cash_flows(i))
      end do
   end function fv_uneven

   elemental pure function pv_simple(r, n, future_value) result(value)
      real(dp), intent(in) :: r, n, future_value
      real(dp) :: value
      value = -future_value / (1.0_dp + r) ** n
   end function pv_simple

   elemental pure function pv_annuity(r, n, payment, payment_type) result(value)
      real(dp), intent(in) :: r, n, payment
      integer, intent(in), optional :: payment_type
      real(dp) :: value
      integer :: kind

      kind = 0
      if (present(payment_type)) kind = payment_type
      if (abs(r) <= sqrt(epsilon(1.0_dp))) then
         value = -payment * n
      else
         value = -payment / r * (1.0_dp - 1.0_dp / (1.0_dp + r) ** n) * (1.0_dp + r) ** kind
      end if
   end function pv_annuity

   elemental pure function pv(r, n, future_value, payment, payment_type) result(value)
      real(dp), intent(in) :: r, n
      real(dp), intent(in), optional :: future_value, payment
      integer, intent(in), optional :: payment_type
      real(dp) :: value, fv_value, pmt_value
      integer :: kind

      fv_value = 0.0_dp
      pmt_value = 0.0_dp
      kind = 0
      if (present(future_value)) fv_value = future_value
      if (present(payment)) pmt_value = payment
      if (present(payment_type)) kind = payment_type
      value = pv_simple(r, n, fv_value) + pv_annuity(r, n, pmt_value, kind)
   end function pv

   pure function pv_uneven(r, cash_flows) result(value)
      real(dp), intent(in) :: r, cash_flows(:)
      real(dp) :: value
      integer :: i

      value = 0.0_dp
      do i = 1, size(cash_flows)
         value = value + pv_simple(r, real(i, dp), cash_flows(i))
      end do
   end function pv_uneven

   pure function npv(r, cash_flows) result(value)
      real(dp), intent(in) :: r, cash_flows(:)
      real(dp) :: value, factor
      integer :: i

      if (size(cash_flows) == 0) then
         value = 0.0_dp
         return
      end if
      value = cash_flows(1)
      factor = 1.0_dp
      do i = 2, size(cash_flows)
         factor = factor / (1.0_dp + r)
         value = value + cash_flows(i) * factor
      end do
   end function npv

   elemental pure function pmt(r, n, present_value, future_value, payment_type) result(value)
      real(dp), intent(in) :: r, n, present_value, future_value
      integer, intent(in), optional :: payment_type
      real(dp) :: value
      integer :: kind

      kind = 0
      if (present(payment_type)) kind = payment_type
      if (abs(r) <= sqrt(epsilon(1.0_dp))) then
         value = -(present_value + future_value) / n
      else
         value = -(present_value + future_value / (1.0_dp + r) ** n) * r / &
            (1.0_dp - 1.0_dp / (1.0_dp + r) ** n) * (1.0_dp + r) ** (-kind)
      end if
   end function pmt

   elemental pure function n_period(r, present_value, future_value, payment, payment_type) result(value)
      real(dp), intent(in) :: r, present_value, future_value, payment
      integer, intent(in), optional :: payment_type
      real(dp) :: value
      integer :: kind

      kind = 0
      if (present(payment_type)) kind = payment_type
      if (abs(r) <= sqrt(epsilon(1.0_dp))) then
         value = -(present_value + future_value) / payment
      else
         value = log(-(future_value * r - payment * (1.0_dp + r) ** kind) / &
            (present_value * r + payment * (1.0_dp + r) ** kind)) / log(1.0_dp + r)
      end if
   end function n_period

   elemental pure function pv_perpetuity(r, payment, growth_rate, payment_type) result(value)
      real(dp), intent(in) :: r, payment
      real(dp), intent(in), optional :: growth_rate
      integer, intent(in), optional :: payment_type
      real(dp) :: value, growth
      integer :: kind

      growth = 0.0_dp
      kind = 0
      if (present(growth_rate)) growth = growth_rate
      if (present(payment_type)) kind = payment_type
      if (growth >= r) then
         value = ieee_value(0.0_dp, ieee_quiet_nan)
      else
         value = -payment / (r - growth) * (1.0_dp + r) ** kind
      end if
   end function pv_perpetuity

   elemental pure function perpetuity_rate(payment, present_value) result(value)
      real(dp), intent(in) :: payment, present_value
      real(dp) :: value
      value = -payment / present_value
   end function perpetuity_rate

   function discount_rate(n, present_value, future_value, payment, payment_type, status) result(rate)
      real(dp), intent(in) :: n, present_value, future_value, payment
      integer, intent(in), optional :: payment_type
      integer, intent(out), optional :: status
      real(dp) :: rate
      type(root_result) :: answer

      answer = solve_discount_rate(n, present_value, future_value, payment, payment_type)
      rate = answer%root
      if (present(status)) status = answer%status
   end function discount_rate

   function solve_discount_rate(n, present_value, future_value, payment, payment_type) result(answer)
      real(dp), intent(in) :: n, present_value, future_value, payment
      integer, intent(in), optional :: payment_type
      type(root_result) :: answer
      integer :: kind

      kind = 0
      if (present(payment_type)) kind = payment_type
      if (n <= 0.0_dp .or. (kind /= 0 .and. kind /= 1)) then
         answer%root = ieee_value(0.0_dp, ieee_quiet_nan)
         answer%function_value = answer%root
         answer%status = fincal_invalid_input
         return
      end if

      call find_rate_root(root_discount_rate, 1.0e-12_dp, 1.0e6_dp, answer, &
         n = n, present_value = present_value, future_value = future_value, &
         payment = payment, payment_type = kind)
   end function solve_discount_rate

   function irr(cash_flows, status, upper_rate) result(rate)
      real(dp), intent(in) :: cash_flows(:)
      integer, intent(out), optional :: status
      real(dp), intent(in), optional :: upper_rate
      real(dp) :: rate
      type(root_result) :: answer

      answer = solve_irr(cash_flows, 0.0_dp, upper_rate)
      rate = answer%root
      if (present(status)) status = answer%status
   end function irr

   function solve_irr(cash_flows, lower_rate, upper_rate) result(answer)
      real(dp), intent(in) :: cash_flows(:)
      real(dp), intent(in), optional :: lower_rate, upper_rate
      type(root_result) :: answer
      real(dp) :: lower, upper

      lower = 0.0_dp
      upper = 1.0e6_dp
      if (present(lower_rate)) lower = lower_rate
      if (present(upper_rate)) upper = upper_rate
      if (size(cash_flows) < 2 .or. lower <= -1.0_dp .or. upper <= lower) then
         answer%root = ieee_value(0.0_dp, ieee_quiet_nan)
         answer%function_value = answer%root
         answer%status = fincal_invalid_input
         return
      end if
      call find_rate_root(root_irr, lower, upper, answer, cash_flows = cash_flows)
   end function solve_irr

   function irr2(cash_flows, cutoff, from_rate, to_rate, step, status) result(rate)
      real(dp), intent(in) :: cash_flows(:)
      real(dp), intent(in), optional :: cutoff, from_rate, to_rate, step
      integer, intent(out), optional :: status
      real(dp) :: rate
      real(dp) :: threshold, first, last, increment, value
      integer :: i, n_steps

      threshold = 0.1_dp
      first = -1.0_dp
      last = 10.0_dp
      increment = 1.0e-6_dp
      if (present(cutoff)) threshold = cutoff
      if (present(from_rate)) first = from_rate
      if (present(to_rate)) last = to_rate
      if (present(step)) increment = step

      if (size(cash_flows) < 2 .or. threshold < 0.0_dp .or. increment <= 0.0_dp .or. last < first) then
         rate = ieee_value(0.0_dp, ieee_quiet_nan)
         if (present(status)) status = fincal_invalid_input
         return
      end if

      n_steps = int((last - first) / increment)
      do i = 0, n_steps
         rate = first + real(i, dp) * increment
         if (rate <= -1.0_dp) cycle
         value = npv(rate, cash_flows)
         if (ieee_is_finite(value)) then
            if (abs(value) < threshold) then
               if (present(status)) status = fincal_ok
               return
            end if
         end if
      end do
      rate = ieee_value(0.0_dp, ieee_quiet_nan)
      if (present(status)) status = fincal_no_root
   end function irr2

   subroutine find_rate_root(mode, lower_rate, upper_rate, answer, n, present_value, &
      future_value, payment, payment_type, cash_flows)
      integer, intent(in) :: mode
      real(dp), intent(in) :: lower_rate, upper_rate
      type(root_result), intent(out) :: answer
      real(dp), intent(in), optional :: n, present_value, future_value, payment
      integer, intent(in), optional :: payment_type
      real(dp), intent(in), optional :: cash_flows(:)
      integer, parameter :: scan_points = 4096, max_iterations = 200
      real(dp), parameter :: function_tolerance = 1.0e-12_dp
      real(dp) :: x0, x1, xa, xb, xm, fa, fb, fm, rate_a, rate_b, rate_m
      integer :: i, iteration

      answer%root = ieee_value(0.0_dp, ieee_quiet_nan)
      answer%function_value = answer%root
      answer%iterations = 0
      answer%status = fincal_no_root

      if (lower_rate <= -1.0_dp .or. upper_rate <= lower_rate) then
         answer%status = fincal_invalid_input
         return
      end if
      if (.not. objective_arguments_valid(mode, n, present_value, future_value, payment, &
         payment_type, cash_flows)) then
         answer%status = fincal_invalid_input
         return
      end if

      x0 = log(1.0_dp + lower_rate)
      x1 = log(1.0_dp + upper_rate)
      xa = x0
      rate_a = exp(xa) - 1.0_dp
      fa = evaluate_objective(mode, rate_a, n, present_value, future_value, payment, &
         payment_type, cash_flows)
      if (ieee_is_finite(fa) .and. abs(fa) <= function_tolerance) then
         answer%root = rate_a
         answer%function_value = fa
         answer%status = fincal_ok
         return
      end if

      do i = 1, scan_points
         xb = x0 + (x1 - x0) * real(i, dp) / real(scan_points, dp)
         rate_b = exp(xb) - 1.0_dp
         fb = evaluate_objective(mode, rate_b, n, present_value, future_value, payment, &
            payment_type, cash_flows)
         if (.not. ieee_is_finite(fb)) then
            xa = xb
            fa = fb
            cycle
         end if
         if (abs(fb) <= function_tolerance) then
            answer%root = rate_b
            answer%function_value = fb
            answer%status = fincal_ok
            answer%iterations = i
            return
         end if
         if (ieee_is_finite(fa) .and. opposite_signs(fa, fb)) then
            do iteration = 1, max_iterations
               answer%iterations = iteration
               xm = 0.5_dp * (xa + xb)
               rate_m = exp(xm) - 1.0_dp
               fm = evaluate_objective(mode, rate_m, n, present_value, future_value, &
                  payment, payment_type, cash_flows)
               if (.not. ieee_is_finite(fm)) then
                  xb = xm
                  cycle
               end if
               if (abs(fm) <= function_tolerance .or. &
                  abs(xb - xa) <= 8.0_dp * epsilon(1.0_dp) * max(1.0_dp, abs(xm))) then
                  answer%root = rate_m
                  answer%function_value = fm
                  answer%status = fincal_ok
                  return
               end if
               if (opposite_signs(fa, fm)) then
                  xb = xm
               else
                  xa = xm
                  fa = fm
               end if
            end do
            answer%root = exp(0.5_dp * (xa + xb)) - 1.0_dp
            answer%function_value = evaluate_objective(mode, answer%root, n, present_value, &
               future_value, payment, payment_type, cash_flows)
            answer%status = fincal_ok
            return
         end if
         xa = xb
         fa = fb
      end do
   end subroutine find_rate_root

   pure function objective_arguments_valid(mode, n, present_value, future_value, payment, &
      payment_type, cash_flows) result(valid)
      integer, intent(in) :: mode
      real(dp), intent(in), optional :: n, present_value, future_value, payment
      integer, intent(in), optional :: payment_type
      real(dp), intent(in), optional :: cash_flows(:)
      logical :: valid

      select case (mode)
      case (root_discount_rate)
         valid = present(n) .and. present(present_value) .and. present(future_value) .and. &
            present(payment) .and. present(payment_type)
      case (root_irr)
         valid = present(cash_flows)
      case default
         valid = .false.
      end select
   end function objective_arguments_valid

   pure function evaluate_objective(mode, rate, n, present_value, future_value, payment, &
      payment_type, cash_flows) result(value)
      integer, intent(in) :: mode
      real(dp), intent(in) :: rate
      real(dp), intent(in), optional :: n, present_value, future_value, payment
      integer, intent(in), optional :: payment_type
      real(dp), intent(in), optional :: cash_flows(:)
      real(dp) :: value

      select case (mode)
      case (root_discount_rate)
         value = discount_rate_objective(rate, n, present_value, future_value, payment, payment_type)
      case (root_irr)
         value = npv(rate, cash_flows)
      case default
         value = ieee_value(0.0_dp, ieee_quiet_nan)
      end select
   end function evaluate_objective


   pure function discount_rate_objective(rate, n, present_value, future_value, payment, &
      payment_type) result(value)
      real(dp), intent(in) :: rate, n, present_value, future_value, payment
      integer, intent(in) :: payment_type
      real(dp) :: value, discount, annuity_factor, exponent

      if (abs(rate) <= sqrt(epsilon(1.0_dp))) then
         discount = 1.0_dp
         annuity_factor = n
      else
         exponent = n * log(1.0_dp + rate)
         if (exponent >= -log(tiny(1.0_dp))) then
            discount = 0.0_dp
         else
            discount = exp(-exponent)
         end if
         annuity_factor = (1.0_dp - discount) / rate * (1.0_dp + rate) ** payment_type
      end if
      value = -present_value - future_value * discount - payment * annuity_factor
   end function discount_rate_objective

   elemental pure function opposite_signs(a, b) result(opposite)
      real(dp), intent(in) :: a, b
      logical :: opposite
      opposite = (a < 0.0_dp .and. b > 0.0_dp) .or. (a > 0.0_dp .and. b < 0.0_dp)
   end function opposite_signs
end module fincal_tvm
