! SPDX-License-Identifier: GPL-3.0-or-later
module qbc_swaps
   use qbc_kinds, only : dp
   use qbc_status, only : qbc_success, qbc_invalid_argument
   use qbc_types, only : qbc_swap, qbc_bond, qbc_curve, qbc_asset_fixed, &
      qbc_leg_fixed_fixed, qbc_leg_fixed_float, qbc_leg_float_fixed, qbc_leg_float_float
   use qbc_bonds, only : valuation_bonds
   implicit none
   private
   public :: valuation_swaps, swap_leg_value

contains

   real(dp) function swap_leg_value(maturity, analysis_date, frequency, rate_type, daycount, &
                                    coupon_rate, spread, principal, is_fixed, curve, status) result(value)
      use qbc_dates, only : qbc_date
      type(qbc_date), intent(in) :: maturity, analysis_date
      integer, intent(in) :: frequency, rate_type
      character(len=*), intent(in) :: daycount
      real(dp), intent(in) :: coupon_rate, spread, principal
      logical, intent(in) :: is_fixed
      type(qbc_curve), intent(in) :: curve
      integer, intent(out), optional :: status
      type(qbc_bond) :: bond
      real(dp) :: c(1)
      integer :: st
      bond%maturity = maturity
      bond%analysis_date = analysis_date
      bond%coupon_rate = coupon_rate
      bond%principal = principal
      bond%spread = spread
      bond%asset_type = qbc_asset_fixed
      bond%frequency = frequency
      bond%rate_type = rate_type
      bond%daycount = daycount
      if (is_fixed) then
         c(1) = coupon_rate
         value = valuation_bonds(bond, c, curve, dirty=.true., status=st)
      else
         c(1) = 0.0_dp
         value = principal + valuation_bonds(bond, c, curve, dirty=.true., spread_only=.true., status=st)
      end if
      if (present(status)) status = st
   end function swap_leg_value

   real(dp) function valuation_swaps(swap, curve1, curve2, basis, status) result(value)
      type(qbc_swap), intent(in) :: swap
      type(qbc_curve), intent(in) :: curve1, curve2
      type(qbc_curve), intent(in), optional :: basis
      integer, intent(out), optional :: status
      logical :: fixed1, fixed2
      real(dp) :: leg1, leg2
      integer :: st
      st = qbc_success
      select case (swap%legs)
      case (qbc_leg_fixed_fixed)
         fixed1 = .true.; fixed2 = .true.
      case (qbc_leg_fixed_float)
         fixed1 = .true.; fixed2 = .false.
      case (qbc_leg_float_fixed)
         fixed1 = .false.; fixed2 = .true.
      case (qbc_leg_float_float)
         fixed1 = .false.; fixed2 = .false.
      case default
         fixed1 = .true.; fixed2 = .true.; st = qbc_invalid_argument
      end select
      leg1 = swap_leg_value(swap%maturity, swap%analysis_date, swap%frequency, swap%rate_type, &
                            swap%daycount, swap%coupon_rate1, swap%spread1, swap%principal1, fixed1, curve1)
      if (present(basis)) then
         leg2 = swap_leg_value(swap%maturity, swap%analysis_date, swap%frequency, swap%rate_type, &
                               swap%daycount, swap%coupon_rate2, swap%spread2, swap%principal2, fixed2, basis)
      else
         leg2 = swap_leg_value(swap%maturity, swap%analysis_date, swap%frequency, swap%rate_type, &
                               swap%daycount, swap%coupon_rate2, swap%spread2, swap%principal2, fixed2, curve2)
      end if
      value = leg1 - swap%exchange_rate * leg2
      if (present(status)) status = st
   end function valuation_swaps

end module qbc_swaps
