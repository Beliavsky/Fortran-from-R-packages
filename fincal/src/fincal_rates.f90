! SPDX-License-Identifier: GPL-2.0-or-later
module fincal_rates
   use fincal_kinds, only : dp
   use fincal_status, only : fincal_ok, fincal_invalid_input
   implicit none
   private

   public :: eir, ear, ear_continuous, ear_to_bey, ear_to_hpr
   public :: bdy, bdy_to_mmy, hpr, hpr_to_bey, hpr_to_ear
   public :: hpr_to_mmy, mmy_to_hpr, continuous_rate, nominal_rate
   public :: r_continuous, r_norminal
   public :: ear2bey, ear2hpr, bdy2mmy, hpr2bey, hpr2ear, hpr2mmy, mmy2hpr

   interface ear2bey
      module procedure ear_to_bey
   end interface ear2bey
   interface ear2hpr
      module procedure ear_to_hpr
   end interface ear2hpr
   interface bdy2mmy
      module procedure bdy_to_mmy
   end interface bdy2mmy
   interface hpr2bey
      module procedure hpr_to_bey
   end interface hpr2bey
   interface hpr2ear
      module procedure hpr_to_ear
   end interface hpr2ear
   interface hpr2mmy
      module procedure hpr_to_mmy
   end interface hpr2mmy
   interface mmy2hpr
      module procedure mmy_to_hpr
   end interface mmy2hpr
contains
   function eir(r, n, p, rate_type, status) result(value)
      real(dp), intent(in) :: r
      integer, intent(in), optional :: n, p
      character(len=*), intent(in), optional :: rate_type
      integer, intent(out), optional :: status
      real(dp) :: value
      integer :: nn, pp
      character(len=1) :: kind

      nn = 1
      if (present(n)) nn = n
      pp = 12
      if (present(p)) pp = p
      kind = 'e'
      if (present(rate_type)) kind = lower_first(rate_type)

      if (nn <= 0 .or. pp <= 0 .or. (kind /= 'e' .and. kind /= 'p')) then
         value = 0.0_dp
         if (present(status)) status = fincal_invalid_input
         return
      end if

      if (kind == 'e') then
         value = (1.0_dp + r / real(nn, dp)) ** (real(nn, dp) / real(pp, dp)) - 1.0_dp
      else
         value = r / real(pp, dp)
      end if
      if (present(status)) status = fincal_ok
   end function eir

   elemental pure function ear(r, m) result(value)
      real(dp), intent(in) :: r
      integer, intent(in) :: m
      real(dp) :: value
      value = (1.0_dp + r / real(m, dp)) ** m - 1.0_dp
   end function ear

   elemental pure function ear_continuous(r) result(value)
      real(dp), intent(in) :: r
      real(dp) :: value
      value = exp(r) - 1.0_dp
   end function ear_continuous

   elemental pure function ear_to_bey(effective_annual_rate) result(value)
      real(dp), intent(in) :: effective_annual_rate
      real(dp) :: value
      value = 2.0_dp * (sqrt(1.0_dp + effective_annual_rate) - 1.0_dp)
   end function ear_to_bey

   elemental pure function ear_to_hpr(effective_annual_rate, days) result(value)
      real(dp), intent(in) :: effective_annual_rate, days
      real(dp) :: value
      value = (1.0_dp + effective_annual_rate) ** (days / 365.0_dp) - 1.0_dp
   end function ear_to_hpr

   elemental pure function bdy(discount, face_value, days) result(value)
      real(dp), intent(in) :: discount, face_value, days
      real(dp) :: value
      value = 360.0_dp * discount / (face_value * days)
   end function bdy

   elemental pure function bdy_to_mmy(bank_discount_yield, days) result(value)
      real(dp), intent(in) :: bank_discount_yield, days
      real(dp) :: value
      value = 360.0_dp * bank_discount_yield / (360.0_dp - days * bank_discount_yield)
   end function bdy_to_mmy

   elemental pure function hpr(ending_value, beginning_value, cash_flow_received) result(value)
      real(dp), intent(in) :: ending_value, beginning_value
      real(dp), intent(in), optional :: cash_flow_received
      real(dp) :: value, cash_flow
      cash_flow = 0.0_dp
      if (present(cash_flow_received)) cash_flow = cash_flow_received
      value = (ending_value - beginning_value + cash_flow) / beginning_value
   end function hpr

   elemental pure function hpr_to_bey(holding_period_return, months) result(value)
      real(dp), intent(in) :: holding_period_return, months
      real(dp) :: value
      value = 2.0_dp * ((1.0_dp + holding_period_return) ** (6.0_dp / months) - 1.0_dp)
   end function hpr_to_bey

   elemental pure function hpr_to_ear(holding_period_return, days) result(value)
      real(dp), intent(in) :: holding_period_return, days
      real(dp) :: value
      value = (1.0_dp + holding_period_return) ** (365.0_dp / days) - 1.0_dp
   end function hpr_to_ear

   elemental pure function hpr_to_mmy(holding_period_return, days) result(value)
      real(dp), intent(in) :: holding_period_return, days
      real(dp) :: value
      value = 360.0_dp * holding_period_return / days
   end function hpr_to_mmy

   elemental pure function mmy_to_hpr(money_market_yield, days) result(value)
      real(dp), intent(in) :: money_market_yield, days
      real(dp) :: value
      value = money_market_yield * days / 360.0_dp
   end function mmy_to_hpr

   elemental pure function continuous_rate(r, m) result(value)
      real(dp), intent(in) :: r
      integer, intent(in) :: m
      real(dp) :: value
      value = real(m, dp) * log(1.0_dp + r / real(m, dp))
   end function continuous_rate

   elemental pure function nominal_rate(rc, m) result(value)
      real(dp), intent(in) :: rc
      integer, intent(in) :: m
      real(dp) :: value
      value = real(m, dp) * (exp(rc / real(m, dp)) - 1.0_dp)
   end function nominal_rate

   elemental pure function r_continuous(r, m) result(value)
      real(dp), intent(in) :: r
      integer, intent(in) :: m
      real(dp) :: value
      value = continuous_rate(r, m)
   end function r_continuous

   elemental pure function r_norminal(rc, m) result(value)
      real(dp), intent(in) :: rc
      integer, intent(in) :: m
      real(dp) :: value
      value = nominal_rate(rc, m)
   end function r_norminal

   pure function lower_first(text) result(ch)
      character(len=*), intent(in) :: text
      character(len=1) :: ch
      integer :: code
      if (len_trim(text) == 0) then
         ch = ' '
      else
         ch = text(1:1)
         code = iachar(ch)
         if (code >= iachar('A') .and. code <= iachar('Z')) ch = achar(code + 32)
      end if
   end function lower_first
end module fincal_rates
