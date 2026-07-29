! bondAnalyst modern Fortran port
! Copyright (C) 2022 MaheshP Kumar
! Copyright (C) 2026 Fortran port contributors
! SPDX-License-Identifier: GPL-3.0-only
module bondanalyst_money_market
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   use bondanalyst_kinds, only : dp
   use bondanalyst_support, only : ba_success, ba_invalid_argument, round_decimal, &
      quiet_nan, set_status
   implicit none
   private

   public :: pricingtbill, pricingmoneymarketinstrusingaor
   public :: fvmoneymarketinstrusingaor, computingaormoneymarketinstr
   public :: pricingcommercialpaper, computingquoteddiscratemmi
   public :: fvmmiusingquoteddiscrate

contains

   function pricingtbill(maturityval, daystomaturity, daysinyear, &
      mmquoteddiscrate, istat) result(value)
      real(dp), intent(in) :: maturityval, mmquoteddiscrate
      integer, intent(in) :: daystomaturity, daysinyear
      integer, optional, intent(out) :: istat
      real(dp) :: value, fraction

      if (.not. valid_day_count(daystomaturity, daysinyear) .or. &
         .not. ieee_is_finite(maturityval) .or. .not. ieee_is_finite(mmquoteddiscrate)) then
         value = quiet_nan()
         call set_status(istat, ba_invalid_argument)
      else
         fraction = real(daystomaturity, dp)/real(daysinyear, dp)
         value = abs(maturityval)*(1.0_dp-fraction*mmquoteddiscrate)
         value = round_decimal(value, 0)
         call set_status(istat, ba_success)
      end if
   end function pricingtbill


   function pricingmoneymarketinstrusingaor(maturityval, daystomaturity, &
      daysinyear, aor, istat) result(value)
      real(dp), intent(in) :: maturityval, aor
      integer, intent(in) :: daystomaturity, daysinyear
      integer, optional, intent(out) :: istat
      real(dp) :: value, denominator, fraction

      if (.not. valid_day_count(daystomaturity, daysinyear) .or. &
         .not. ieee_is_finite(maturityval) .or. .not. ieee_is_finite(aor)) then
         value = quiet_nan()
         call set_status(istat, ba_invalid_argument)
      else
         fraction = real(daystomaturity, dp)/real(daysinyear, dp)
         denominator = 1.0_dp + fraction*aor
         if (denominator <= 0.0_dp) then
            value = quiet_nan()
            call set_status(istat, ba_invalid_argument)
         else
            value = round_decimal(abs(maturityval)/denominator, 0)
            call set_status(istat, ba_success)
         end if
      end if
   end function pricingmoneymarketinstrusingaor


   function fvmoneymarketinstrusingaor(pvmmi, daystomaturity, daysinyear, &
      aor, istat) result(value)
      real(dp), intent(in) :: pvmmi, aor
      integer, intent(in) :: daystomaturity, daysinyear
      integer, optional, intent(out) :: istat
      real(dp) :: value, fraction

      if (.not. valid_day_count(daystomaturity, daysinyear) .or. &
         .not. ieee_is_finite(pvmmi) .or. .not. ieee_is_finite(aor)) then
         value = quiet_nan()
         call set_status(istat, ba_invalid_argument)
      else
         fraction = real(daystomaturity, dp)/real(daysinyear, dp)
         value = abs(pvmmi)*(1.0_dp+fraction*aor)
         value = round_decimal(value, 0)
         call set_status(istat, ba_success)
      end if
   end function fvmoneymarketinstrusingaor


   function computingaormoneymarketinstr(pvmmi, fvmmi, daystomaturity, &
      daysinyear, istat) result(value)
      real(dp), intent(in) :: pvmmi, fvmmi
      integer, intent(in) :: daystomaturity, daysinyear
      integer, optional, intent(out) :: istat
      real(dp) :: value

      if (.not. valid_day_count(daystomaturity, daysinyear) .or. abs(pvmmi) <= tiny(1.0_dp) .or. &
         .not. ieee_is_finite(pvmmi) .or. .not. ieee_is_finite(fvmmi)) then
         value = quiet_nan()
         call set_status(istat, ba_invalid_argument)
      else
         value = real(daysinyear, dp)/real(daystomaturity, dp)
         value = value*(fvmmi-pvmmi)/pvmmi
         value = round_decimal(value, 5)
         call set_status(istat, ba_success)
      end if
   end function computingaormoneymarketinstr


   function pricingcommercialpaper(maturityval, daystomaturity, daysinyear, &
      mmquoteddiscrate, istat) result(value)
      real(dp), intent(in) :: maturityval, mmquoteddiscrate
      integer, intent(in) :: daystomaturity, daysinyear
      integer, optional, intent(out) :: istat
      real(dp) :: value, fraction

      if (.not. valid_day_count(daystomaturity, daysinyear) .or. &
         .not. ieee_is_finite(maturityval) .or. .not. ieee_is_finite(mmquoteddiscrate)) then
         value = quiet_nan()
         call set_status(istat, ba_invalid_argument)
      else
         fraction = real(daystomaturity, dp)/real(daysinyear, dp)
         value = abs(maturityval)*(1.0_dp-fraction*mmquoteddiscrate)
         value = round_decimal(value, 2)
         call set_status(istat, ba_success)
      end if
   end function pricingcommercialpaper


   function computingquoteddiscratemmi(pvmmi, fvmmi, daystomaturity, &
      daysinyear, istat) result(value)
      real(dp), intent(in) :: pvmmi, fvmmi
      integer, intent(in) :: daystomaturity, daysinyear
      integer, optional, intent(out) :: istat
      real(dp) :: value

      if (.not. valid_day_count(daystomaturity, daysinyear) .or. abs(fvmmi) <= tiny(1.0_dp) .or. &
         .not. ieee_is_finite(pvmmi) .or. .not. ieee_is_finite(fvmmi)) then
         value = quiet_nan()
         call set_status(istat, ba_invalid_argument)
      else
         value = real(daysinyear, dp)/real(daystomaturity, dp)
         value = value*(fvmmi-pvmmi)/fvmmi
         value = round_decimal(value, 5)
         call set_status(istat, ba_success)
      end if
   end function computingquoteddiscratemmi


   function fvmmiusingquoteddiscrate(pvmmi, daystomaturity, daysinyear, &
      mmquoteddiscrate, istat) result(value)
      real(dp), intent(in) :: pvmmi, mmquoteddiscrate
      integer, intent(in) :: daystomaturity, daysinyear
      integer, optional, intent(out) :: istat
      real(dp) :: value, denominator, fraction

      if (.not. valid_day_count(daystomaturity, daysinyear) .or. &
         .not. ieee_is_finite(pvmmi) .or. .not. ieee_is_finite(mmquoteddiscrate)) then
         value = quiet_nan()
         call set_status(istat, ba_invalid_argument)
      else
         fraction = real(daystomaturity, dp)/real(daysinyear, dp)
         denominator = 1.0_dp-fraction*mmquoteddiscrate
         if (denominator <= 0.0_dp) then
            value = quiet_nan()
            call set_status(istat, ba_invalid_argument)
         else
            value = round_decimal(abs(pvmmi)/denominator, 0)
            call set_status(istat, ba_success)
         end if
      end if
   end function fvmmiusingquoteddiscrate


   pure function valid_day_count(daystomaturity, daysinyear) result(ok)
      integer, intent(in) :: daystomaturity, daysinyear
      logical :: ok

      ok = daystomaturity > 0 .and. daysinyear > 0
   end function valid_day_count

end module bondanalyst_money_market
