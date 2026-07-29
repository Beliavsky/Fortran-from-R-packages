! bondAnalyst modern Fortran port
! Copyright (C) 2022 MaheshP Kumar
! Copyright (C) 2026 Fortran port contributors
! SPDX-License-Identifier: GPL-3.0-only
module bondanalyst_duration
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   use bondanalyst_kinds, only : dp
   use bondanalyst_support, only : ba_success, ba_invalid_argument, round_decimal, &
      quiet_nan, set_status
   implicit none
   private

   public :: aiactdtcon, airoundeddaysconv
   public :: macduration, pvfullprice, macdurationonfp
   public :: macdurationoncouponrate, modifduration
   public :: modifdurationusingmacduration
   public :: estimatedpercentchangepvfullprice, approxmodifduration
   public :: approxmacdurationusingapprmodifduration, effdurtncallablebond
   public :: moneyduration, changepvfullbondprice, computingbondpvbp
   public :: conventionalpercentpricechange

contains

   function aiactdtcon(cpmt, dt1, dt2, stdt, istat) result(value)
      real(dp), intent(in) :: cpmt
      integer, intent(in) :: dt1, dt2, stdt
      integer, optional, intent(out) :: istat
      real(dp) :: value
      integer :: days_from_last_coupon, days_between_coupons

      days_from_last_coupon = stdt-dt1
      days_between_coupons = dt2-dt1
      if (days_between_coupons == 0 .or. .not. ieee_is_finite(cpmt)) then
         value = quiet_nan()
         call set_status(istat, ba_invalid_argument)
      else
         value = 0.5_dp*cpmt*real(days_from_last_coupon, dp)/ &
            real(days_between_coupons, dp)
         value = round_decimal(value, 6)
         call set_status(istat, ba_success)
      end if
   end function aiactdtcon


   function airoundeddaysconv(cpmt, bfrstldt, stldt, elpsmnths, &
      daysbtwncpns, istat) result(value)
      real(dp), intent(in) :: cpmt
      integer, intent(in) :: bfrstldt, stldt, elpsmnths, daysbtwncpns
      integer, optional, intent(out) :: istat
      real(dp) :: value
      integer :: days_from_last_coupon

      days_from_last_coupon = stldt-bfrstldt
      if (daysbtwncpns == 0 .or. .not. ieee_is_finite(cpmt)) then
         value = quiet_nan()
         call set_status(istat, ba_invalid_argument)
      else
         value = 0.5_dp*cpmt*real(days_from_last_coupon+30*elpsmnths, dp)/ &
            real(daysbtwncpns, dp)
         value = round_decimal(value, 6)
         call set_status(istat, ba_success)
      end if
   end function airoundeddaysconv


   function macduration(n, ytm, coupon, maturityval, dayscpntosettle, &
      dayscouponperiod, istat) result(value)
      integer, intent(in) :: n, dayscpntosettle, dayscouponperiod
      real(dp), intent(in) :: ytm, coupon, maturityval
      integer, optional, intent(out) :: istat
      real(dp) :: value, numerator, denominator, tt
      integer :: status

      call duration_components(n, ytm, coupon, maturityval, dayscpntosettle, &
         dayscouponperiod, numerator, denominator, tt, status)
      if (status /= ba_success .or. abs(denominator) <= tiny(1.0_dp)) then
         value = quiet_nan()
         if (status == ba_success) status = ba_invalid_argument
      else
         value = round_decimal(numerator/denominator, 4)
      end if
      call set_status(istat, status)
   end function macduration


   function pvfullprice(n, ytm, coupon, maturityval, dayscpntosettle, &
      dayscouponperiod, istat) result(value)
      integer, intent(in) :: n, dayscpntosettle, dayscouponperiod
      real(dp), intent(in) :: ytm, coupon, maturityval
      integer, optional, intent(out) :: istat
      real(dp) :: value, numerator, denominator, tt
      integer :: status

      call duration_components(n, ytm, coupon, maturityval, dayscpntosettle, &
         dayscouponperiod, numerator, denominator, tt, status)
      if (status /= ba_success) then
         value = quiet_nan()
      else
         value = round_decimal(denominator, 4)
      end if
      call set_status(istat, status)
   end function pvfullprice


   function macdurationonfp(fp, n, ytm, cpn, mv, dayscpntosettle, &
      dayscouponperiod, istat) result(value)
      real(dp), intent(in) :: fp, ytm, cpn, mv
      integer, intent(in) :: n, dayscpntosettle, dayscouponperiod
      integer, optional, intent(out) :: istat
      real(dp) :: value, numerator, denominator, tt
      integer :: status

      call duration_components(n, ytm, cpn, mv, dayscpntosettle, &
         dayscouponperiod, numerator, denominator, tt, status)
      if (status /= ba_success .or. abs(fp) <= tiny(1.0_dp) .or. .not. ieee_is_finite(fp)) then
         value = quiet_nan()
         if (status == ba_success) status = ba_invalid_argument
      else
         value = round_decimal(numerator/fp, 4)
      end if
      call set_status(istat, status)
   end function macdurationonfp


   function macdurationoncouponrate(couponrate, n, ytm, tt, istat) result(value)
      real(dp), intent(in) :: couponrate, ytm, tt
      integer, intent(in) :: n
      integer, optional, intent(out) :: istat
      real(dp) :: value, leftside, right_num, right_den

      if (n <= 0 .or. abs(ytm) <= tiny(1.0_dp) .or. 1.0_dp+ytm <= 0.0_dp .or. &
         .not. ieee_is_finite(couponrate) .or. .not. ieee_is_finite(ytm) .or. &
         .not. ieee_is_finite(tt)) then
         value = quiet_nan()
         call set_status(istat, ba_invalid_argument)
         return
      end if
      leftside = (1.0_dp+ytm)/ytm
      right_num = 1.0_dp+ytm+real(n, dp)*(couponrate-ytm)
      right_den = couponrate*((1.0_dp+ytm)**n-1.0_dp)+ytm
      if (abs(right_den) <= tiny(1.0_dp)) then
         value = quiet_nan()
         call set_status(istat, ba_invalid_argument)
      else
         value = round_decimal(leftside-right_num/right_den-tt, 4)
         call set_status(istat, ba_success)
      end if
   end function macdurationoncouponrate


   function modifduration(n, ytm, coupon, maturityval, dayscpntosettle, &
      dayscouponperiod, istat) result(value)
      integer, intent(in) :: n, dayscpntosettle, dayscouponperiod
      real(dp), intent(in) :: ytm, coupon, maturityval
      integer, optional, intent(out) :: istat
      real(dp) :: value, numerator, denominator, tt
      integer :: status

      call duration_components(n, ytm, coupon, maturityval, dayscpntosettle, &
         dayscouponperiod, numerator, denominator, tt, status)
      if (status /= ba_success .or. abs(denominator) <= tiny(1.0_dp) .or. 1.0_dp+abs(ytm) <= tiny(1.0_dp)) then
         value = quiet_nan()
         if (status == ba_success) status = ba_invalid_argument
      else
         value = numerator/denominator/(1.0_dp+ytm)
         value = round_decimal(value, 4)
      end if
      call set_status(istat, status)
   end function modifduration


   function modifdurationusingmacduration(macduration_value, ytm, istat) result(value)
      real(dp), intent(in) :: macduration_value, ytm
      integer, optional, intent(out) :: istat
      real(dp) :: value

      if (1.0_dp+abs(ytm) <= tiny(1.0_dp) .or. .not. ieee_is_finite(macduration_value) .or. &
         .not. ieee_is_finite(ytm)) then
         value = quiet_nan()
         call set_status(istat, ba_invalid_argument)
      else
         value = round_decimal(macduration_value/(1.0_dp+ytm), 4)
         call set_status(istat, ba_success)
      end if
   end function modifdurationusingmacduration


   function estimatedpercentchangepvfullprice(annualmodifduration, &
      changeinannualytm, istat) result(value)
      real(dp), intent(in) :: annualmodifduration, changeinannualytm
      integer, optional, intent(out) :: istat
      real(dp) :: value

      if (.not. ieee_is_finite(annualmodifduration) .or. &
         .not. ieee_is_finite(changeinannualytm)) then
         value = quiet_nan()
         call set_status(istat, ba_invalid_argument)
      else
         ! Source-compatible sign convention from the R package.
         value = round_decimal(annualmodifduration*changeinannualytm, 6)
         call set_status(istat, ba_success)
      end if
   end function estimatedpercentchangepvfullprice


   function conventionalpercentpricechange(modified_duration, yield_change, istat) &
      result(value)
      real(dp), intent(in) :: modified_duration, yield_change
      integer, optional, intent(out) :: istat
      real(dp) :: value

      if (.not. ieee_is_finite(modified_duration) .or. &
         .not. ieee_is_finite(yield_change)) then
         value = quiet_nan()
         call set_status(istat, ba_invalid_argument)
      else
         value = -modified_duration*yield_change
         call set_status(istat, ba_success)
      end if
   end function conventionalpercentpricechange


   function approxmodifduration(pvbase, pvplus, pvminus, percentchangeytm, &
      istat) result(value)
      real(dp), intent(in) :: pvbase, pvplus, pvminus, percentchangeytm
      integer, optional, intent(out) :: istat
      real(dp) :: value, denominator

      denominator = 2.0_dp*percentchangeytm*pvbase
      if (abs(denominator) <= tiny(1.0_dp) .or. .not. ieee_is_finite(pvbase) .or. &
         .not. ieee_is_finite(pvplus) .or. .not. ieee_is_finite(pvminus) .or. &
         .not. ieee_is_finite(percentchangeytm)) then
         value = quiet_nan()
         call set_status(istat, ba_invalid_argument)
      else
         value = (pvminus-pvplus)/denominator
         call set_status(istat, ba_success)
      end if
   end function approxmodifduration


   function approxmacdurationusingapprmodifduration(approxmodifduration_value, &
      periodicytm, istat) result(value)
      real(dp), intent(in) :: approxmodifduration_value, periodicytm
      integer, optional, intent(out) :: istat
      real(dp) :: value

      if (.not. ieee_is_finite(approxmodifduration_value) .or. &
         .not. ieee_is_finite(periodicytm)) then
         value = quiet_nan()
         call set_status(istat, ba_invalid_argument)
      else
         value = approxmodifduration_value*(1.0_dp+periodicytm)
         value = round_decimal(value, 3)
         call set_status(istat, ba_success)
      end if
   end function approxmacdurationusingapprmodifduration


   function effdurtncallablebond(pvbase, pvplus, pvminus, perchangebenchytm, &
      istat) result(value)
      real(dp), intent(in) :: pvbase, pvplus, pvminus, perchangebenchytm
      integer, optional, intent(out) :: istat
      real(dp) :: value, denominator

      denominator = 2.0_dp*perchangebenchytm*pvbase
      if (abs(denominator) <= tiny(1.0_dp) .or. .not. ieee_is_finite(pvbase) .or. &
         .not. ieee_is_finite(pvplus) .or. .not. ieee_is_finite(pvminus) .or. &
         .not. ieee_is_finite(perchangebenchytm)) then
         value = quiet_nan()
         call set_status(istat, ba_invalid_argument)
      else
         value = round_decimal((pvminus-pvplus)/denominator, 4)
         call set_status(istat, ba_success)
      end if
   end function effdurtncallablebond


   function moneyduration(macduration_value, ytm, pvfullbondprice, istat) result(value)
      real(dp), intent(in) :: macduration_value, ytm, pvfullbondprice
      integer, optional, intent(out) :: istat
      real(dp) :: value

      if (1.0_dp+abs(ytm) <= tiny(1.0_dp) .or. .not. ieee_is_finite(macduration_value) .or. &
         .not. ieee_is_finite(ytm) .or. .not. ieee_is_finite(pvfullbondprice)) then
         value = quiet_nan()
         call set_status(istat, ba_invalid_argument)
      else
         value = macduration_value/(1.0_dp+ytm)*pvfullbondprice
         value = round_decimal(value, 4)
         call set_status(istat, ba_success)
      end if
   end function moneyduration


   function changepvfullbondprice(moneyduration_value, changeytm, istat) result(value)
      real(dp), intent(in) :: moneyduration_value, changeytm
      integer, optional, intent(out) :: istat
      real(dp) :: value

      if (.not. ieee_is_finite(moneyduration_value) .or. .not. ieee_is_finite(changeytm)) then
         value = quiet_nan()
         call set_status(istat, ba_invalid_argument)
      else
         value = round_decimal(-moneyduration_value*changeytm, 4)
         call set_status(istat, ba_success)
      end if
   end function changepvfullbondprice


   function computingbondpvbp(pvplus, pvminus, istat) result(value)
      real(dp), intent(in) :: pvplus, pvminus
      integer, optional, intent(out) :: istat
      real(dp) :: value

      if (.not. ieee_is_finite(pvplus) .or. .not. ieee_is_finite(pvminus)) then
         value = quiet_nan()
         call set_status(istat, ba_invalid_argument)
      else
         value = round_decimal(0.5_dp*(pvminus-pvplus), 4)
         call set_status(istat, ba_success)
      end if
   end function computingbondpvbp


   subroutine duration_components(n, ytm, coupon, maturityval, dayscpntosettle, &
      dayscouponperiod, numerator, denominator, tt, status)
      integer, intent(in) :: n, dayscpntosettle, dayscouponperiod
      real(dp), intent(in) :: ytm, coupon, maturityval
      real(dp), intent(out) :: numerator, denominator, tt
      integer, intent(out) :: status
      integer :: i
      real(dp) :: time_i

      numerator = quiet_nan()
      denominator = quiet_nan()
      tt = 0.0_dp
      if (n <= 0 .or. ytm < 0.0_dp .or. .not. ieee_is_finite(ytm) .or. &
         .not. ieee_is_finite(coupon) .or. .not. ieee_is_finite(maturityval)) then
         status = ba_invalid_argument
         return
      end if
      if (dayscpntosettle /= 0 .and. dayscouponperiod /= 0) then
         tt = real(dayscpntosettle, dp)/real(dayscouponperiod, dp)
      end if

      time_i = real(n, dp)-tt
      numerator = time_i*(coupon+maturityval)/(1.0_dp+ytm)**time_i
      denominator = (coupon+maturityval)/(1.0_dp+ytm)**time_i
      do i = 1, n-1
         time_i = real(i, dp)-tt
         numerator = numerator + time_i*coupon/(1.0_dp+ytm)**time_i
         denominator = denominator + coupon/(1.0_dp+ytm)**time_i
      end do
      status = ba_success
   end subroutine duration_components

end module bondanalyst_duration
