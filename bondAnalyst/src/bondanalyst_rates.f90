! bondAnalyst modern Fortran port
! Copyright (C) 2022 MaheshP Kumar
! Copyright (C) 2026 Fortran port contributors
! SPDX-License-Identifier: GPL-3.0-only
module bondanalyst_rates
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   use bondanalyst_kinds, only : dp
   use bondanalyst_support, only : ba_success, ba_invalid_argument, ba_size_mismatch, &
      ba_no_root, ba_out_of_range, round_decimal, quiet_nan, set_status, &
      all_finite
   implicit none
   private

   public :: ytmzerocouponbond
   public :: computingbondytmratefivedecimalplaces
   public :: computingbondytmratesixdecimalplaces
   public :: convertaprtodifferentperiodcity, extracompensationforhigherrisk
   public :: annualytmzcbforperiodicity, earzcbvariousperiodicity
   public :: returnincomefrn, periodicdiscratefrn, discmarginfrn
   public :: computingytc, computingparrate, forwards, saforwards
   public :: computinggspread, computingzspread
   public :: effectiveannualratezcb

contains

   function ytmzerocouponbond(maturityval, n, price, istat) result(value)
      real(dp), intent(in) :: maturityval, price
      integer, intent(in) :: n
      integer, optional, intent(out) :: istat
      real(dp) :: value

      if (maturityval <= 0.0_dp .or. price <= 0.0_dp .or. n <= 0 .or. &
         .not. ieee_is_finite(maturityval) .or. .not. ieee_is_finite(price)) then
         value = quiet_nan()
         call set_status(istat, ba_invalid_argument)
      else
         value = (maturityval/price)**(1.0_dp/real(n, dp)) - 1.0_dp
         value = round_decimal(value, 5)
         call set_status(istat, ba_success)
      end if
   end function ytmzerocouponbond


   function computingbondytmratefivedecimalplaces(couponpmt, mv, bondpv, period, istat) &
      result(value)
      real(dp), intent(in) :: couponpmt, mv, bondpv
      integer, intent(in) :: period
      integer, optional, intent(out) :: istat
      real(dp) :: value, root
      integer :: status

      if (.not. valid_level_cashflow_inputs(couponpmt, mv, bondpv, period)) then
         value = quiet_nan()
         call set_status(istat, ba_invalid_argument)
         return
      end if
      status = solve_level_cashflow_root(couponpmt, mv, bondpv, period, root)
      if (status == ba_success) then
         value = round_decimal(root, 5)
      else
         value = quiet_nan()
      end if
      call set_status(istat, status)
   end function computingbondytmratefivedecimalplaces


   function computingbondytmratesixdecimalplaces(couponpmt, mv, bondpv, period, istat) &
      result(value)
      real(dp), intent(in) :: couponpmt, mv, bondpv
      integer, intent(in) :: period
      integer, optional, intent(out) :: istat
      real(dp) :: value, root
      integer :: status

      if (.not. valid_level_cashflow_inputs(couponpmt, mv, bondpv, period)) then
         value = quiet_nan()
         call set_status(istat, ba_invalid_argument)
         return
      end if
      status = solve_level_cashflow_root(couponpmt, mv, bondpv, period, root)
      if (status == ba_success) then
         value = round_decimal(root, 6)
      else
         value = quiet_nan()
      end if
      call set_status(istat, status)
   end function computingbondytmratesixdecimalplaces


   function convertaprtodifferentperiodcity(givenapr, givenperiodicity, &
      desiredperiodicity, istat) result(value)
      real(dp), intent(in) :: givenapr
      integer, intent(in) :: givenperiodicity, desiredperiodicity
      integer, optional, intent(out) :: istat
      real(dp) :: value

      if (givenperiodicity <= 0 .or. desiredperiodicity <= 0 .or. &
         1.0_dp+givenapr/real(givenperiodicity, dp) <= 0.0_dp .or. &
         .not. ieee_is_finite(givenapr)) then
         value = quiet_nan()
         call set_status(istat, ba_invalid_argument)
      else
         value = (1.0_dp+givenapr/real(givenperiodicity, dp))** &
            (real(givenperiodicity, dp)/real(desiredperiodicity, dp))
         value = (value-1.0_dp)*real(desiredperiodicity, dp)
         value = round_decimal(value, 4)
         call set_status(istat, ba_success)
      end if
   end function convertaprtodifferentperiodcity


   function extracompensationforhigherrisk(aprofriskybond, &
      aprofcomparablebond, istat) result(value)
      real(dp), intent(in) :: aprofriskybond, aprofcomparablebond
      integer, optional, intent(out) :: istat
      real(dp) :: value

      if (.not. ieee_is_finite(aprofriskybond) .or. &
         .not. ieee_is_finite(aprofcomparablebond)) then
         value = quiet_nan()
         call set_status(istat, ba_invalid_argument)
      else
         value = (aprofriskybond-aprofcomparablebond)*10000.0_dp
         call set_status(istat, ba_success)
      end if
   end function extracompensationforhigherrisk


   function annualytmzcbforperiodicity(maturityval, yearstomaturity, zcbprice, &
      desiredperiodicity, istat) result(value)
      real(dp), intent(in) :: maturityval, yearstomaturity, zcbprice
      integer, intent(in) :: desiredperiodicity
      integer, optional, intent(out) :: istat
      real(dp) :: value, periods

      periods = yearstomaturity*real(desiredperiodicity, dp)
      if (maturityval <= 0.0_dp .or. zcbprice <= 0.0_dp .or. periods <= 0.0_dp .or. &
         .not. ieee_is_finite(maturityval) .or. .not. ieee_is_finite(yearstomaturity) .or. &
         .not. ieee_is_finite(zcbprice)) then
         value = quiet_nan()
         call set_status(istat, ba_invalid_argument)
      else
         value = ((maturityval/zcbprice)**(1.0_dp/periods)-1.0_dp)* &
            real(desiredperiodicity, dp)
         value = round_decimal(value, 6)
         call set_status(istat, ba_success)
      end if
   end function annualytmzcbforperiodicity


   function earzcbvariousperiodicity(maturityval, yearstomaturity, zcbprice, &
      desiredperiodicity, istat) result(value)
      real(dp), intent(in) :: maturityval, yearstomaturity, zcbprice
      integer, intent(in) :: desiredperiodicity
      integer, optional, intent(out) :: istat
      real(dp) :: value, periods

      periods = yearstomaturity*real(desiredperiodicity, dp)
      if (maturityval <= 0.0_dp .or. zcbprice <= 0.0_dp .or. periods <= 0.0_dp .or. &
         .not. ieee_is_finite(maturityval) .or. .not. ieee_is_finite(yearstomaturity) .or. &
         .not. ieee_is_finite(zcbprice)) then
         value = quiet_nan()
         call set_status(istat, ba_invalid_argument)
      else
         ! Source-compatible behavior: this is the periodic rate, despite the R name.
         value = (maturityval/zcbprice)**(1.0_dp/periods)-1.0_dp
         value = round_decimal(value, 6)
         call set_status(istat, ba_success)
      end if
   end function earzcbvariousperiodicity


   function effectiveannualratezcb(maturityval, yearstomaturity, zcbprice, istat) &
      result(value)
      real(dp), intent(in) :: maturityval, yearstomaturity, zcbprice
      integer, optional, intent(out) :: istat
      real(dp) :: value

      if (maturityval <= 0.0_dp .or. zcbprice <= 0.0_dp .or. yearstomaturity <= 0.0_dp .or. &
         .not. ieee_is_finite(maturityval) .or. .not. ieee_is_finite(yearstomaturity) .or. &
         .not. ieee_is_finite(zcbprice)) then
         value = quiet_nan()
         call set_status(istat, ba_invalid_argument)
      else
         value = (maturityval/zcbprice)**(1.0_dp/yearstomaturity)-1.0_dp
         call set_status(istat, ba_success)
      end if
   end function effectiveannualratezcb


   function returnincomefrn(index, qtdmargin, maturityval, periodicity, istat) result(value)
      real(dp), intent(in) :: index, qtdmargin, maturityval
      integer, intent(in) :: periodicity
      integer, optional, intent(out) :: istat
      real(dp) :: value

      if (periodicity <= 0 .or. .not. ieee_is_finite(index) .or. &
         .not. ieee_is_finite(qtdmargin) .or. .not. ieee_is_finite(maturityval)) then
         value = quiet_nan()
         call set_status(istat, ba_invalid_argument)
      else
         value = (index+qtdmargin)*maturityval/real(periodicity, dp)
         value = round_decimal(value, 6)
         call set_status(istat, ba_success)
      end if
   end function returnincomefrn


   function periodicdiscratefrn(estrtrn, mvfrn, pricefrn, maturityyears, &
      periodicity, istat) result(value)
      real(dp), intent(in) :: estrtrn, mvfrn, pricefrn
      integer, intent(in) :: maturityyears, periodicity
      integer, optional, intent(out) :: istat
      real(dp) :: value, root
      integer :: periods, status

      periods = maturityyears*periodicity
      if (.not. valid_level_cashflow_inputs(estrtrn, mvfrn, pricefrn, periods)) then
         value = quiet_nan()
         call set_status(istat, ba_invalid_argument)
         return
      end if
      status = solve_level_cashflow_root(estrtrn, mvfrn, pricefrn, periods, root)
      if (status == ba_success) then
         value = round_decimal(root, 6)
      else
         value = quiet_nan()
      end if
      call set_status(istat, status)
   end function periodicdiscratefrn


   function discmarginfrn(index, estrtrn, mvfrn, pricefrn, maturityyears, &
      periodicity, istat) result(value)
      real(dp), intent(in) :: index, estrtrn, mvfrn, pricefrn
      integer, intent(in) :: maturityyears, periodicity
      integer, optional, intent(out) :: istat
      real(dp) :: value, root
      integer :: periods, status

      periods = maturityyears*periodicity
      if (.not. ieee_is_finite(index) .or. &
         .not. valid_level_cashflow_inputs(estrtrn, mvfrn, pricefrn, periods)) then
         value = quiet_nan()
         call set_status(istat, ba_invalid_argument)
         return
      end if
      status = solve_level_cashflow_root(estrtrn, mvfrn, pricefrn, periods, root)
      if (status == ba_success) then
         root = round_decimal(root, 6)
         value = round_decimal(root*real(periodicity, dp)-index, 6)
      else
         value = quiet_nan()
      end if
      call set_status(istat, status)
   end function discmarginfrn


   function computingytc(couponpmt, callableval, bondpv, maturityyears, &
      ytcyears, istat) result(value)
      real(dp), intent(in) :: couponpmt, callableval, bondpv
      integer, intent(in) :: maturityyears, ytcyears
      integer, optional, intent(out) :: istat
      real(dp) :: value, root
      integer :: status

      if (ytcyears > maturityyears) then
         value = quiet_nan()
         call set_status(istat, ba_out_of_range)
         return
      end if
      if (.not. valid_level_cashflow_inputs(couponpmt, callableval, bondpv, ytcyears)) then
         value = quiet_nan()
         call set_status(istat, ba_invalid_argument)
         return
      end if
      status = solve_level_cashflow_root(couponpmt, callableval, bondpv, ytcyears, root)
      if (status == ba_success) then
         value = round_decimal(root, 5)
      else
         value = quiet_nan()
      end if
      call set_status(istat, status)
   end function computingytc


   function computingparrate(spotrates, times, mv, pv, n, istat) result(value)
      real(dp), intent(in) :: spotrates(:), times(:), mv, pv
      integer, intent(in) :: n
      integer, optional, intent(out) :: istat
      real(dp) :: value, numerator, denominator

      if (size(spotrates) == 0 .or. size(spotrates) /= size(times)) then
         value = quiet_nan()
         call set_status(istat, ba_size_mismatch)
      else if (n < 1 .or. n > size(spotrates) .or. .not. all_finite(spotrates) .or. &
         .not. all_finite(times) .or. any(times <= 0.0_dp) .or. &
         any(1.0_dp+spotrates <= 0.0_dp) .or. .not. ieee_is_finite(mv) .or. &
         .not. ieee_is_finite(pv)) then
         value = quiet_nan()
         call set_status(istat, ba_invalid_argument)
      else
         numerator = 1.0_dp - 1.0_dp/(1.0_dp+spotrates(n))**n
         denominator = sum(1.0_dp/(1.0_dp+spotrates)**times)
         if (abs(denominator) <= tiny(1.0_dp)) then
            value = quiet_nan()
            call set_status(istat, ba_invalid_argument)
         else
            value = round_decimal(100.0_dp*numerator/denominator, 3)
            call set_status(istat, ba_success)
         end if
      end if
   end function computingparrate


   function forwards(spots, yrsfrbegins, yrsfrapplies, t, n, istat) result(value)
      real(dp), intent(in) :: spots(:), t(:)
      integer, intent(in) :: yrsfrbegins, yrsfrapplies, n
      integer, optional, intent(out) :: istat
      real(dp) :: value, numerator, denominator
      integer :: ending

      ending = yrsfrbegins + yrsfrapplies
      if (size(spots) == 0 .or. size(spots) /= size(t)) then
         value = quiet_nan()
         call set_status(istat, ba_size_mismatch)
      else if (yrsfrbegins <= 0 .or. yrsfrapplies <= 0 .or. n /= size(spots) .or. &
         ending > size(spots) .or. .not. all_finite(spots) .or. .not. all_finite(t) .or. &
         any(t <= 0.0_dp) .or. any(1.0_dp+spots <= 0.0_dp)) then
         value = quiet_nan()
         call set_status(istat, ba_invalid_argument)
      else
         numerator = (1.0_dp+spots(ending))**ending
         denominator = (1.0_dp+spots(yrsfrbegins))**yrsfrbegins
         value = (numerator/denominator)**(1.0_dp/real(yrsfrapplies, dp))-1.0_dp
         value = round_decimal(value, 4)
         call set_status(istat, ba_success)
      end if
   end function forwards


   function saforwards(spots, bgn, aply, times, n, istat) result(value)
      real(dp), intent(in) :: spots(:), times(:)
      integer, intent(in) :: bgn, aply, n
      integer, optional, intent(out) :: istat
      real(dp) :: value, numerator, denominator
      integer :: ending

      ending = bgn + aply
      if (size(spots) == 0 .or. size(spots) /= size(times)) then
         value = quiet_nan()
         call set_status(istat, ba_size_mismatch)
      else if (bgn <= 0 .or. aply <= 0 .or. n /= size(spots) .or. ending > size(spots) .or. &
         .not. all_finite(spots) .or. .not. all_finite(times) .or. any(times <= 0.0_dp) .or. &
         any(1.0_dp+spots/2.0_dp <= 0.0_dp)) then
         value = quiet_nan()
         call set_status(istat, ba_invalid_argument)
      else
         numerator = (1.0_dp+spots(ending)/2.0_dp)**ending
         denominator = (1.0_dp+spots(bgn)/2.0_dp)**bgn
         value = 2.0_dp*((numerator/denominator)**(1.0_dp/real(aply, dp))-1.0_dp)
         value = round_decimal(value, 6)
         call set_status(istat, ba_success)
      end if
   end function saforwards


   function computinggspread(ytmcorpbond, ytmbenchgovtbond, istat) result(value)
      real(dp), intent(in) :: ytmcorpbond, ytmbenchgovtbond
      integer, optional, intent(out) :: istat
      real(dp) :: value

      if (.not. ieee_is_finite(ytmcorpbond) .or. .not. ieee_is_finite(ytmbenchgovtbond)) then
         value = quiet_nan()
         call set_status(istat, ba_invalid_argument)
      else
         value = ytmcorpbond-ytmbenchgovtbond
         call set_status(istat, ba_success)
      end if
   end function computinggspread


   function computingzspread(coupons, mv, bondpv, n, spots, istat) result(value)
      real(dp), intent(in) :: coupons, mv, bondpv, spots(:)
      integer, intent(in) :: n
      integer, optional, intent(out) :: istat
      real(dp) :: value, root
      integer :: status

      if (n < 1 .or. n > size(spots) .or. coupons < 0.0_dp .or. mv <= 0.0_dp .or. &
         bondpv <= 0.0_dp .or. .not. all_finite(spots) .or. any(1.0_dp+spots <= 0.0_dp) .or. &
         .not. ieee_is_finite(coupons) .or. .not. ieee_is_finite(mv) .or. &
         .not. ieee_is_finite(bondpv)) then
         value = quiet_nan()
         call set_status(istat, ba_invalid_argument)
         return
      end if
      status = solve_zspread_root(coupons, mv, bondpv, n, spots, root)
      if (status == ba_success) then
         value = round_decimal(root, 4)
      else
         value = quiet_nan()
      end if
      call set_status(istat, status)
   end function computingzspread


   function solve_level_cashflow_root(coupon, principal, target_price, periods, root) &
      result(status)
      real(dp), intent(in) :: coupon, principal, target_price
      integer, intent(in) :: periods
      real(dp), intent(out) :: root
      integer :: status
      integer :: iter
      real(dp) :: a, b, c, fa, fb, fc, tolerance

      tolerance = 1.0e-12_dp
      a = 0.0_dp
      fa = level_cashflow_price(coupon, principal, periods, a)-target_price
      if (abs(fa) <= tolerance) then
         root = 0.0_dp
         status = ba_success
         return
      end if
      if (fa < 0.0_dp) then
         root = quiet_nan()
         status = ba_no_root
         return
      end if
      b = 0.01_dp
      fb = level_cashflow_price(coupon, principal, periods, b)-target_price
      do while (fb > 0.0_dp .and. b < 1.0e8_dp)
         b = 2.0_dp*b+0.01_dp
         fb = level_cashflow_price(coupon, principal, periods, b)-target_price
      end do
      if (fb > 0.0_dp) then
         root = quiet_nan()
         status = ba_no_root
         return
      end if
      do iter = 1, 256
         c = 0.5_dp*(a+b)
         fc = level_cashflow_price(coupon, principal, periods, c)-target_price
         if (abs(fc) <= tolerance .or. 0.5_dp*(b-a) <= &
            tolerance*max(1.0_dp, abs(c))) then
            root = c
            status = ba_success
            return
         end if
         if (fc > 0.0_dp) then
            a = c
            fa = fc
         else
            b = c
            fb = fc
         end if
      end do
      root = 0.5_dp*(a+b)
      status = ba_no_root
   end function solve_level_cashflow_root


   function solve_zspread_root(coupon, principal, target_price, periods, &
      spots, root) result(status)
      real(dp), intent(in) :: coupon, principal, target_price, spots(:)
      integer, intent(in) :: periods
      real(dp), intent(out) :: root
      integer :: status
      integer :: iter
      real(dp) :: a, b, c, fa, fb, fc, tolerance

      tolerance = 1.0e-12_dp
      a = 0.0_dp
      fa = zspread_price(coupon, principal, periods, spots, a)-target_price
      if (abs(fa) <= tolerance) then
         root = 0.0_dp
         status = ba_success
         return
      end if
      if (fa < 0.0_dp) then
         root = quiet_nan()
         status = ba_no_root
         return
      end if
      b = 0.01_dp
      fb = zspread_price(coupon, principal, periods, spots, b)-target_price
      do while (fb > 0.0_dp .and. b < 1.0e8_dp)
         b = 2.0_dp*b+0.01_dp
         fb = zspread_price(coupon, principal, periods, spots, b)-target_price
      end do
      if (fb > 0.0_dp) then
         root = quiet_nan()
         status = ba_no_root
         return
      end if
      do iter = 1, 256
         c = 0.5_dp*(a+b)
         fc = zspread_price(coupon, principal, periods, spots, c)-target_price
         if (abs(fc) <= tolerance .or. 0.5_dp*(b-a) <= &
            tolerance*max(1.0_dp, abs(c))) then
            root = c
            status = ba_success
            return
         end if
         if (fc > 0.0_dp) then
            a = c
            fa = fc
         else
            b = c
            fb = fc
         end if
      end do
      root = 0.5_dp*(a+b)
      status = ba_no_root
   end function solve_zspread_root


   pure function zspread_price(coupon, principal, periods, spots, spread) result(price)
      real(dp), intent(in) :: coupon, principal, spots(:), spread
      integer, intent(in) :: periods
      real(dp) :: price
      integer :: i

      price = 0.0_dp
      do i = 1, periods-1
         price = price + coupon/(1.0_dp+spots(i)+spread)**i
      end do
      price = price + (coupon+principal)/(1.0_dp+spots(periods)+spread)**periods
   end function zspread_price


   pure function valid_level_cashflow_inputs(coupon, principal, price, periods) result(ok)
      real(dp), intent(in) :: coupon, principal, price
      integer, intent(in) :: periods
      logical :: ok

      ok = periods > 0 .and. coupon >= 0.0_dp .and. principal > 0.0_dp .and. &
         price > 0.0_dp .and. ieee_is_finite(coupon) .and. &
         ieee_is_finite(principal) .and. ieee_is_finite(price)
   end function valid_level_cashflow_inputs


   pure function level_cashflow_price(coupon, principal, periods, rate) result(price)
      real(dp), intent(in) :: coupon, principal, rate
      integer, intent(in) :: periods
      real(dp) :: price
      integer :: i

      price = 0.0_dp
      do i = 1, periods-1
         price = price + coupon/(1.0_dp+rate)**i
      end do
      price = price + (coupon+principal)/(1.0_dp+rate)**periods
   end function level_cashflow_price

end module bondanalyst_rates
