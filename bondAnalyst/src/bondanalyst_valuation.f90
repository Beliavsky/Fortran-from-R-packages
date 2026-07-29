! bondAnalyst modern Fortran port
! Copyright (C) 2022 MaheshP Kumar
! Copyright (C) 2026 Fortran port contributors
! SPDX-License-Identifier: GPL-3.0-only
module bondanalyst_valuation
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   use bondanalyst_kinds, only : dp
   use bondanalyst_support, only : ba_success, ba_invalid_argument, ba_size_mismatch, &
      round_decimal, quiet_nan, set_status, all_finite, valid_discount_base
   implicit none
   private

   public :: discouponpmtsbond, dismaturityvalbond, bondpriceyearlycoupons
   public :: pvcoupondeficiency, bondpricedefcoupon, pvexcesscoupon
   public :: bondpriceexcesscoupon, pricingzerocouponbond
   public :: pricingsacpnbond, pricingqtrlycpnbond, pricingwithspots
   public :: pricingwithsptseq, matrixmethod, pricingfrn, frpricing
   public :: pricingwithzspread, pricingwithgspread

contains

   function discouponpmtsbond(couponpmt, times, r, istat) result(value)
      real(dp), intent(in) :: couponpmt(:), times(:), r
      integer, optional, intent(out) :: istat
      real(dp) :: value

      if (size(couponpmt) == 0 .or. size(couponpmt) /= size(times)) then
         value = quiet_nan()
         call set_status(istat, ba_size_mismatch)
      else if (.not. all_finite(couponpmt) .or. .not. all_finite(times) .or. &
         any(times <= 0.0_dp) .or. r < 0.0_dp .or. .not. ieee_is_finite(r)) then
         value = quiet_nan()
         call set_status(istat, ba_invalid_argument)
      else
         value = round_decimal(sum(couponpmt/(1.0_dp+r)**times), 3)
         call set_status(istat, ba_success)
      end if
   end function discouponpmtsbond


   function dismaturityvalbond(bondmaturityval, n, r, istat) result(value)
      real(dp), intent(in) :: bondmaturityval, r
      integer, intent(in) :: n
      integer, optional, intent(out) :: istat
      real(dp) :: value

      if (n < 0 .or. .not. valid_discount_base(r) .or. &
         .not. ieee_is_finite(bondmaturityval)) then
         value = quiet_nan()
         call set_status(istat, ba_invalid_argument)
      else
         value = round_decimal(bondmaturityval/(1.0_dp+r)**n, 3)
         call set_status(istat, ba_success)
      end if
   end function dismaturityvalbond


   function bondpriceyearlycoupons(couponpmt, times, bondmaturityval, n, r, istat) &
      result(value)
      real(dp), intent(in) :: couponpmt(:), times(:), bondmaturityval, r
      integer, intent(in) :: n
      integer, optional, intent(out) :: istat
      real(dp) :: value

      if (size(couponpmt) == 0 .or. size(couponpmt) /= size(times)) then
         value = quiet_nan()
         call set_status(istat, ba_size_mismatch)
      else if (n < 0 .or. r < 0.0_dp .or. .not. all_finite(couponpmt) .or. &
         .not. all_finite(times) .or. any(times <= 0.0_dp) .or. &
         .not. ieee_is_finite(bondmaturityval) .or. .not. ieee_is_finite(r)) then
         value = quiet_nan()
         call set_status(istat, ba_invalid_argument)
      else
         value = sum(couponpmt/(1.0_dp+r)**times)
         value = value + bondmaturityval/(1.0_dp+r)**n
         value = round_decimal(value, 3)
         call set_status(istat, ba_success)
      end if
   end function bondpriceyearlycoupons


   function pvcoupondeficiency(coupondeficiency, times, r, istat) result(value)
      real(dp), intent(in) :: coupondeficiency(:), times(:), r
      integer, optional, intent(out) :: istat
      real(dp) :: value

      if (size(coupondeficiency) == 0 .or. size(coupondeficiency) /= size(times)) then
         value = quiet_nan()
         call set_status(istat, ba_size_mismatch)
      else if (.not. all_finite(coupondeficiency) .or. .not. all_finite(times) .or. &
         .not. valid_discount_base(r)) then
         value = quiet_nan()
         call set_status(istat, ba_invalid_argument)
      else
         value = round_decimal(sum(coupondeficiency/(1.0_dp+r)**times), 3)
         call set_status(istat, ba_success)
      end if
   end function pvcoupondeficiency


   function bondpricedefcoupon(parvalue, coupondeficiency, times, r, istat) result(value)
      real(dp), intent(in) :: parvalue, coupondeficiency(:), times(:), r
      integer, optional, intent(out) :: istat
      real(dp) :: value, deficiency

      if (size(coupondeficiency) == 0 .or. size(coupondeficiency) /= size(times)) then
         value = quiet_nan()
         call set_status(istat, ba_size_mismatch)
      else if (.not. all_finite(coupondeficiency) .or. .not. all_finite(times) .or. &
         .not. valid_discount_base(r) .or. .not. ieee_is_finite(parvalue)) then
         value = quiet_nan()
         call set_status(istat, ba_invalid_argument)
      else
         deficiency = sum(coupondeficiency/(1.0_dp+r)**times)
         value = round_decimal(parvalue-abs(deficiency), 3)
         call set_status(istat, ba_success)
      end if
   end function bondpricedefcoupon


   function pvexcesscoupon(couponexcess, times, r, istat) result(value)
      real(dp), intent(in) :: couponexcess(:), times(:), r
      integer, optional, intent(out) :: istat
      real(dp) :: value

      if (size(couponexcess) == 0 .or. size(couponexcess) /= size(times)) then
         value = quiet_nan()
         call set_status(istat, ba_size_mismatch)
      else if (.not. all_finite(couponexcess) .or. .not. all_finite(times) .or. &
         .not. valid_discount_base(r)) then
         value = quiet_nan()
         call set_status(istat, ba_invalid_argument)
      else
         value = round_decimal(sum(couponexcess/(1.0_dp+r)**times), 3)
         call set_status(istat, ba_success)
      end if
   end function pvexcesscoupon


   function bondpriceexcesscoupon(couponexcess, times, r, istat) result(value)
      real(dp), intent(in) :: couponexcess(:), times(:), r
      integer, optional, intent(out) :: istat
      real(dp) :: value

      if (size(couponexcess) == 0 .or. size(couponexcess) /= size(times)) then
         value = quiet_nan()
         call set_status(istat, ba_size_mismatch)
      else if (.not. all_finite(couponexcess) .or. .not. all_finite(times) .or. &
         .not. valid_discount_base(r)) then
         value = quiet_nan()
         call set_status(istat, ba_invalid_argument)
      else
         value = 100.0_dp + sum(couponexcess/(1.0_dp+r)**times)
         value = round_decimal(value, 3)
         call set_status(istat, ba_success)
      end if
   end function bondpriceexcesscoupon


   function pricingzerocouponbond(maturityval, n, r, istat) result(value)
      real(dp), intent(in) :: maturityval, r
      integer, intent(in) :: n
      integer, optional, intent(out) :: istat
      real(dp) :: value

      if (n < 0 .or. .not. valid_discount_base(r) .or. &
         .not. ieee_is_finite(maturityval)) then
         value = quiet_nan()
         call set_status(istat, ba_invalid_argument)
      else
         value = round_decimal(maturityval/(1.0_dp+r)**n, 3)
         call set_status(istat, ba_success)
      end if
   end function pricingzerocouponbond


   function pricingsacpnbond(sacoupons, times, maturityval, n, r, istat) result(value)
      real(dp), intent(in) :: sacoupons(:), times(:), maturityval, r
      integer, intent(in) :: n
      integer, optional, intent(out) :: istat
      real(dp) :: value, periodic_rate

      if (size(sacoupons) == 0 .or. size(sacoupons) /= size(times)) then
         value = quiet_nan()
         call set_status(istat, ba_size_mismatch)
      else if (n < 0 .or. r < 0.0_dp .or. .not. all_finite(sacoupons) .or. &
         .not. all_finite(times) .or. any(times <= 0.0_dp) .or. &
         .not. ieee_is_finite(maturityval) .or. .not. ieee_is_finite(r)) then
         value = quiet_nan()
         call set_status(istat, ba_invalid_argument)
      else
         periodic_rate = r/2.0_dp
         value = sum(sacoupons/(1.0_dp+periodic_rate)**times)
         value = value + maturityval/(1.0_dp+periodic_rate)**(2*n)
         value = round_decimal(value, 2)
         call set_status(istat, ba_success)
      end if
   end function pricingsacpnbond


   function pricingqtrlycpnbond(qcoupons, times, mv, n, r, istat) result(value)
      real(dp), intent(in) :: qcoupons(:), times(:), mv, r
      integer, intent(in) :: n
      integer, optional, intent(out) :: istat
      real(dp) :: value, periodic_rate

      if (size(qcoupons) == 0 .or. size(qcoupons) /= size(times)) then
         value = quiet_nan()
         call set_status(istat, ba_size_mismatch)
      else if (n < 0 .or. r < 0.0_dp .or. .not. all_finite(qcoupons) .or. &
         .not. all_finite(times) .or. any(times <= 0.0_dp) .or. &
         .not. ieee_is_finite(mv) .or. .not. ieee_is_finite(r)) then
         value = quiet_nan()
         call set_status(istat, ba_invalid_argument)
      else
         periodic_rate = r/4.0_dp
         value = sum(qcoupons/(1.0_dp+periodic_rate)**times)
         value = value + mv/(1.0_dp+periodic_rate)**(4*n)
         value = round_decimal(value, 2)
         call set_status(istat, ba_success)
      end if
   end function pricingqtrlycpnbond


   function pricingwithspots(coupons, spots, times, mv, n, istat) result(value)
      real(dp), intent(in) :: coupons(:), spots(:), times(:), mv
      integer, intent(in) :: n
      integer, optional, intent(out) :: istat
      real(dp) :: value

      if (size(coupons) == 0 .or. size(coupons) /= size(times) .or. &
         size(spots) /= size(coupons)) then
         value = quiet_nan()
         call set_status(istat, ba_size_mismatch)
      else if (n < 1 .or. n > size(spots) .or. .not. all_finite(coupons) .or. &
         .not. all_finite(spots) .or. .not. all_finite(times) .or. &
         any(times <= 0.0_dp) .or. any(1.0_dp+spots <= 0.0_dp) .or. &
         .not. ieee_is_finite(mv)) then
         value = quiet_nan()
         call set_status(istat, ba_invalid_argument)
      else
         value = sum(coupons/(1.0_dp+spots)**times)
         value = value + mv/(1.0_dp+spots(n))**n
         value = round_decimal(value, 2)
         call set_status(istat, ba_success)
      end if
   end function pricingwithspots


   function pricingwithsptseq(cpns, sp, t, mv, n, istat) result(value)
      real(dp), intent(in) :: cpns(:), sp(:), t(:), mv
      integer, intent(in) :: n
      integer, optional, intent(out) :: istat
      real(dp) :: value

      if (size(cpns) == 0 .or. size(cpns) /= size(t) .or. size(sp) /= size(cpns)) then
         value = quiet_nan()
         call set_status(istat, ba_size_mismatch)
      else if (n < 1 .or. n > size(sp) .or. .not. all_finite(cpns) .or. &
         .not. all_finite(sp) .or. .not. all_finite(t) .or. any(t <= 0.0_dp) .or. &
         any(1.0_dp+sp <= 0.0_dp) .or. .not. ieee_is_finite(mv)) then
         value = quiet_nan()
         call set_status(istat, ba_invalid_argument)
      else
         value = sum(cpns/(1.0_dp+sp)**t)
         value = value + mv/(1.0_dp+sp(n))**n
         value = round_decimal(value, 3)
         call set_status(istat, ba_success)
      end if
   end function pricingwithsptseq


   function matrixmethod(couponpmt, times, maturityval, n, r1, r2, istat) result(value)
      real(dp), intent(in) :: couponpmt(:), times(:), maturityval, r1, r2
      integer, intent(in) :: n
      integer, optional, intent(out) :: istat
      real(dp) :: value, r

      if (size(couponpmt) == 0 .or. size(couponpmt) /= size(times)) then
         value = quiet_nan()
         call set_status(istat, ba_size_mismatch)
      else
         r = 0.5_dp*(r1+r2)
         if (n < 0 .or. .not. all_finite(couponpmt) .or. .not. all_finite(times) .or. &
            any(times <= 0.0_dp) .or. .not. valid_discount_base(r) .or. &
            .not. ieee_is_finite(maturityval)) then
            value = quiet_nan()
            call set_status(istat, ba_invalid_argument)
         else
            value = sum(couponpmt/(1.0_dp+r)**times)
            value = value + maturityval/(1.0_dp+r)**n
            value = round_decimal(value, 3)
            call set_status(istat, ba_success)
         end if
      end if
   end function matrixmethod


   function pricingfrn(estrtrn, t, mv, maturityperiod, estdisc, istat) result(value)
      real(dp), intent(in) :: estrtrn(:), t(:), mv, estdisc
      integer, intent(in) :: maturityperiod
      integer, optional, intent(out) :: istat
      real(dp) :: value

      if (size(estrtrn) == 0 .or. size(estrtrn) /= size(t)) then
         value = quiet_nan()
         call set_status(istat, ba_size_mismatch)
      else if (maturityperiod < 0 .or. estdisc < 0.0_dp .or. &
         .not. all_finite(estrtrn) .or. .not. all_finite(t) .or. &
         any(t <= 0.0_dp) .or. .not. ieee_is_finite(mv) .or. &
         .not. ieee_is_finite(estdisc)) then
         value = quiet_nan()
         call set_status(istat, ba_invalid_argument)
      else
         value = sum(estrtrn/(1.0_dp+estdisc)**t)
         value = value + mv/(1.0_dp+estdisc)**maturityperiod
         value = round_decimal(value, 3)
         call set_status(istat, ba_success)
      end if
   end function pricingfrn


   function frpricing(cpns, fri, mv, n, istat) result(value)
      real(dp), intent(in) :: cpns(:), fri(:), mv
      integer, intent(in) :: n
      integer, optional, intent(out) :: istat
      real(dp) :: value

      if (size(cpns) == 0 .or. size(cpns) /= size(fri)) then
         value = quiet_nan()
         call set_status(istat, ba_size_mismatch)
      else if (n < 1 .or. n > size(fri) .or. .not. all_finite(cpns) .or. &
         .not. all_finite(fri) .or. any(abs(fri) <= tiny(1.0_dp)) .or. &
         .not. ieee_is_finite(mv)) then
         value = quiet_nan()
         call set_status(istat, ba_invalid_argument)
      else
         value = sum(cpns/fri) + mv/fri(n)
         value = round_decimal(value, 3)
         call set_status(istat, ba_success)
      end if
   end function frpricing


   function pricingwithzspread(cpns, spots, t, mv, n, zsprd, istat) result(value)
      real(dp), intent(in) :: cpns(:), spots(:), t(:), mv, zsprd
      integer, intent(in) :: n
      integer, optional, intent(out) :: istat
      real(dp) :: value

      if (size(cpns) == 0 .or. size(cpns) /= size(t) .or. size(spots) /= size(cpns)) then
         value = quiet_nan()
         call set_status(istat, ba_size_mismatch)
      else if (n < 1 .or. n > size(spots) .or. .not. all_finite(cpns) .or. &
         .not. all_finite(spots) .or. .not. all_finite(t) .or. any(t <= 0.0_dp) .or. &
         any(1.0_dp+zsprd+spots <= 0.0_dp) .or. .not. ieee_is_finite(mv) .or. &
         .not. ieee_is_finite(zsprd)) then
         value = quiet_nan()
         call set_status(istat, ba_invalid_argument)
      else
         value = sum(cpns/(1.0_dp+zsprd+spots)**t)
         value = value + mv/(1.0_dp+zsprd+spots(n))**n
         value = round_decimal(value, 3)
         call set_status(istat, ba_success)
      end if
   end function pricingwithzspread


   function pricingwithgspread(coupons, t, mv, n, ytmbenchgovtbond, gspread, istat) &
      result(value)
      real(dp), intent(in) :: coupons(:), t(:), mv, ytmbenchgovtbond, gspread
      integer, intent(in) :: n
      integer, optional, intent(out) :: istat
      real(dp) :: value, rate

      rate = ytmbenchgovtbond + gspread
      if (size(coupons) == 0 .or. size(coupons) /= size(t)) then
         value = quiet_nan()
         call set_status(istat, ba_size_mismatch)
      else if (n < 0 .or. ytmbenchgovtbond < 0.0_dp .or. .not. all_finite(coupons) .or. &
         .not. all_finite(t) .or. any(t <= 0.0_dp) .or. &
         .not. valid_discount_base(rate) .or. .not. ieee_is_finite(mv)) then
         value = quiet_nan()
         call set_status(istat, ba_invalid_argument)
      else
         value = sum(coupons/(1.0_dp+rate)**t) + mv/(1.0_dp+rate)**n
         value = round_decimal(value, 2)
         call set_status(istat, ba_success)
      end if
   end function pricingwithgspread

end module bondanalyst_valuation
