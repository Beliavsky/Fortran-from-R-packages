! Copyright (C) 1998-2025 T. W. Yee, University of Auckland.
! Modern Fortran computational translation, 2026.
! SPDX-License-Identifier: GPL-3.0-only
module vgam_inflated
   use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
   use vgam_kinds, only : dp
   use vgam_distributions, only : dbeta_v, pbeta_v, qbeta_v, dbetabinom_ab, pbetabinom, rbetabinom
   implicit none
   private
   public :: dzoabeta, pzoabeta, qzoabeta, rzoabeta
   public :: dzoibetabinom_ab, pzoibetabinom_ab, rzoibetabinom_ab
contains

   elemental real(dp) function dzoabeta(x, shape1, shape2, pobs0, pobs1) result(d)
      real(dp), intent(in) :: x, shape1, shape2, pobs0, pobs1
      if (.not. valid_probabilities(pobs0, pobs1) .or. shape1 <= 0.0_dp .or. shape2 <= 0.0_dp) then
         d = ieee_value(0.0_dp, ieee_quiet_nan)
      else if (x == 0.0_dp .and. pobs0 > 0.0_dp) then
         d = pobs0
      else if (x == 1.0_dp .and. pobs1 > 0.0_dp) then
         d = pobs1
      else
         d = (1.0_dp - pobs0 - pobs1)*dbeta_v(x, shape1, shape2)
      end if
   end function dzoabeta

   elemental real(dp) function pzoabeta(q, shape1, shape2, pobs0, pobs1) result(p)
      real(dp), intent(in) :: q, shape1, shape2, pobs0, pobs1
      if (.not. valid_probabilities(pobs0, pobs1) .or. shape1 <= 0.0_dp .or. shape2 <= 0.0_dp) then
         p = ieee_value(0.0_dp, ieee_quiet_nan)
      else if (q < 0.0_dp) then
         p = 0.0_dp
      else if (q >= 1.0_dp) then
         p = 1.0_dp
      else
         p = pobs0 + (1.0_dp - pobs0 - pobs1)*pbeta_v(q, shape1, shape2)
         p = min(1.0_dp, max(0.0_dp, p))
      end if
   end function pzoabeta

   real(dp) function qzoabeta(probability, shape1, shape2, pobs0, pobs1) result(q)
      real(dp), intent(in) :: probability, shape1, shape2, pobs0, pobs1
      real(dp) :: base
      if (.not. valid_probabilities(pobs0, pobs1) .or. shape1 <= 0.0_dp .or. shape2 <= 0.0_dp .or. &
          probability < 0.0_dp .or. probability > 1.0_dp) then
         q = ieee_value(0.0_dp, ieee_quiet_nan)
      else if (probability <= pobs0) then
         q = 0.0_dp
      else if (probability >= 1.0_dp - pobs1) then
         q = 1.0_dp
      else
         base = (probability - pobs0)/(1.0_dp - pobs0 - pobs1)
         q = qbeta_v(base, shape1, shape2)
      end if
   end function qzoabeta

   real(dp) function rzoabeta(shape1, shape2, pobs0, pobs1) result(x)
      real(dp), intent(in) :: shape1, shape2, pobs0, pobs1
      real(dp) :: u
      call random_number(u)
      x = qzoabeta(u, shape1, shape2, pobs0, pobs1)
   end function rzoabeta

   elemental real(dp) function dzoibetabinom_ab(x, size, shape1, shape2, pstr0, pstrsize) result(d)
      integer, intent(in) :: x, size
      real(dp), intent(in) :: shape1, shape2, pstr0, pstrsize
      if (size < 0 .or. shape1 <= 0.0_dp .or. shape2 <= 0.0_dp .or. &
          .not. valid_probabilities(pstr0, pstrsize)) then
         d = ieee_value(0.0_dp, ieee_quiet_nan)
      else
         d = (1.0_dp - pstr0 - pstrsize)*dbetabinom_ab(x, size, shape1, shape2)
         if (x == 0) d = d + pstr0
         if (x == size) d = d + pstrsize
      end if
   end function dzoibetabinom_ab

   real(dp) function pzoibetabinom_ab(q, size, shape1, shape2, pstr0, pstrsize) result(p)
      integer, intent(in) :: q, size
      real(dp), intent(in) :: shape1, shape2, pstr0, pstrsize
      if (size < 0 .or. shape1 <= 0.0_dp .or. shape2 <= 0.0_dp .or. &
          .not. valid_probabilities(pstr0, pstrsize)) then
         p = ieee_value(0.0_dp, ieee_quiet_nan)
      else if (q < 0) then
         p = 0.0_dp
      else if (q >= size) then
         p = 1.0_dp
      else
         p = (1.0_dp - pstr0 - pstrsize)*pbetabinom(q, size, shape1, shape2) + pstr0
         p = min(1.0_dp, max(0.0_dp, p))
      end if
   end function pzoibetabinom_ab

   integer function rzoibetabinom_ab(size, shape1, shape2, pstr0, pstrsize) result(x)
      integer, intent(in) :: size
      real(dp), intent(in) :: shape1, shape2, pstr0, pstrsize
      real(dp) :: u
      if (size < 0 .or. shape1 <= 0.0_dp .or. shape2 <= 0.0_dp .or. &
          .not. valid_probabilities(pstr0, pstrsize)) then
         x = -1
         return
      end if
      call random_number(u)
      if (u < pstr0) then
         x = 0
      else if (u > 1.0_dp - pstrsize) then
         x = size
      else
         x = rbetabinom(size, shape1, shape2)
      end if
   end function rzoibetabinom_ab

   elemental logical function valid_probabilities(p0, p1) result(ok)
      real(dp), intent(in) :: p0, p1
      ok = p0 >= 0.0_dp .and. p1 >= 0.0_dp .and. p0 + p1 <= 1.0_dp
   end function valid_probabilities

end module vgam_inflated
