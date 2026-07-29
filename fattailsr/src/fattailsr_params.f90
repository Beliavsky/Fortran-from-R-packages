! FatTailsR modern Fortran translation
! Copyright (C) 2014-2026 Patrice Kiener
! Licensed under GPL-2.0-only. See COPYING.
module fattailsr_params
   use fattailsr_kinds, only : dp
   implicit none
   private

   type, public :: kiener_parameters
      real(dp) :: m = 0.0_dp
      real(dp) :: g = 1.0_dp
      real(dp) :: a = 3.2_dp
      real(dp) :: k = 3.2_dp
      real(dp) :: w = 3.2_dp
      real(dp) :: d = 0.0_dp
      real(dp) :: e = 0.0_dp
   end type kiener_parameters

   public :: aw2k, aw2d, aw2e, ad2e, ad2k, ad2w, ae2d, ae2k, ae2w
   public :: ak2d, ak2e, ak2w, de2a, de2k, de2w, dk2a, dk2e, dk2w
   public :: dw2a, dw2e, dw2k, ek2a, ek2d, ek2w, ew2a, ew2d, ew2k
   public :: kd2a, kd2e, kd2w, ke2a, ke2d, ke2w, kw2a, kw2d, kw2e
   public :: make_k1, make_k2, make_k3, make_k4, parameters_valid

contains

   elemental pure function aw2k(a, w) result(k)
      real(dp), intent(in) :: a, w
      real(dp) :: k
      k = 2.0_dp/(1.0_dp/a + 1.0_dp/w)
   end function aw2k

   elemental pure function aw2d(a, w) result(d)
      real(dp), intent(in) :: a, w
      real(dp) :: d
      d = (-1.0_dp/a + 1.0_dp/w)/2.0_dp
   end function aw2d

   elemental pure function aw2e(a, w) result(e)
      real(dp), intent(in) :: a, w
      real(dp) :: e
      e = (a - w)/(a + w)
   end function aw2e


   elemental pure function ad2e(a, d) result(e)
      real(dp), intent(in) :: a, d
      real(dp) :: e
      e = d*a/(1.0_dp + d*a)
   end function ad2e

   elemental pure function ad2k(a, d) result(k)
      real(dp), intent(in) :: a, d
      real(dp) :: k
      k = a/(1.0_dp + d*a)
   end function ad2k

   elemental pure function ad2w(a, d) result(w)
      real(dp), intent(in) :: a, d
      real(dp) :: w
      w = 1.0_dp/(2.0_dp*d + 1.0_dp/a)
   end function ad2w

   elemental pure function ae2d(a, e) result(d)
      real(dp), intent(in) :: a, e
      real(dp) :: d
      d = e/((1.0_dp-e)*a)
   end function ae2d

   elemental pure function ae2k(a, e) result(k)
      real(dp), intent(in) :: a, e
      real(dp) :: k
      k = (1.0_dp-e)*a
   end function ae2k

   elemental pure function ae2w(a, e) result(w)
      real(dp), intent(in) :: a, e
      real(dp) :: w
      w = a*(1.0_dp-e)/(1.0_dp+e)
   end function ae2w

   elemental pure function ak2d(a, k) result(d)
      real(dp), intent(in) :: a, k
      real(dp) :: d
      d = (a-k)/(a*k)
   end function ak2d

   elemental pure function ak2e(a, k) result(e)
      real(dp), intent(in) :: a, k
      real(dp) :: e
      e = (a-k)/a
   end function ak2e

   elemental pure function ak2w(a, k) result(w)
      real(dp), intent(in) :: a, k
      real(dp) :: w
      w = 1.0_dp/(2.0_dp/k - 1.0_dp/a)
   end function ak2w

   elemental pure function de2a(d, e) result(a)
      real(dp), intent(in) :: d, e
      real(dp) :: a
      a = e/(d*(1.0_dp-e))
   end function de2a

   elemental pure function de2k(d, e) result(k)
      real(dp), intent(in) :: d, e
      real(dp) :: k
      k = e/d
   end function de2k

   elemental pure function de2w(d, e) result(w)
      real(dp), intent(in) :: d, e
      real(dp) :: w
      w = e/(d*(1.0_dp+e))
   end function de2w

   elemental pure function dk2a(d, k) result(a)
      real(dp), intent(in) :: d, k
      real(dp) :: a
      a = 1.0_dp/(1.0_dp/k-d)
   end function dk2a

   elemental pure function dk2e(d, k) result(e)
      real(dp), intent(in) :: d, k
      real(dp) :: e
      e = d*k
   end function dk2e

   elemental pure function dk2w(d, k) result(w)
      real(dp), intent(in) :: d, k
      real(dp) :: w
      w = 1.0_dp/(1.0_dp/k+d)
   end function dk2w

   elemental pure function dw2a(d, w) result(a)
      real(dp), intent(in) :: d, w
      real(dp) :: a
      a = -1.0_dp/(2.0_dp*d-1.0_dp/w)
   end function dw2a

   elemental pure function dw2e(d, w) result(e)
      real(dp), intent(in) :: d, w
      real(dp) :: e
      e = d*w/(1.0_dp-d*w)
   end function dw2e

   elemental pure function dw2k(d, w) result(k)
      real(dp), intent(in) :: d, w
      real(dp) :: k
      k = w/(1.0_dp-d*w)
   end function dw2k

   elemental pure function ek2a(e, k) result(a)
      real(dp), intent(in) :: e, k
      real(dp) :: a
      a = k/(1.0_dp-e)
   end function ek2a

   elemental pure function ek2d(e, k) result(d)
      real(dp), intent(in) :: e, k
      real(dp) :: d
      d = e/k
   end function ek2d

   elemental pure function ek2w(e, k) result(w)
      real(dp), intent(in) :: e, k
      real(dp) :: w
      w = k/(1.0_dp+e)
   end function ek2w

   elemental pure function ew2a(e, w) result(a)
      real(dp), intent(in) :: e, w
      real(dp) :: a
      a = w*(1.0_dp+e)/(1.0_dp-e)
   end function ew2a

   elemental pure function ew2d(e, w) result(d)
      real(dp), intent(in) :: e, w
      real(dp) :: d
      d = e/((1.0_dp+e)*w)
   end function ew2d

   elemental pure function ew2k(e, w) result(k)
      real(dp), intent(in) :: e, w
      real(dp) :: k
      k = (1.0_dp+e)*w
   end function ew2k

   elemental pure function kd2a(k, d) result(a)
      real(dp), intent(in) :: k, d
      real(dp) :: a
      a = 1.0_dp/(1.0_dp/k - d)
   end function kd2a

   elemental pure function kd2e(k, d) result(e)
      real(dp), intent(in) :: k, d
      real(dp) :: e
      e = k*d
   end function kd2e

   elemental pure function kd2w(k, d) result(w)
      real(dp), intent(in) :: k, d
      real(dp) :: w
      w = 1.0_dp/(1.0_dp/k + d)
   end function kd2w

   elemental pure function ke2a(k, e) result(a)
      real(dp), intent(in) :: k, e
      real(dp) :: a
      a = k/(1.0_dp - e)
   end function ke2a

   elemental pure function ke2d(k, e) result(d)
      real(dp), intent(in) :: k, e
      real(dp) :: d
      d = e/k
   end function ke2d

   elemental pure function ke2w(k, e) result(w)
      real(dp), intent(in) :: k, e
      real(dp) :: w
      w = k/(1.0_dp + e)
   end function ke2w


   elemental pure function kw2a(k, w) result(a)
      real(dp), intent(in) :: k, w
      real(dp) :: a
      a = 1.0_dp/(2.0_dp/k - 1.0_dp/w)
   end function kw2a

   elemental pure function kw2d(k, w) result(d)
      real(dp), intent(in) :: k, w
      real(dp) :: d
      d = (k-w)/(w*k)
   end function kw2d

   elemental pure function kw2e(k, w) result(e)
      real(dp), intent(in) :: k, w
      real(dp) :: e
      e = (k-w)/w
   end function kw2e

   pure function make_k1(m, g, k) result(par)
      real(dp), intent(in) :: m, g, k
      type(kiener_parameters) :: par
      par = kiener_parameters(m=m, g=g, a=k, k=k, w=k, d=0.0_dp, e=0.0_dp)
   end function make_k1

   pure function make_k2(m, g, a, w) result(par)
      real(dp), intent(in) :: m, g, a, w
      type(kiener_parameters) :: par
      par = kiener_parameters(m=m, g=g, a=a, k=aw2k(a,w), w=w, &
                              d=aw2d(a,w), e=aw2e(a,w))
   end function make_k2

   pure function make_k3(m, g, k, d) result(par)
      real(dp), intent(in) :: m, g, k, d
      type(kiener_parameters) :: par
      par = kiener_parameters(m=m, g=g, a=kd2a(k,d), k=k, w=kd2w(k,d), &
                              d=d, e=kd2e(k,d))
   end function make_k3

   pure function make_k4(m, g, k, e) result(par)
      real(dp), intent(in) :: m, g, k, e
      type(kiener_parameters) :: par
      par = kiener_parameters(m=m, g=g, a=ke2a(k,e), k=k, w=ke2w(k,e), &
                              d=ke2d(k,e), e=e)
   end function make_k4

   elemental pure function parameters_valid(par) result(ok)
      type(kiener_parameters), intent(in) :: par
      logical :: ok
      ok = par%g > 0.0_dp .and. par%a > 0.0_dp .and. par%k > 0.0_dp .and. &
           par%w > 0.0_dp .and. abs(par%e) < 1.0_dp
   end function parameters_valid

end module fattailsr_params
