! FatTailsR modern Fortran translation
! Copyright (C) 2014-2026 Patrice Kiener
! Licensed under GPL-2.0-only. See COPYING.
module fattailsr_distributions
   use fattailsr_kinds, only : dp
   use fattailsr_math, only : pi, sqrt3, kiener_scale, logit, invlogit, &
      clamp_probability, incomplete_beta
   use fattailsr_params, only : kiener_parameters, make_k1, make_k2, make_k3, make_k4
   implicit none
   private

   public :: dlogisst, plogisst, qlogisst, rlogisst, dplogisst, dqlogisst
   public :: llogisst, dllogisst, qllogisst, ltmlogisst, rtmlogisst, eslogisst
   public :: qkiener, pkiener, dkiener, rkiener, qlkiener, lkiener
   public :: dpkiener, dqkiener, dlkiener
   public :: ltmkiener, rtmkiener, dtmlqkiener, eskiener, varkiener
   public :: ckiener, hkiener
   public :: qkiener1, qkiener2, qkiener3, qkiener4
   public :: pkiener1, pkiener2, pkiener3, pkiener4
   public :: dkiener1, dkiener2, dkiener3, dkiener4

contains

   elemental pure function safe_exp(x) result(v)
      real(dp), intent(in) :: x
      real(dp) :: v
      v = exp(min(max(x, -700.0_dp), 700.0_dp))
   end function safe_exp

   elemental pure function dlogisst(x, m, g) result(v)
      real(dp), intent(in) :: x
      real(dp), intent(in), optional :: m, g
      real(dp) :: v, mu, gamma, z, ez
      mu = 0.0_dp
      gamma = 1.0_dp
      if (present(m)) mu = m
      if (present(g)) gamma = g
      z = (x - mu)/(gamma*kiener_scale)
      if (z >= 0.0_dp) then
         ez = exp(-z)
         v = ez/((1.0_dp + ez)**2*gamma*kiener_scale)
      else
         ez = exp(z)
         v = ez/((1.0_dp + ez)**2*gamma*kiener_scale)
      end if
   end function dlogisst

   elemental pure function plogisst(x, m, g) result(v)
      real(dp), intent(in) :: x
      real(dp), intent(in), optional :: m, g
      real(dp) :: v, mu, gamma
      mu = 0.0_dp
      gamma = 1.0_dp
      if (present(m)) mu = m
      if (present(g)) gamma = g
      v = invlogit((x - mu)/(gamma*kiener_scale))
   end function plogisst

   elemental pure function qlogisst(p, m, g) result(v)
      real(dp), intent(in) :: p
      real(dp), intent(in), optional :: m, g
      real(dp) :: v, mu, gamma
      mu = 0.0_dp
      gamma = 1.0_dp
      if (present(m)) mu = m
      if (present(g)) gamma = g
      v = mu + gamma*kiener_scale*logit(p)
   end function qlogisst

   subroutine rlogisst(x, m, g)
      real(dp), intent(out) :: x(:)
      real(dp), intent(in), optional :: m, g
      real(dp) :: u(size(x)), mu, gamma
      mu = 0.0_dp
      gamma = 1.0_dp
      if (present(m)) mu = m
      if (present(g)) gamma = g
      call random_number(u)
      x = qlogisst(u, mu, gamma)
   end subroutine rlogisst

   elemental pure function dplogisst(p, m, g) result(v)
      real(dp), intent(in) :: p
      real(dp), intent(in), optional :: m, g
      real(dp) :: v, gamma
      gamma = 1.0_dp
      if (present(g)) gamma = g
      v = clamp_probability(p)*(1.0_dp - clamp_probability(p))/(gamma*kiener_scale)
      if (present(m)) v = v + 0.0_dp*m
   end function dplogisst

   elemental pure function dqlogisst(p, m, g) result(v)
      real(dp), intent(in) :: p
      real(dp), intent(in), optional :: m, g
      real(dp) :: v, gamma, q
      gamma = 1.0_dp
      if (present(g)) gamma = g
      q = clamp_probability(p)
      v = gamma*kiener_scale/(q*(1.0_dp - q))
      if (present(m)) v = v + 0.0_dp*m
   end function dqlogisst

   elemental pure function llogisst(x, m, g) result(v)
      real(dp), intent(in) :: x
      real(dp), intent(in), optional :: m, g
      real(dp) :: v, mu, gamma
      mu = 0.0_dp
      gamma = 1.0_dp
      if (present(m)) mu = m
      if (present(g)) gamma = g
      v = (x - mu)/(gamma*kiener_scale)
   end function llogisst

   elemental pure function dllogisst(lp, m, g) result(v)
      real(dp), intent(in) :: lp
      real(dp), intent(in), optional :: m, g
      real(dp) :: v, gamma, p
      gamma = 1.0_dp
      if (present(g)) gamma = g
      p = invlogit(lp)
      v = p*(1.0_dp - p)/(gamma*kiener_scale)
      if (present(m)) v = v + 0.0_dp*m
   end function dllogisst

   elemental pure function qllogisst(lp, m, g) result(v)
      real(dp), intent(in) :: lp
      real(dp), intent(in), optional :: m, g
      real(dp) :: v, mu, gamma
      mu = 0.0_dp
      gamma = 1.0_dp
      if (present(m)) mu = m
      if (present(g)) gamma = g
      v = mu + gamma*kiener_scale*lp
   end function qllogisst

   elemental pure function ltmlogisst(p, m, g) result(v)
      real(dp), intent(in) :: p
      real(dp), intent(in), optional :: m, g
      real(dp) :: v, mu, gamma, q
      mu = 0.0_dp
      gamma = 1.0_dp
      if (present(m)) mu = m
      if (present(g)) gamma = g
      q = clamp_probability(p)
      v = mu + gamma*kiener_scale/q*((1.0_dp-q)*log(1.0_dp-q) + q*log(q))
   end function ltmlogisst

   elemental pure function rtmlogisst(p, m, g) result(v)
      real(dp), intent(in) :: p
      real(dp), intent(in), optional :: m, g
      real(dp) :: v, mu, gamma, q
      mu = 0.0_dp
      gamma = 1.0_dp
      if (present(m)) mu = m
      if (present(g)) gamma = g
      q = clamp_probability(p)
      v = mu - gamma*kiener_scale/(1.0_dp-q)*&
          ((1.0_dp-q)*log(1.0_dp-q) + q*log(q))
   end function rtmlogisst

   elemental pure function eslogisst(p, m, g, signed_es) result(v)
      real(dp), intent(in) :: p
      real(dp), intent(in), optional :: m, g
      logical, intent(in), optional :: signed_es
      real(dp) :: v
      logical :: signed
      signed = .false.
      if (present(signed_es)) signed = signed_es
      if (p <= 0.5_dp) then
         v = ltmlogisst(p, m, g)
      else
         v = rtmlogisst(p, m, g)
      end if
      if (.not. signed) v = abs(v)
   end function eslogisst

   elemental pure function qlkiener(lp, par) result(v)
      real(dp), intent(in) :: lp
      type(kiener_parameters), intent(in) :: par
      real(dp) :: v
      v = par%m + kiener_scale*par%g*par%k*&
          (-safe_exp(-lp/par%a) + safe_exp(lp/par%w))/2.0_dp
   end function qlkiener

   elemental pure function qkiener(p, par) result(v)
      real(dp), intent(in) :: p
      type(kiener_parameters), intent(in) :: par
      real(dp) :: v
      v = qlkiener(logit(p), par)
   end function qkiener

   elemental pure function dq_dlp(lp, par) result(v)
      real(dp), intent(in) :: lp
      type(kiener_parameters), intent(in) :: par
      real(dp) :: v
      v = kiener_scale*par%g*par%k*&
          (safe_exp(-lp/par%a)/par%a + safe_exp(lp/par%w)/par%w)/2.0_dp
   end function dq_dlp

   pure function lkiener(x, par) result(lp)
      real(dp), intent(in) :: x
      type(kiener_parameters), intent(in) :: par
      real(dp) :: lp
      real(dp) :: lo, hi, mid, f, deriv, trial
      integer :: iter

      if (abs(x - par%m) <= tiny(1.0_dp)) then
         lp = 0.0_dp
         return
      end if
      lo = -1.0_dp
      hi = 1.0_dp
      do while (qlkiener(lo, par) > x .and. lo > -1.0e6_dp)
         lo = 2.0_dp*lo
      end do
      do while (qlkiener(hi, par) < x .and. hi < 1.0e6_dp)
         hi = 2.0_dp*hi
      end do
      lp = min(max((x - par%m)/(par%g*kiener_scale), lo), hi)
      do iter = 1, 100
         f = qlkiener(lp, par) - x
         if (abs(f) <= 32.0_dp*epsilon(1.0_dp)*max(1.0_dp, abs(x))) exit
         if (f > 0.0_dp) then
            hi = lp
         else
            lo = lp
         end if
         deriv = dq_dlp(lp, par)
         trial = lp - f/max(deriv, tiny(1.0_dp))
         if (trial <= lo .or. trial >= hi) trial = 0.5_dp*(lo + hi)
         lp = trial
      end do
      mid = 0.5_dp*(lo + hi)
      if (abs(qlkiener(mid, par) - x) < abs(qlkiener(lp, par) - x)) lp = mid
   end function lkiener

   pure function pkiener(x, par) result(v)
      real(dp), intent(in) :: x
      type(kiener_parameters), intent(in) :: par
      real(dp) :: v
      v = invlogit(lkiener(x, par))
   end function pkiener

   elemental pure function dlkiener(lp, par) result(v)
      real(dp), intent(in) :: lp
      type(kiener_parameters), intent(in) :: par
      real(dp) :: v, p
      p = invlogit(lp)
      v = p*(1.0_dp - p)/dq_dlp(lp, par)
   end function dlkiener

   pure function dkiener(x, par) result(v)
      real(dp), intent(in) :: x
      type(kiener_parameters), intent(in) :: par
      real(dp) :: v
      v = dlkiener(lkiener(x, par), par)
   end function dkiener

   subroutine rkiener(x, par)
      real(dp), intent(out) :: x(:)
      type(kiener_parameters), intent(in) :: par
      real(dp) :: u(size(x))
      call random_number(u)
      x = qkiener(u, par)
   end subroutine rkiener

   elemental pure function dpkiener(p, par) result(v)
      real(dp), intent(in) :: p
      type(kiener_parameters), intent(in) :: par
      real(dp) :: v, q
      q = clamp_probability(p)
      v = q*(1.0_dp-q)/dq_dlp(logit(q), par)
   end function dpkiener

   elemental pure function dqkiener(p, par) result(v)
      real(dp), intent(in) :: p
      type(kiener_parameters), intent(in) :: par
      real(dp) :: v, q
      q = clamp_probability(p)
      v = dq_dlp(logit(q), par)/(q*(1.0_dp-q))
   end function dqkiener

   elemental pure function ltmkiener(p, par) result(v)
      real(dp), intent(in) :: p
      type(kiener_parameters), intent(in) :: par
      real(dp) :: v, q, b1, b2
      q = clamp_probability(p)
      b1 = incomplete_beta(q, 1.0_dp - 1.0_dp/par%a, 1.0_dp + 1.0_dp/par%a)
      b2 = incomplete_beta(q, 1.0_dp + 1.0_dp/par%w, 1.0_dp - 1.0_dp/par%w)
      v = par%m + kiener_scale*par%g*par%k/(2.0_dp*q)*(-b1 + b2)
   end function ltmkiener

   elemental pure function rtmkiener(p, par) result(v)
      real(dp), intent(in) :: p
      type(kiener_parameters), intent(in) :: par
      real(dp) :: v, q, b1, b2, u
      q = clamp_probability(p)
      u = 1.0_dp - q
      b1 = incomplete_beta(u, 1.0_dp + 1.0_dp/par%a, 1.0_dp - 1.0_dp/par%a)
      b2 = incomplete_beta(u, 1.0_dp - 1.0_dp/par%w, 1.0_dp + 1.0_dp/par%w)
      v = par%m + kiener_scale*par%g*par%k/(2.0_dp*u)*(-b1 + b2)
   end function rtmkiener

   elemental pure function dtmlqkiener(p, par) result(v)
      real(dp), intent(in) :: p
      type(kiener_parameters), intent(in) :: par
      real(dp) :: v
      if (p <= 0.5_dp) then
         v = ltmkiener(p, par) - qkiener(p, par)
      else
         v = rtmkiener(p, par) - qkiener(p, par)
      end if
   end function dtmlqkiener

   elemental pure function eskiener(p, par, signed_es) result(v)
      real(dp), intent(in) :: p
      type(kiener_parameters), intent(in) :: par
      logical, intent(in), optional :: signed_es
      real(dp) :: v
      logical :: signed
      signed = .false.
      if (present(signed_es)) signed = signed_es
      if (p <= 0.5_dp) then
         v = ltmkiener(p, par)
      else
         v = rtmkiener(p, par)
      end if
      if (.not. signed) v = abs(v)
   end function eskiener

   elemental pure function varkiener(p, par) result(v)
      real(dp), intent(in) :: p
      type(kiener_parameters), intent(in) :: par
      real(dp) :: v
      if (p <= 0.5_dp) then
         v = -qkiener(p, par)
      else
         v = qkiener(p, par)
      end if
   end function varkiener

   elemental pure function ckiener(p, par) result(v)
      real(dp), intent(in) :: p
      type(kiener_parameters), intent(in) :: par
      real(dp) :: v, lp
      lp = logit(p)
      if (abs(lp) <= sqrt(epsilon(1.0_dp))) then
         v = 1.0_dp
      else
         v = par%k/lp*sinh(lp/par%k)*safe_exp(par%e*lp/par%k)
      end if
   end function ckiener

   elemental pure function hkiener(p, par) result(v)
      real(dp), intent(in) :: p
      type(kiener_parameters), intent(in) :: par
      real(dp) :: v
      if (p <= 0.5_dp) then
         v = (ltmkiener(p, par) - par%m)/(ltmlogisst(p, par%m, par%g) - par%m)
      else
         v = (rtmkiener(p, par) - par%m)/(rtmlogisst(p, par%m, par%g) - par%m)
      end if
   end function hkiener

   elemental pure function qkiener1(p, m, g, k) result(v)
      real(dp), intent(in) :: p, m, g, k
      real(dp) :: v
      v = qkiener(p, make_k1(m,g,k))
   end function qkiener1

   elemental pure function qkiener2(p, m, g, a, w) result(v)
      real(dp), intent(in) :: p, m, g, a, w
      real(dp) :: v
      v = qkiener(p, make_k2(m,g,a,w))
   end function qkiener2

   elemental pure function qkiener3(p, m, g, k, d) result(v)
      real(dp), intent(in) :: p, m, g, k, d
      real(dp) :: v
      v = qkiener(p, make_k3(m,g,k,d))
   end function qkiener3

   elemental pure function qkiener4(p, m, g, k, e) result(v)
      real(dp), intent(in) :: p, m, g, k, e
      real(dp) :: v
      v = qkiener(p, make_k4(m,g,k,e))
   end function qkiener4

   pure function pkiener1(x, m, g, k) result(v)
      real(dp), intent(in) :: x, m, g, k
      real(dp) :: v
      v = pkiener(x, make_k1(m,g,k))
   end function pkiener1

   pure function pkiener2(x, m, g, a, w) result(v)
      real(dp), intent(in) :: x, m, g, a, w
      real(dp) :: v
      v = pkiener(x, make_k2(m,g,a,w))
   end function pkiener2

   pure function pkiener3(x, m, g, k, d) result(v)
      real(dp), intent(in) :: x, m, g, k, d
      real(dp) :: v
      v = pkiener(x, make_k3(m,g,k,d))
   end function pkiener3

   pure function pkiener4(x, m, g, k, e) result(v)
      real(dp), intent(in) :: x, m, g, k, e
      real(dp) :: v
      v = pkiener(x, make_k4(m,g,k,e))
   end function pkiener4

   pure function dkiener1(x, m, g, k) result(v)
      real(dp), intent(in) :: x, m, g, k
      real(dp) :: v
      v = dkiener(x, make_k1(m,g,k))
   end function dkiener1

   pure function dkiener2(x, m, g, a, w) result(v)
      real(dp), intent(in) :: x, m, g, a, w
      real(dp) :: v
      v = dkiener(x, make_k2(m,g,a,w))
   end function dkiener2

   pure function dkiener3(x, m, g, k, d) result(v)
      real(dp), intent(in) :: x, m, g, k, d
      real(dp) :: v
      v = dkiener(x, make_k3(m,g,k,d))
   end function dkiener3

   pure function dkiener4(x, m, g, k, e) result(v)
      real(dp), intent(in) :: x, m, g, k, e
      real(dp) :: v
      v = dkiener(x, make_k4(m,g,k,e))
   end function dkiener4

end module fattailsr_distributions
