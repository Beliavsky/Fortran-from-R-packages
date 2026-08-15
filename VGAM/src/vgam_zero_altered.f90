! Copyright (C) 1998-2025 T. W. Yee, University of Auckland.
! Modern Fortran computational translation, 2026.
! SPDX-License-Identifier: GPL-3.0-only
module vgam_zero_altered
   use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
   use vgam_kinds, only : dp
   use vgam_special, only : log1p_v
   use vgam_distributions, only : dnorm_v, pnorm_v, qnorm_v, dpois_v, ppois_v, qpois_v, &
      dbinom_v, pbinom_v, qbinom_v, dnbinom_v, pnbinom_v, qnbinom_v
   implicit none
   private

   public :: dposnorm_v, pposnorm_v, qposnorm_v, rposnorm_v
   public :: dposgeom_v, pposgeom_v, qposgeom_v, rposgeom_v
   public :: dzapois_v, pzapois_v, qzapois_v, rzapois_v
   public :: dzanbinom_v, pzanbinom_v, qzanbinom_v, rzanbinom_v
   public :: dzageom_v, pzageom_v, qzageom_v, rzageom_v
   public :: dzabinom_v, pzabinom_v, qzabinom_v, rzabinom_v
   public :: dzibinom_v, pzibinom_v, qzibinom_v, rzibinom_v
   public :: dzigeom_v, pzigeom_v, qzigeom_v, rzigeom_v

contains

   elemental real(dp) function dposnorm_v(x, mean, sd, log_density) result(v)
      real(dp), intent(in) :: x, mean, sd
      logical, intent(in), optional :: log_density
      logical :: lg
      real(dp) :: norm, ld
      lg = .false.; if (present(log_density)) lg = log_density
      if (sd <= 0.0_dp) then
         v = qnan(); return
      end if
      norm = pnorm_v(mean/sd, 0.0_dp, 1.0_dp)
      if (x < 0.0_dp .or. norm <= 0.0_dp) then
         ld = -huge(1.0_dp)
      else
         ld = dnorm_v(x, mean, sd, .true.) - log(norm)
      end if
      v = merge(ld, exp(ld), lg)
   end function dposnorm_v

   elemental real(dp) function pposnorm_v(q, mean, sd) result(p)
      real(dp), intent(in) :: q, mean, sd
      real(dp) :: p0, den
      if (sd <= 0.0_dp) then
         p = qnan(); return
      end if
      if (q <= 0.0_dp) then
         p = 0.0_dp; return
      end if
      p0 = pnorm_v(0.0_dp, mean, sd); den = 1.0_dp - p0
      if (den <= 0.0_dp) then
         p = qnan()
      else
         p = min(1.0_dp, max(0.0_dp, (pnorm_v(q, mean, sd) - p0)/den))
      end if
   end function pposnorm_v

   elemental real(dp) function qposnorm_v(p, mean, sd) result(q)
      real(dp), intent(in) :: p, mean, sd
      real(dp) :: p0
      if (sd <= 0.0_dp .or. p < 0.0_dp .or. p > 1.0_dp) then
         q = qnan(); return
      end if
      p0 = pnorm_v(0.0_dp, mean, sd)
      q = qnorm_v(p0 + p*(1.0_dp - p0), mean, sd)
   end function qposnorm_v

   real(dp) function rposnorm_v(mean, sd) result(x)
      real(dp), intent(in) :: mean, sd
      real(dp) :: u
      call random_number(u); x = qposnorm_v(u, mean, sd)
   end function rposnorm_v

   elemental real(dp) function dposgeom_v(x, prob, log_density) result(v)
      integer, intent(in) :: x
      real(dp), intent(in) :: prob
      logical, intent(in), optional :: log_density
      logical :: lg
      real(dp) :: ld
      lg = .false.; if (present(log_density)) lg = log_density
      if (prob <= 0.0_dp .or. prob >= 1.0_dp .or. x < 1) then
         ld = -huge(1.0_dp)
      else
         ld = log(prob) + real(x - 1, dp)*log1p_v(-prob)
      end if
      v = merge(ld, exp(ld), lg)
   end function dposgeom_v

   elemental real(dp) function pposgeom_v(q, prob) result(p)
      integer, intent(in) :: q
      real(dp), intent(in) :: prob
      if (prob <= 0.0_dp .or. prob >= 1.0_dp) then
         p = qnan()
      else if (q < 1) then
         p = 0.0_dp
      else
         p = 1.0_dp - (1.0_dp - prob)**q
      end if
   end function pposgeom_v

   integer function qposgeom_v(p, prob) result(q)
      real(dp), intent(in) :: p, prob
      if (prob <= 0.0_dp .or. prob >= 1.0_dp .or. p <= 0.0_dp .or. p > 1.0_dp) then
         q = -1
      else if (p >= 1.0_dp) then
         q = huge(0)
      else
         q = max(1, ceiling(log1p_v(-p)/log1p_v(-prob)))
      end if
   end function qposgeom_v

   integer function rposgeom_v(prob) result(x)
      real(dp), intent(in) :: prob
      real(dp) :: u
      call random_number(u); u = max(u, epsilon(1.0_dp)); x = qposgeom_v(u, prob)
   end function rposgeom_v

   elemental real(dp) function dzapois_v(x, lambda, pobs0, log_density) result(v)
      integer, intent(in) :: x
      real(dp), intent(in) :: lambda, pobs0
      logical, intent(in), optional :: log_density
      real(dp) :: p0, d, ld
      logical :: lg
      lg = .false.; if (present(log_density)) lg = log_density
      if (lambda <= 0.0_dp .or. pobs0 < 0.0_dp .or. pobs0 > 1.0_dp .or. x < 0) then
         v = qnan(); return
      end if
      p0 = dpois_v(0, lambda)
      if (x == 0) then
         d = pobs0
      else
         d = (1.0_dp - pobs0)*dpois_v(x, lambda)/max(1.0_dp - p0, tiny(1.0_dp))
      end if
      ld = merge(log(d), -huge(1.0_dp), d > 0.0_dp); v = merge(ld, d, lg)
   end function dzapois_v

   elemental real(dp) function pzapois_v(q, lambda, pobs0) result(p)
      integer, intent(in) :: q
      real(dp), intent(in) :: lambda, pobs0
      real(dp) :: p0
      if (lambda <= 0.0_dp .or. pobs0 < 0.0_dp .or. pobs0 > 1.0_dp) then
         p = qnan(); return
      end if
      if (q < 0) then; p = 0.0_dp; return; end if
      if (q == 0) then; p = pobs0; return; end if
      p0 = dpois_v(0, lambda)
      p = pobs0 + (1.0_dp - pobs0)*(ppois_v(q, lambda) - p0)/max(1.0_dp - p0, tiny(1.0_dp))
      p = min(1.0_dp, max(0.0_dp, p))
   end function pzapois_v

   integer function qzapois_v(p, lambda, pobs0) result(q)
      real(dp), intent(in) :: p, lambda, pobs0
      real(dp) :: p0, target
      if (lambda <= 0.0_dp .or. pobs0 < 0.0_dp .or. pobs0 > 1.0_dp .or. p < 0.0_dp .or. p > 1.0_dp) then
         q = -1; return
      end if
      if (p <= pobs0) then; q = 0; return; end if
      p0 = dpois_v(0, lambda)
      target = p0 + (1.0_dp - p0)*(p - pobs0)/max(1.0_dp - pobs0, tiny(1.0_dp))
      q = max(1, qpois_v(min(1.0_dp, target), lambda))
   end function qzapois_v

   integer function rzapois_v(lambda, pobs0) result(x)
      real(dp), intent(in) :: lambda, pobs0
      real(dp) :: u
      call random_number(u); x = qzapois_v(u, lambda, pobs0)
   end function rzapois_v

   elemental real(dp) function dzanbinom_v(x, mu, sizev, pobs0, log_density) result(v)
      integer, intent(in) :: x
      real(dp), intent(in) :: mu, sizev, pobs0
      logical, intent(in), optional :: log_density
      real(dp) :: prob, p0, d, ld
      logical :: lg
      lg = .false.; if (present(log_density)) lg = log_density
      if (min(mu, sizev) <= 0.0_dp .or. pobs0 < 0.0_dp .or. pobs0 > 1.0_dp .or. x < 0) then
         v = qnan(); return
      end if
      prob = sizev/(sizev + mu); p0 = dnbinom_v(0, sizev, prob)
      if (x == 0) then; d = pobs0
      else; d = (1.0_dp - pobs0)*dnbinom_v(x, sizev, prob)/max(1.0_dp - p0, tiny(1.0_dp)); end if
      ld = merge(log(d), -huge(1.0_dp), d > 0.0_dp); v = merge(ld, d, lg)
   end function dzanbinom_v

   elemental real(dp) function pzanbinom_v(q, mu, sizev, pobs0) result(p)
      integer, intent(in) :: q
      real(dp), intent(in) :: mu, sizev, pobs0
      real(dp) :: prob, p0
      if (min(mu, sizev) <= 0.0_dp .or. pobs0 < 0.0_dp .or. pobs0 > 1.0_dp) then
         p = qnan(); return
      end if
      if (q < 0) then; p = 0.0_dp; return; end if
      if (q == 0) then; p = pobs0; return; end if
      prob = sizev/(sizev + mu); p0 = dnbinom_v(0, sizev, prob)
      p = pobs0 + (1.0_dp - pobs0)*(pnbinom_v(q, sizev, prob) - p0)/max(1.0_dp - p0, tiny(1.0_dp))
      p = min(1.0_dp, max(0.0_dp, p))
   end function pzanbinom_v

   integer function qzanbinom_v(p, mu, sizev, pobs0) result(q)
      real(dp), intent(in) :: p, mu, sizev, pobs0
      real(dp) :: prob, p0, target
      if (min(mu, sizev) <= 0.0_dp .or. pobs0 < 0.0_dp .or. pobs0 > 1.0_dp .or. p < 0.0_dp .or. p > 1.0_dp) then
         q = -1; return
      end if
      if (p <= pobs0) then; q = 0; return; end if
      prob = sizev/(sizev + mu); p0 = dnbinom_v(0, sizev, prob)
      target = p0 + (1.0_dp - p0)*(p - pobs0)/max(1.0_dp - pobs0, tiny(1.0_dp))
      q = max(1, qnbinom_v(min(1.0_dp, target), sizev, prob))
   end function qzanbinom_v

   integer function rzanbinom_v(mu, sizev, pobs0) result(x)
      real(dp), intent(in) :: mu, sizev, pobs0
      real(dp) :: u
      call random_number(u); x = qzanbinom_v(u, mu, sizev, pobs0)
   end function rzanbinom_v

   elemental real(dp) function dzageom_v(x, prob, pobs0, log_density) result(v)
      integer, intent(in) :: x
      real(dp), intent(in) :: prob, pobs0
      logical, intent(in), optional :: log_density
      real(dp) :: d, ld
      logical :: lg
      lg = .false.; if (present(log_density)) lg = log_density
      if (prob <= 0.0_dp .or. prob >= 1.0_dp .or. pobs0 < 0.0_dp .or. pobs0 > 1.0_dp .or. x < 0) then
         v = qnan(); return
      end if
      if (x == 0) then; d = pobs0; else; d = (1.0_dp - pobs0)*dposgeom_v(x, prob); end if
      ld = merge(log(d), -huge(1.0_dp), d > 0.0_dp); v = merge(ld, d, lg)
   end function dzageom_v

   elemental real(dp) function pzageom_v(q, prob, pobs0) result(p)
      integer, intent(in) :: q
      real(dp), intent(in) :: prob, pobs0
      if (prob <= 0.0_dp .or. prob >= 1.0_dp .or. pobs0 < 0.0_dp .or. pobs0 > 1.0_dp) then
         p = qnan(); return
      end if
      if (q < 0) then; p = 0.0_dp
      else if (q == 0) then; p = pobs0
      else; p = pobs0 + (1.0_dp - pobs0)*pposgeom_v(q, prob); end if
   end function pzageom_v

   integer function qzageom_v(p, prob, pobs0) result(q)
      real(dp), intent(in) :: p, prob, pobs0
      if (prob <= 0.0_dp .or. prob >= 1.0_dp .or. pobs0 < 0.0_dp .or. pobs0 > 1.0_dp .or. &
          p < 0.0_dp .or. p > 1.0_dp) then
         q = -1
      else if (p <= pobs0) then
         q = 0
      else
         q = qposgeom_v((p - pobs0)/max(1.0_dp - pobs0, tiny(1.0_dp)), prob)
      end if
   end function qzageom_v

   integer function rzageom_v(prob, pobs0) result(x)
      real(dp), intent(in) :: prob, pobs0
      real(dp) :: u
      call random_number(u); x = qzageom_v(u, prob, pobs0)
   end function rzageom_v

   elemental real(dp) function dzabinom_v(x, n, prob, pobs0, log_density) result(v)
      integer, intent(in) :: x, n
      real(dp), intent(in) :: prob, pobs0
      logical, intent(in), optional :: log_density
      real(dp) :: p0, d, ld
      logical :: lg
      lg = .false.; if (present(log_density)) lg = log_density
      if (n < 1 .or. prob <= 0.0_dp .or. prob >= 1.0_dp .or. pobs0 < 0.0_dp .or. pobs0 > 1.0_dp .or. x < 0) then
         v = qnan(); return
      end if
      p0 = dbinom_v(0, n, prob)
      if (x == 0) then; d = pobs0
      else; d = (1.0_dp - pobs0)*dbinom_v(x, n, prob)/max(1.0_dp - p0, tiny(1.0_dp)); end if
      ld = merge(log(d), -huge(1.0_dp), d > 0.0_dp); v = merge(ld, d, lg)
   end function dzabinom_v

   elemental real(dp) function pzabinom_v(q, n, prob, pobs0) result(p)
      integer, intent(in) :: q, n
      real(dp), intent(in) :: prob, pobs0
      real(dp) :: p0
      if (n < 1 .or. prob <= 0.0_dp .or. prob >= 1.0_dp .or. pobs0 < 0.0_dp .or. pobs0 > 1.0_dp) then
         p = qnan(); return
      end if
      if (q < 0) then; p = 0.0_dp; return; end if
      if (q == 0) then; p = pobs0; return; end if
      if (q >= n) then; p = 1.0_dp; return; end if
      p0 = dbinom_v(0, n, prob)
      p = pobs0 + (1.0_dp - pobs0)*(pbinom_v(q, n, prob) - p0)/max(1.0_dp - p0, tiny(1.0_dp))
      p = min(1.0_dp, max(0.0_dp, p))
   end function pzabinom_v

   integer function qzabinom_v(p, n, prob, pobs0) result(q)
      real(dp), intent(in) :: p, prob, pobs0
      integer, intent(in) :: n
      real(dp) :: p0, target
      if (n < 1 .or. prob <= 0.0_dp .or. prob >= 1.0_dp .or. pobs0 < 0.0_dp .or. pobs0 > 1.0_dp .or. &
          p < 0.0_dp .or. p > 1.0_dp) then
         q = -1; return
      end if
      if (p <= pobs0) then; q = 0; return; end if
      p0 = dbinom_v(0, n, prob)
      target = p0 + (1.0_dp - p0)*(p - pobs0)/max(1.0_dp - pobs0, tiny(1.0_dp))
      q = max(1, qbinom_v(min(1.0_dp, target), n, prob))
   end function qzabinom_v

   integer function rzabinom_v(n, prob, pobs0) result(x)
      integer, intent(in) :: n
      real(dp), intent(in) :: prob, pobs0
      real(dp) :: u
      call random_number(u); x = qzabinom_v(u, n, prob, pobs0)
   end function rzabinom_v

   elemental real(dp) function dzibinom_v(x, n, prob, pstr0, log_density) result(v)
      integer, intent(in) :: x, n
      real(dp), intent(in) :: prob, pstr0
      logical, intent(in), optional :: log_density
      real(dp) :: p0, limit, d, ld
      logical :: lg
      lg = .false.; if (present(log_density)) lg = log_density
      if (n < 0 .or. prob < 0.0_dp .or. prob > 1.0_dp .or. x < 0 .or. x > n) then
         v = qnan(); return
      end if
      p0 = dbinom_v(0, n, prob); limit = -p0/max(1.0_dp - p0, tiny(1.0_dp))
      if (pstr0 < limit .or. pstr0 > 1.0_dp) then; v = qnan(); return; end if
      if (x == 0) then; d = pstr0 + (1.0_dp - pstr0)*p0
      else; d = (1.0_dp - pstr0)*dbinom_v(x, n, prob); end if
      ld = merge(log(d), -huge(1.0_dp), d > 0.0_dp); v = merge(ld, d, lg)
   end function dzibinom_v

   elemental real(dp) function pzibinom_v(q, n, prob, pstr0) result(p)
      integer, intent(in) :: q, n
      real(dp), intent(in) :: prob, pstr0
      real(dp) :: p0, limit
      if (n < 0 .or. prob < 0.0_dp .or. prob > 1.0_dp) then; p = qnan(); return; end if
      p0 = dbinom_v(0, n, prob); limit = -p0/max(1.0_dp - p0, tiny(1.0_dp))
      if (pstr0 < limit .or. pstr0 > 1.0_dp) then; p = qnan(); return; end if
      if (q < 0) then; p = 0.0_dp
      else if (q >= n) then; p = 1.0_dp
      else; p = pstr0 + (1.0_dp - pstr0)*pbinom_v(q, n, prob); end if
      p = min(1.0_dp, max(0.0_dp, p))
   end function pzibinom_v

   integer function qzibinom_v(p, n, prob, pstr0) result(q)
      real(dp), intent(in) :: p, prob, pstr0
      integer, intent(in) :: n
      integer :: k
      if (p < 0.0_dp .or. p > 1.0_dp) then; q = -1; return; end if
      do k = 0, n
         if (pzibinom_v(k, n, prob, pstr0) >= p) then; q = k; return; end if
      end do
      q = n
   end function qzibinom_v

   integer function rzibinom_v(n, prob, pstr0) result(x)
      integer, intent(in) :: n
      real(dp), intent(in) :: prob, pstr0
      real(dp) :: u
      call random_number(u); x = qzibinom_v(u, n, prob, pstr0)
   end function rzibinom_v

   elemental real(dp) function dzigeom_v(x, prob, pstr0, log_density) result(v)
      integer, intent(in) :: x
      real(dp), intent(in) :: prob, pstr0
      logical, intent(in), optional :: log_density
      real(dp) :: p0, limit, d, ld
      logical :: lg
      lg = .false.; if (present(log_density)) lg = log_density
      if (prob <= 0.0_dp .or. prob > 1.0_dp .or. x < 0) then; v = qnan(); return; end if
      p0 = prob; limit = -p0/max(1.0_dp - p0, tiny(1.0_dp))
      if (pstr0 < limit .or. pstr0 > 1.0_dp) then; v = qnan(); return; end if
      if (x == 0) then
         d = pstr0 + (1.0_dp - pstr0)*prob
      else
         d = (1.0_dp - pstr0)*prob*(1.0_dp - prob)**x
      end if
      ld = merge(log(d), -huge(1.0_dp), d > 0.0_dp); v = merge(ld, d, lg)
   end function dzigeom_v

   elemental real(dp) function pzigeom_v(q, prob, pstr0) result(p)
      integer, intent(in) :: q
      real(dp), intent(in) :: prob, pstr0
      real(dp) :: p0, limit, basecdf
      if (prob <= 0.0_dp .or. prob > 1.0_dp) then; p = qnan(); return; end if
      p0 = prob; limit = -p0/max(1.0_dp - p0, tiny(1.0_dp))
      if (pstr0 < limit .or. pstr0 > 1.0_dp) then; p = qnan(); return; end if
      if (q < 0) then; p = 0.0_dp; return; end if
      basecdf = 1.0_dp - (1.0_dp - prob)**(q + 1)
      p = pstr0 + (1.0_dp - pstr0)*basecdf; p = min(1.0_dp, max(0.0_dp, p))
   end function pzigeom_v

   integer function qzigeom_v(p, prob, pstr0) result(q)
      real(dp), intent(in) :: p, prob, pstr0
      integer :: k
      if (p < 0.0_dp .or. p > 1.0_dp) then; q = -1; return; end if
      do k = 0, 100000
         if (pzigeom_v(k, prob, pstr0) >= p) then; q = k; return; end if
      end do
      q = -1
   end function qzigeom_v

   integer function rzigeom_v(prob, pstr0) result(x)
      real(dp), intent(in) :: prob, pstr0
      real(dp) :: u
      call random_number(u); x = qzigeom_v(u, prob, pstr0)
   end function rzigeom_v

   elemental real(dp) function qnan() result(x)
      x = ieee_value(0.0_dp, ieee_quiet_nan)
   end function qnan
end module vgam_zero_altered
