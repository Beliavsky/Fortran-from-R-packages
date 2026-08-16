! Argus distribution functions and random variate generators.
! Ported from argus 0.1.1 by Wolfgang Hoermann and Christoph Baumgarten.
! SPDX-License-Identifier: GPL-2.0-or-later
module argus_distribution
   use argus_kinds, only : dp, i8
   use argus_special, only : log_gamma_p_3half, inv_gamma_p_3half_log, &
      log1mexp, expm1_stable, nan_dp, neg_inf_dp
   implicit none
   private

   integer, parameter, public :: ARGUS_INVERSION = 1
   integer, parameter, public :: ARGUS_ROU = 2

   real(dp), parameter :: log_two_pi_half = 0.5_dp*log(2.0_dp*acos(-1.0_dp))
   real(dp), parameter :: log_half = log(0.5_dp)

   public :: dargus, pargus, qargus
   public :: dargus_recycle, pargus_recycle, qargus_recycle
   public :: rargus, rargus_varying, seed_argus_rng

contains

   elemental real(dp) function dargus(x, chi, log_pdf) result(ans)
      real(dp), intent(in) :: x, chi
      logical, intent(in), optional :: log_pdf
      logical :: want_log
      real(dp) :: logc, logf

      want_log = .false.
      if (present(log_pdf)) want_log = log_pdf

      if (chi <= 0.0_dp) then
         ans = nan_dp()
         return
      end if
      if (x <= 0.0_dp .or. x >= 1.0_dp) then
         if (want_log) then
            ans = neg_inf_dp()
         else
            ans = 0.0_dp
         end if
         return
      end if

      logc = log_gamma_p_3half(0.5_dp*chi*chi)
      logf = 3.0_dp*log(chi) - log_two_pi_half - logc - log_half + log(x) + &
         0.5_dp*(log(1.0_dp-x) + log(1.0_dp+x)) - 0.5_dp*chi*chi*(1.0_dp-x*x)

      if (want_log) then
         ans = logf
      else
         ans = exp(logf)
      end if
   end function dargus

   elemental real(dp) function pargus(x, chi, lower, log_p) result(ans)
      real(dp), intent(in) :: x, chi
      logical, intent(in), optional :: lower, log_p
      logical :: want_lower, want_log
      real(dp) :: log_upper, logc, z

      want_lower = .true.
      want_log = .false.
      if (present(lower)) want_lower = lower
      if (present(log_p)) want_log = log_p

      if (chi <= 0.0_dp) then
         ans = nan_dp()
         return
      end if

      if (x <= 0.0_dp) then
         if (want_lower) then
            if (want_log) then
               ans = neg_inf_dp()
            else
               ans = 0.0_dp
            end if
         else
            if (want_log) then
               ans = 0.0_dp
            else
               ans = 1.0_dp
            end if
         end if
         return
      else if (x >= 1.0_dp) then
         if (want_lower) then
            if (want_log) then
               ans = 0.0_dp
            else
               ans = 1.0_dp
            end if
         else
            if (want_log) then
               ans = neg_inf_dp()
            else
               ans = 0.0_dp
            end if
         end if
         return
      end if

      logc = log_gamma_p_3half(0.5_dp*chi*chi)
      z = 0.5_dp*chi*chi*(1.0_dp-x*x)
      log_upper = log_gamma_p_3half(z) - logc
      if (log_upper > 0.0_dp) log_upper = 0.0_dp

      if (.not. want_lower) then
         if (want_log) then
            ans = log_upper
         else
            ans = exp(log_upper)
         end if
      else
         if (want_log) then
            ans = log1mexp(log_upper)
         else
            ans = -expm1_stable(log_upper)
         end if
      end if
   end function pargus

   elemental real(dp) function qargus(p, chi, lower, log_p) result(ans)
      real(dp), intent(in) :: p, chi
      logical, intent(in), optional :: lower, log_p
      logical :: want_lower, input_log
      real(dp) :: logc, logtarget, logprob, zmax, y

      want_lower = .true.
      input_log = .false.
      if (present(lower)) want_lower = lower
      if (present(log_p)) input_log = log_p

      if (chi <= 0.0_dp) then
         ans = nan_dp()
         return
      end if

      if (input_log) then
         if (p > 0.0_dp) then
            ans = nan_dp()
            return
         end if
         logprob = p
      else
         if (p < 0.0_dp .or. p > 1.0_dp) then
            ans = nan_dp()
            return
         end if
         if (p <= 0.0_dp) then
            logprob = neg_inf_dp()
         else
            logprob = log(p)
         end if
      end if

      if (want_lower) then
         if ((.not. input_log .and. p <= 0.0_dp) .or. &
             (input_log .and. logprob < -huge(1.0_dp))) then
            ans = 0.0_dp
            return
         else if ((.not. input_log .and. p >= 1.0_dp) .or. &
                  (input_log .and. logprob >= 0.0_dp)) then
            ans = 1.0_dp
            return
         end if
      else
         if ((.not. input_log .and. p <= 0.0_dp) .or. &
             (input_log .and. logprob < -huge(1.0_dp))) then
            ans = 1.0_dp
            return
         else if ((.not. input_log .and. p >= 1.0_dp) .or. &
                  (input_log .and. logprob >= 0.0_dp)) then
            ans = 0.0_dp
            return
         end if
      end if

      zmax = 0.5_dp*chi*chi
      logc = log_gamma_p_3half(zmax)
      if (want_lower) then
         if (input_log) then
            logtarget = logc + log1mexp(logprob)
         else
            logtarget = logc + log1mexp(logprob)
         end if
      else
         logtarget = logc + logprob
      end if

      y = inv_gamma_p_3half_log(logtarget, zmax)
      ans = sqrt(max(0.0_dp, 1.0_dp - y/zmax))
   end function qargus

   subroutine dargus_recycle(x, chi, ans, log_pdf)
      real(dp), intent(in) :: x(:), chi(:)
      real(dp), intent(out) :: ans(:)
      logical, intent(in), optional :: log_pdf
      integer :: i, nx, nc

      nx = size(x); nc = size(chi)
      if (size(ans) /= max(nx,nc) .or. nx == 0 .or. nc == 0) error stop "dargus_recycle: size mismatch"
      do i = 1, size(ans)
         ans(i) = dargus(x(1+mod(i-1,nx)), chi(1+mod(i-1,nc)), log_pdf)
      end do
   end subroutine dargus_recycle

   subroutine pargus_recycle(x, chi, ans, lower, log_p)
      real(dp), intent(in) :: x(:), chi(:)
      real(dp), intent(out) :: ans(:)
      logical, intent(in), optional :: lower, log_p
      integer :: i, nx, nc

      nx = size(x); nc = size(chi)
      if (size(ans) /= max(nx,nc) .or. nx == 0 .or. nc == 0) error stop "pargus_recycle: size mismatch"
      do i = 1, size(ans)
         ans(i) = pargus(x(1+mod(i-1,nx)), chi(1+mod(i-1,nc)), lower, log_p)
      end do
   end subroutine pargus_recycle

   subroutine qargus_recycle(p, chi, ans, lower, log_p)
      real(dp), intent(in) :: p(:), chi(:)
      real(dp), intent(out) :: ans(:)
      logical, intent(in), optional :: lower, log_p
      integer :: i, np, nc

      np = size(p); nc = size(chi)
      if (size(ans) /= max(np,nc) .or. np == 0 .or. nc == 0) error stop "qargus_recycle: size mismatch"
      do i = 1, size(ans)
         ans(i) = qargus(p(1+mod(i-1,np)), chi(1+mod(i-1,nc)), lower, log_p)
      end do
   end subroutine qargus_recycle

   subroutine rargus(sample, chi, method)
      real(dp), intent(out) :: sample(:)
      real(dp), intent(in) :: chi
      integer, intent(in), optional :: method
      integer :: meth, i
      real(dp) :: u

      if (chi <= 0.0_dp) error stop "rargus: chi must be positive"
      meth = ARGUS_INVERSION
      if (present(method)) meth = method

      select case (meth)
      case (ARGUS_INVERSION)
         do i = 1, size(sample)
            call random_number(u)
            sample(i) = qargus(u, chi)
         end do
      case (ARGUS_ROU)
         do i = 1, size(sample)
            sample(i) = rargus_rou_one(chi)
         end do
      case default
         error stop "rargus: unknown method"
      end select
   end subroutine rargus

   subroutine rargus_varying(sample, chi, method)
      real(dp), intent(out) :: sample(:)
      real(dp), intent(in) :: chi(:)
      integer, intent(in), optional :: method
      integer :: meth, i, nc
      real(dp) :: u, c

      nc = size(chi)
      if (nc == 0) error stop "rargus_varying: chi must not be empty"
      if (any(chi <= 0.0_dp)) error stop "rargus_varying: chi must be positive"
      meth = ARGUS_INVERSION
      if (present(method)) meth = method

      do i = 1, size(sample)
         c = chi(1+mod(i-1,nc))
         select case (meth)
         case (ARGUS_INVERSION)
            call random_number(u)
            sample(i) = qargus(u, c)
         case (ARGUS_ROU)
            sample(i) = rargus_rou_one(c)
         case default
            error stop "rargus_varying: unknown method"
         end select
      end do
   end subroutine rargus_varying

   real(dp) function rargus_rou_one(chi) result(x)
      real(dp), intent(in) :: chi
      real(dp) :: m, b, ap, am, xp, xm, u, v, y, r1, r2, t

      t = chi*chi
      if (chi <= 1.0_dp) then
         m = 0.5_dp*t
         b = sqrt(sqrt(0.5_dp)*chi)*exp(-0.25_dp*t)
         ap = 0.0_dp
         ! Algebraically equivalent to the upstream expression but stable
         ! for very small chi (avoids subtracting nearly equal numbers).
         xm = t/(t + 5.0_dp + sqrt(t*(t+6.0_dp) + 25.0_dp))
         am = (xm-0.5_dp*t)*sqrt(sqrt(xm)*exp(-xm))
      else
         m = 0.5_dp
         b = exp(-(1.0_dp+log(2.0_dp))/4.0_dp)
         xp = min(1.5_dp+sqrt(2.0_dp), 0.5_dp*t)
         xm = 1.5_dp-sqrt(2.0_dp)
         ap = (xp-0.5_dp)*sqrt(sqrt(xp)*exp(-xp))
         am = (xm-0.5_dp)*sqrt(sqrt(xm)*exp(-xm))
      end if

      do
         call random_number(r1)
         call random_number(r2)
         u = b*r1
         if (u <= 0.0_dp) cycle
         v = am + r2*(ap-am)
         y = v/u + m
         if (y < 0.0_dp .or. y > 0.5_dp*t) cycle
         if (u*u > sqrt(y)*exp(-y)) cycle
         x = sqrt(max(0.0_dp, 1.0_dp-2.0_dp*y/t))
         return
      end do
   end function rargus_rou_one

   subroutine seed_argus_rng(seed)
      integer, intent(in) :: seed
      integer, allocatable :: put(:)
      integer :: i, n
      integer(i8) :: s

      call random_seed(size=n)
      allocate(put(n))
      s = int(seed,kind=i8)
      if (s == 0_i8) s = 104729_i8
      do i = 1, n
         s = modulo(1664525_i8*s + 1013904223_i8 + 69069_i8*int(i,i8), 2147483647_i8)
         put(i) = int(max(1_i8,s),kind(put))
      end do
      call random_seed(put=put)
   end subroutine seed_argus_rng

end module argus_distribution
