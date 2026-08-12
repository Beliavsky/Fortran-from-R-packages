! SPDX-License-Identifier: GPL-2.0-or-later
module ceoptim_sampling
   use ceoptim_kinds, only : dp
   use ceoptim_rng, only : rng_state, rng_uniform, rng_normal, rng_gamma
   use ceoptim_linalg, only : covariance_factor, symmetric_pinv
   implicit none
   private

   real(dp), parameter :: sqrt2 = 1.4142135623730950488016887242096981_dp

   type, public :: tmvn_result
      real(dp), allocatable :: x(:, :)
      real(dp) :: rho = 1.0_dp
      integer :: nar = 0
      integer :: ngibbs = 0
      integer :: status = 0
      character(len=:), allocatable :: message
   end type tmvn_result

   public :: rtmvnorm, dirichlet_rand, normal_cdf, normal_quantile, truncated_normal

contains

   pure function normal_cdf(x) result(p)
      real(dp), intent(in) :: x
      real(dp) :: p
      p = 0.5_dp * erfc(-x / sqrt2)
   end function normal_cdf

   pure function normal_quantile(p) result(x)
      real(dp), intent(in) :: p
      real(dp) :: x
      real(dp) :: q, r
      real(dp), parameter :: a1 = -3.969683028665376e+01_dp
      real(dp), parameter :: a2 =  2.209460984245205e+02_dp
      real(dp), parameter :: a3 = -2.759285104469687e+02_dp
      real(dp), parameter :: a4 =  1.383577518672690e+02_dp
      real(dp), parameter :: a5 = -3.066479806614716e+01_dp
      real(dp), parameter :: a6 =  2.506628277459239e+00_dp
      real(dp), parameter :: b1 = -5.447609879822406e+01_dp
      real(dp), parameter :: b2 =  1.615858368580409e+02_dp
      real(dp), parameter :: b3 = -1.556989798598866e+02_dp
      real(dp), parameter :: b4 =  6.680131188771972e+01_dp
      real(dp), parameter :: b5 = -1.328068155288572e+01_dp
      real(dp), parameter :: c1 = -7.784894002430293e-03_dp
      real(dp), parameter :: c2 = -3.223964580411365e-01_dp
      real(dp), parameter :: c3 = -2.400758277161838e+00_dp
      real(dp), parameter :: c4 = -2.549732539343734e+00_dp
      real(dp), parameter :: c5 =  4.374664141464968e+00_dp
      real(dp), parameter :: c6 =  2.938163982698783e+00_dp
      real(dp), parameter :: d1 =  7.784695709041462e-03_dp
      real(dp), parameter :: d2 =  3.224671290700398e-01_dp
      real(dp), parameter :: d3 =  2.445134137142996e+00_dp
      real(dp), parameter :: d4 =  3.754408661907416e+00_dp
      real(dp), parameter :: plow = 0.02425_dp
      real(dp), parameter :: phigh = 1.0_dp - plow

      if (p <= 0.0_dp) then
         x = -huge(1.0_dp)
      else if (p >= 1.0_dp) then
         x = huge(1.0_dp)
      else if (p < plow) then
         q = sqrt(-2.0_dp * log(p))
         x = (((((c1*q + c2)*q + c3)*q + c4)*q + c5)*q + c6) / &
             ((((d1*q + d2)*q + d3)*q + d4)*q + 1.0_dp)
      else if (p <= phigh) then
         q = p - 0.5_dp
         r = q*q
         x = (((((a1*r + a2)*r + a3)*r + a4)*r + a5)*r + a6) * q / &
             (((((b1*r + b2)*r + b3)*r + b4)*r + b5)*r + 1.0_dp)
      else
         q = sqrt(-2.0_dp * log(1.0_dp - p))
         x = -(((((c1*q + c2)*q + c3)*q + c4)*q + c5)*q + c6) / &
              ((((d1*q + d2)*q + d3)*q + d4)*q + 1.0_dp)
      end if
   end function normal_quantile

   function truncated_normal(rng, mean, sd, lower, upper) result(x)
      type(rng_state), intent(inout) :: rng
      real(dp), intent(in) :: mean, sd, lower, upper
      real(dp) :: x
      real(dp) :: a, b, z

      if (sd <= 0.0_dp .or. lower > upper) then
         x = mean
         return
      end if
      a = (lower - mean) / sd
      b = (upper - mean) / sd
      z = trunc_standard(rng, a, b)
      x = mean + sd * z
   end function truncated_normal

   recursive function trunc_standard(rng, a, b) result(z)
      type(rng_state), intent(inout) :: rng
      real(dp), intent(in) :: a, b
      real(dp) :: z
      real(dp) :: alpha, pa, pb, u, u2, z0, log_ratio, mode_sq
      logical :: afinite, bfinite

      afinite = abs(a) < 0.5_dp * huge(1.0_dp)
      bfinite = abs(b) < 0.5_dp * huge(1.0_dp)

      if (a >= b) then
         z = a
         return
      end if

      if (afinite .and. a > 0.45_dp) then
         alpha = 0.5_dp * (a + sqrt(a*a + 4.0_dp))
         do
            z0 = a - log(max(rng_uniform(rng), tiny(1.0_dp))) / alpha
            if (bfinite .and. z0 > b) cycle
            u2 = rng_uniform(rng)
            if (u2 <= exp(-0.5_dp * (z0 - alpha)**2)) exit
         end do
         z = z0
         return
      end if

      if (bfinite .and. b < -0.45_dp) then
         z = -trunc_standard(rng, -b, -a)
         return
      end if

      if (afinite) then
         pa = normal_cdf(a)
      else
         pa = 0.0_dp
      end if
      if (bfinite) then
         pb = normal_cdf(b)
      else
         pb = 1.0_dp
      end if

      if (pb > pa .and. pb - pa > 32.0_dp * epsilon(1.0_dp) * max(1.0_dp, pa, pb)) then
         u = pa + rng_uniform(rng) * (pb - pa)
         u = min(max(u, tiny(1.0_dp)), 1.0_dp - epsilon(1.0_dp))
         z = normal_quantile(u)
         return
      end if

      if (afinite .and. bfinite) then
         if (a <= 0.0_dp .and. b >= 0.0_dp) then
            mode_sq = 0.0_dp
         else
            mode_sq = min(a*a, b*b)
         end if
         do
            z0 = a + (b - a) * rng_uniform(rng)
            log_ratio = -0.5_dp * (z0*z0 - mode_sq)
            if (log(max(rng_uniform(rng), tiny(1.0_dp))) <= log_ratio) exit
         end do
         z = z0
      else
         do
            z0 = rng_normal(rng)
            if ((.not. afinite .or. z0 >= a) .and. (.not. bfinite .or. z0 <= b)) exit
         end do
         z = z0
      end if
   end function trunc_standard

   subroutine dirichlet_rand(alpha, n, rng, x, status)
      real(dp), intent(in) :: alpha(:)
      integer, intent(in) :: n
      type(rng_state), intent(inout) :: rng
      real(dp), allocatable, intent(out) :: x(:, :)
      integer, intent(out), optional :: status
      integer :: i, j, p
      real(dp) :: total

      p = size(alpha)
      allocate(x(n, p))
      if (n < 0 .or. p < 1 .or. any(alpha <= 0.0_dp)) then
         x = 0.0_dp
         if (present(status)) status = 1
         return
      end if
      do i = 1, n
         total = 0.0_dp
         do j = 1, p
            x(i, j) = rng_gamma(rng, alpha(j))
            total = total + x(i, j)
         end do
         x(i, :) = x(i, :) / total
      end do
      if (present(status)) status = 0
   end subroutine dirichlet_rand

   subroutine rtmvnorm(n, mu, sigma, rng, out, a, b, rho_thr, max_sample)
      integer, intent(in) :: n
      real(dp), intent(in) :: mu(:), sigma(:, :)
      type(rng_state), intent(inout) :: rng
      type(tmvn_result), intent(out) :: out
      real(dp), intent(in), optional :: a(:, :), b(:)
      real(dp), intent(in), optional :: rho_thr
      integer, intent(in), optional :: max_sample

      real(dp), allocatable :: factor(:, :), z(:), draw(:), precision(:, :)
      real(dp) :: rt, rho, lower, upper, rhs, ai, conditional_mean, conditional_var
      integer :: p, m, info, maxs, batch, j, i, k, stored, valid_total, trials, nr
      logical :: constrained, ok, impossible

      p = size(mu)
      allocate(out%x(max(0, n), p))
      out%x = 0.0_dp
      out%rho = 1.0_dp
      out%nar = 0
      out%ngibbs = 0
      out%status = 0
      out%message = 'ok'

      if (n < 0 .or. p < 1 .or. size(sigma, 1) /= p .or. size(sigma, 2) /= p) then
         out%status = 1
         out%message = 'invalid dimensions in rtmvnorm'
         return
      end if
      constrained = present(a) .or. present(b)
      if (present(a) .neqv. present(b)) then
         out%status = 2
         out%message = 'A and b must be supplied together'
         return
      end if
      if (constrained) then
         if (size(a, 2) /= p .or. size(a, 1) /= size(b)) then
            out%status = 3
            out%message = 'A and b are not conformable'
            return
         end if
         m = size(a, 1)
      else
         m = 0
      end if

      call covariance_factor(sigma, factor, info)
      if (info /= 0) then
         out%status = 4
         out%message = 'covariance matrix is not positive semidefinite'
         return
      end if
      allocate(z(p), draw(p))

      if (.not. constrained .or. m == 0) then
         do i = 1, n
            call draw_mvn(mu, factor, rng, z, draw)
            out%x(i, :) = draw
         end do
         out%nar = n
         out%rho = 1.0_dp
         return
      end if

      rt = 1.0e-4_dp
      if (present(rho_thr)) then
         if (rho_thr >= 0.0_dp) rt = rho_thr
      end if
      maxs = 1000000
      if (present(max_sample)) then
         if (max_sample >= 0) maxs = max_sample
      end if

      stored = 0
      valid_total = 0
      trials = 0
      rho = 1.0_dp
      batch = max(1, n)

      if (rt < 1.0_dp) then
         do while (stored < n .and. (rho > rt .or. batch < maxs))
            nr = 0
            do j = 1, batch
               call draw_mvn(mu, factor, rng, z, draw)
               ok = feasible(draw, a, b)
               if (ok) then
                  nr = nr + 1
                  if (stored < n) then
                     stored = stored + 1
                     out%x(stored, :) = draw
                  end if
               end if
            end do
            valid_total = valid_total + nr
            trials = trials + batch
            rho = real(valid_total, dp) / real(trials, dp)
            if (stored >= n) exit
            if (rho > 0.0_dp) then
               batch = min(maxs, max(1, ceiling(real(n - stored, dp) / rho)), 10 * batch)
            else
               batch = min(maxs, 10 * batch)
            end if
         end do
      end if
      out%nar = stored
      out%rho = rho

      if (stored >= n) return

      if (stored > 0) then
         draw = out%x(stored, :)
      else
         draw = mu
      end if

      call symmetric_pinv(sigma, precision, info)
      if (info /= 0) then
         out%status = 5
         out%message = 'failed to compute covariance pseudoinverse for Gibbs sampler'
         return
      end if

      do while (stored < n)
         do i = 1, p
            if (precision(i, i) <= 0.0_dp) then
               out%status = 6
               out%message = 'nonpositive conditional precision in Gibbs sampler'
               return
            end if
            conditional_var = 1.0_dp / precision(i, i)
            conditional_mean = mu(i)
            do k = 1, p
               if (k /= i) conditional_mean = conditional_mean - &
                  precision(i, k) * (draw(k) - mu(k)) / precision(i, i)
            end do

            lower = -huge(1.0_dp)
            upper = huge(1.0_dp)
            impossible = .false.
            do j = 1, m
               ai = a(j, i)
               rhs = b(j)
               do k = 1, p
                  if (k /= i) rhs = rhs - a(j, k) * draw(k)
               end do
               if (ai > 0.0_dp) then
                  upper = min(upper, rhs / ai)
               else if (ai < 0.0_dp) then
                  lower = max(lower, rhs / ai)
               else if (rhs < -100.0_dp * epsilon(1.0_dp) * max(1.0_dp, abs(b(j)))) then
                  impossible = .true.
               end if
            end do
            if (impossible .or. lower > upper) then
               out%status = 7
               out%message = 'empty conditional interval in Gibbs sampler'
               return
            end if
            draw(i) = truncated_normal(rng, conditional_mean, sqrt(conditional_var), lower, upper)
         end do
         stored = stored + 1
         out%x(stored, :) = draw
         out%ngibbs = out%ngibbs + 1
      end do
   end subroutine rtmvnorm

   subroutine draw_mvn(mu, factor, rng, z, x)
      real(dp), intent(in) :: mu(:), factor(:, :)
      type(rng_state), intent(inout) :: rng
      real(dp), intent(out) :: z(:), x(:)
      integer :: j
      do j = 1, size(z)
         z(j) = rng_normal(rng)
      end do
      x = mu + matmul(factor, z)
   end subroutine draw_mvn

   pure logical function feasible(x, a, b) result(ok)
      real(dp), intent(in) :: x(:), a(:, :), b(:)
      integer :: j
      real(dp) :: lhs, tol
      ok = .true.
      do j = 1, size(b)
         lhs = dot_product(a(j, :), x)
         tol = 100.0_dp * epsilon(1.0_dp) * max(1.0_dp, abs(lhs), abs(b(j)))
         if (lhs > b(j) + tol) then
            ok = .false.
            return
         end if
      end do
   end function feasible

end module ceoptim_sampling
