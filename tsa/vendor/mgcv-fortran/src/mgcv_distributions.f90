! SPDX-License-Identifier: GPL-2.0-or-later
module mgcv_distributions
   use mgcv_kinds, only : dp, pi_dp
   use mgcv_linalg, only : cholesky_upper, spd_solve
   use mgcv_utils, only : log_expm1, log1p_stable
   implicit none
   private
   public :: normal_cdf, normal_pdf, dpnorm
   public :: rmvn, dmvn_log, rmvt, dmvt_log
   public :: random_normal, random_gamma, random_poisson
   public :: tweedie_deviance, tweedie_variance, rtweedie
   public :: weighted_chisq_cdf

contains

   elemental function normal_pdf(x) result(y)
      real(dp), intent(in) :: x
      real(dp) :: y
      y = exp(-0.5_dp * x * x) / sqrt(2.0_dp * pi_dp)
   end function normal_pdf

   elemental function normal_cdf(x) result(y)
      real(dp), intent(in) :: x
      real(dp) :: y
      y = 0.5_dp * erfc(-x / sqrt(2.0_dp))
   end function normal_cdf

   elemental function dpnorm(x0, x1, log_p) result(value)
      real(dp), intent(in) :: x0, x1
      logical, intent(in), optional :: log_p
      real(dp) :: value, a, b, p0, p1, lp
      logical :: want_log
      want_log = .false.; if (present(log_p)) want_log = log_p
      if (x0 >= x1) then
         value = merge(-huge(1.0_dp), 0.0_dp, want_log); return
      end if
      a = x0; b = x1
      if (a > 0.0_dp .and. b > 0.0_dp) then
         a = -x1; b = -x0
      end if
      p0 = normal_cdf(a); p1 = normal_cdf(b)
      if (p0 <= 0.0_dp) then
         lp = log(max(tiny(1.0_dp), p1))
      else
         lp = log(p0) + log_expm1(log(p1) - log(p0))
      end if
      value = merge(lp, exp(lp), want_log)
   end function dpnorm

   subroutine set_seed(seed)
      integer, intent(in) :: seed
      integer :: n, i
      integer, allocatable :: values(:)
      call random_seed(size=n)
      allocate(values(n))
      do i = 1, n
         values(i) = modulo(seed + 104729 * i, huge(1) - 1)
         if (values(i) <= 0) values(i) = i
      end do
      call random_seed(put=values)
   end subroutine set_seed

   real(dp) function random_normal() result(z)
      real(dp) :: u1, u2
      call random_number(u1); call random_number(u2)
      u1 = max(u1, tiny(1.0_dp))
      z = sqrt(-2.0_dp * log(u1)) * cos(2.0_dp * pi_dp * u2)
   end function random_normal

   recursive real(dp) function random_gamma(shape, scale) result(x)
      real(dp), intent(in) :: shape, scale
      real(dp) :: d, c, z, v, u
      if (shape <= 0.0_dp .or. scale <= 0.0_dp) then
         x = 0.0_dp; return
      end if
      if (shape < 1.0_dp) then
         call random_number(u)
         x = random_gamma(shape + 1.0_dp, scale) * max(u, tiny(1.0_dp))**(1.0_dp / shape)
         return
      end if
      d = shape - 1.0_dp / 3.0_dp
      c = 1.0_dp / sqrt(9.0_dp * d)
      do
         z = random_normal()
         v = 1.0_dp + c * z
         if (v <= 0.0_dp) cycle
         v = v**3
         call random_number(u)
         if (u < 1.0_dp - 0.0331_dp * z**4) exit
         if (log(max(u, tiny(1.0_dp))) < 0.5_dp * z * z + d * (1.0_dp - v + log(v))) exit
      end do
      x = scale * d * v
   end function random_gamma

   integer function random_poisson(lambda) result(k)
      real(dp), intent(in) :: lambda
      real(dp) :: l, p, u, z
      integer :: n
      if (lambda <= 0.0_dp) then
         k = 0; return
      end if
      if (lambda < 30.0_dp) then
         l = exp(-lambda); p = 1.0_dp; n = 0
         do
            n = n + 1; call random_number(u); p = p * u
            if (p <= l) exit
         end do
         k = n - 1
      else
         do
            z = lambda + sqrt(lambda) * random_normal()
            if (z >= 0.0_dp) exit
         end do
         k = nint(z)
      end if
   end function random_poisson

   subroutine rmvn(n, mu, covariance, sample, status, seed)
      integer, intent(in) :: n
      real(dp), intent(in) :: mu(:), covariance(:, :)
      real(dp), allocatable, intent(out) :: sample(:, :)
      integer, intent(out) :: status
      integer, intent(in), optional :: seed
      real(dp), allocatable :: r(:, :), z(:)
      integer :: i, j, p
      if (present(seed)) call set_seed(seed)
      p = size(mu)
      if (size(covariance, 1) /= p .or. size(covariance, 2) /= p .or. n < 1) then
         allocate(sample(0, 0)); status = 1; return
      end if
      call cholesky_upper(covariance, r, status, 1.0e-12_dp)
      if (status /= 0) then
         allocate(sample(0, 0)); return
      end if
      allocate(sample(n, p), z(p))
      do i = 1, n
         do j = 1, p
            z(j) = random_normal()
         end do
         sample(i, :) = mu + matmul(z, r)
      end do
   end subroutine rmvn

   function dmvn_log(x, mu, covariance, status) result(value)
      real(dp), intent(in) :: x(:), mu(:), covariance(:, :)
      integer, intent(out) :: status
      real(dp) :: value
      real(dp), allocatable :: r(:, :), rhs(:, :), z(:, :)
      integer :: i, p
      p = size(mu)
      if (size(x) /= p .or. size(covariance, 1) /= p .or. size(covariance, 2) /= p) then
         status = 1; value = -huge(1.0_dp); return
      end if
      call cholesky_upper(covariance, r, status, 1.0e-12_dp)
      if (status /= 0) then
         value = -huge(1.0_dp); return
      end if
      allocate(rhs(p, 1)); rhs(:, 1) = x - mu
      call spd_solve(covariance, rhs, z, status, 1.0e-12_dp)
      if (status /= 0) then
         value = -huge(1.0_dp); return
      end if
      value = -0.5_dp * real(p, dp) * log(2.0_dp * pi_dp)
      do i = 1, p
         value = value - log(r(i, i))
      end do
      value = value - 0.5_dp * dot_product(x - mu, z(:, 1))
   end function dmvn_log

   subroutine rmvt(n, mu, covariance, df, sample, status, seed)
      integer, intent(in) :: n
      real(dp), intent(in) :: mu(:), covariance(:, :), df
      real(dp), allocatable, intent(out) :: sample(:, :)
      integer, intent(out) :: status
      integer, intent(in), optional :: seed
      real(dp), allocatable :: z(:, :)
      real(dp) :: chi
      integer :: i
      if (df <= 0.0_dp) then
         allocate(sample(0, 0)); status = 1; return
      end if
      call rmvn(n, 0.0_dp * mu, covariance, z, status, seed)
      if (status /= 0) then
         allocate(sample(0, 0)); return
      end if
      allocate(sample(n, size(mu)))
      do i = 1, n
         chi = random_gamma(df / 2.0_dp, 2.0_dp)
         sample(i, :) = mu + z(i, :) * sqrt(df / max(chi, tiny(1.0_dp)))
      end do
   end subroutine rmvt

   function dmvt_log(x, mu, covariance, df, status) result(value)
      real(dp), intent(in) :: x(:), mu(:), covariance(:, :), df
      integer, intent(out) :: status
      real(dp) :: value, quad
      real(dp), allocatable :: r(:, :), rhs(:, :), z(:, :)
      integer :: i, p
      p = size(mu)
      if (df <= 0.0_dp .or. size(x) /= p .or. size(covariance, 1) /= p .or. &
          size(covariance, 2) /= p) then
         status = 1; value = -huge(1.0_dp); return
      end if
      call cholesky_upper(covariance, r, status, 1.0e-12_dp)
      if (status /= 0) then
         value = -huge(1.0_dp); return
      end if
      allocate(rhs(p, 1)); rhs(:, 1) = x - mu
      call spd_solve(covariance, rhs, z, status, 1.0e-12_dp)
      if (status /= 0) then
         value = -huge(1.0_dp); return
      end if
      quad = dot_product(x - mu, z(:, 1))
      value = log_gamma((df + real(p, dp)) / 2.0_dp) - log_gamma(df / 2.0_dp)
      value = value - 0.5_dp * real(p, dp) * log(df * pi_dp)
      do i = 1, p
         value = value - log(r(i, i))
      end do
      value = value - 0.5_dp * (df + real(p, dp)) * log1p_stable(quad / df)
   end function dmvt_log

   elemental function tweedie_variance(mu, p) result(value)
      real(dp), intent(in) :: mu, p
      real(dp) :: value
      value = max(mu, tiny(1.0_dp))**p
   end function tweedie_variance

   elemental function tweedie_deviance(y, mu, p) result(value)
      real(dp), intent(in) :: y, mu, p
      real(dp) :: value, yy, mm
      yy = max(y, 0.0_dp); mm = max(mu, tiny(1.0_dp))
      if (abs(p - 1.0_dp) < 1.0e-8_dp) then
         if (yy <= tiny(1.0_dp)) then
            value = 2.0_dp * mm
         else
            value = 2.0_dp * (yy * log(yy / mm) - (yy - mm))
         end if
      else if (abs(p - 2.0_dp) < 1.0e-8_dp) then
         if (yy <= 0.0_dp) then
            value = huge(1.0_dp)
         else
            value = 2.0_dp * (yy / mm - log(yy / mm) - 1.0_dp)
         end if
      else
         value = 2.0_dp * (yy**(2.0_dp - p) / ((1.0_dp - p) * (2.0_dp - p)) - &
                 yy * mm**(1.0_dp - p) / (1.0_dp - p) + mm**(2.0_dp - p) / (2.0_dp - p))
      end if
   end function tweedie_deviance

   subroutine rtweedie(n, mu, phi, p, x, status, seed)
      integer, intent(in) :: n
      real(dp), intent(in) :: mu, phi, p
      real(dp), allocatable, intent(out) :: x(:)
      integer, intent(out) :: status
      integer, intent(in), optional :: seed
      real(dp) :: lambda, alpha, gamma_scale
      integer :: i, count
      if (present(seed)) call set_seed(seed)
      if (n < 1 .or. mu < 0.0_dp .or. phi <= 0.0_dp .or. p <= 1.0_dp .or. p >= 2.0_dp) then
         allocate(x(0)); status = 1; return
      end if
      lambda = mu**(2.0_dp - p) / (phi * (2.0_dp - p))
      alpha = (2.0_dp - p) / (p - 1.0_dp)
      gamma_scale = phi * (p - 1.0_dp) * mu**(p - 1.0_dp)
      allocate(x(n)); status = 0
      do i = 1, n
         count = random_poisson(lambda)
         if (count == 0) then
            x(i) = 0.0_dp
         else
            x(i) = random_gamma(real(count, dp) * alpha, gamma_scale)
         end if
      end do
   end subroutine rtweedie

   function weighted_chisq_cdf(q, lambda, df, noncentral, status, n_grid) result(prob)
      real(dp), intent(in) :: q, lambda(:)
      real(dp), intent(in), optional :: df(:), noncentral(:)
      integer, intent(out) :: status
      integer, intent(in), optional :: n_grid
      real(dp) :: prob, h, t, phase, log_amp, integrand, sumv, dfi, nci, lam
      integer :: i, j, n
      status = 0
      if (size(lambda) == 0 .or. any(lambda < 0.0_dp)) then
         status = 1; prob = 0.0_dp; return
      end if
      if (q <= 0.0_dp) then
         prob = 0.0_dp; return
      end if
      n = 4000; if (present(n_grid)) n = max(200, n_grid)
      if (mod(n, 2) == 1) n = n + 1
      h = 40.0_dp / real(n, dp)
      sumv = 0.0_dp
      do i = 1, n
         t = (real(i, dp) - 0.5_dp) * h
         phase = -t * q
         log_amp = 0.0_dp
         do j = 1, size(lambda)
            lam = lambda(j)
            dfi = 1.0_dp; if (present(df)) dfi = df(j)
            nci = 0.0_dp; if (present(noncentral)) nci = noncentral(j)
            phase = phase + 0.5_dp * dfi * atan(2.0_dp * t * lam)
            phase = phase + nci * t * lam / (1.0_dp + 4.0_dp * t * t * lam * lam)
            log_amp = log_amp - 0.25_dp * dfi * log1p_stable(4.0_dp * t * t * lam * lam)
            log_amp = log_amp - 2.0_dp * nci * t * t * lam * lam / &
                      (1.0_dp + 4.0_dp * t * t * lam * lam)
         end do
         integrand = exp(log_amp) * sin(phase) / t
         sumv = sumv + integrand
      end do
      prob = max(0.0_dp, min(1.0_dp, 0.5_dp - h * sumv / pi_dp))
   end function weighted_chisq_cdf

end module mgcv_distributions
