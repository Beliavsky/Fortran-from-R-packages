! SPDX-License-Identifier: GPL-2.0-or-later
module flexmix_numeric
   use flexmix_kinds, only : dp
   implicit none
   private
   real(dp), parameter :: pi = acos(-1.0_dp)
   real(dp), parameter :: tiny_prob = 1.0e-300_dp
   public :: normalize_log_probabilities
   public :: weighted_mean_vector
   public :: weighted_covariance
   public :: weighted_least_squares
   public :: solve_linear_system
   public :: normal_logpdf
   public :: poisson_logpmf
   public :: binomial_logpmf
   public :: mvnormal_logpdf_rows
   public :: inverse_logdet_spd
   public :: symmetric_eigen_jacobi
   public :: digamma_approx
   public :: trigamma_approx
   public :: clip_probability
   public :: weighted_variance_unbiased

contains

   pure elemental function clip_probability(x) result(y)
      real(dp), intent(in) :: x !! Probability-like value to constrain away from exact zero and one.
      real(dp) :: y
      real(dp), parameter :: eps = 1.0e-12_dp
      y = min(1.0_dp - eps, max(eps, x))
   end function clip_probability

   pure elemental function normal_logpdf(x, mean_value, sd_value) result(value)
      real(dp), intent(in) :: x !! Observation whose Gaussian log density is required.
      real(dp), intent(in) :: mean_value !! Gaussian mean in the same units as `x`.
      real(dp), intent(in) :: sd_value !! Positive Gaussian standard deviation.
      real(dp) :: value
      real(dp) :: sd_safe
      sd_safe = max(sd_value, sqrt(tiny(1.0_dp)))
      value = -0.5_dp * log(2.0_dp * pi) - log(sd_safe) - 0.5_dp * ((x - mean_value) / sd_safe)**2
   end function normal_logpdf

   pure elemental function poisson_logpmf(y, lambda_value) result(value)
      real(dp), intent(in) :: y !! Nonnegative count, represented in working precision.
      real(dp), intent(in) :: lambda_value !! Positive Poisson mean parameter.
      real(dp) :: value
      real(dp) :: lam
      lam = max(lambda_value, tiny_prob)
      if (y < 0.0_dp) then
         value = -huge(1.0_dp)
      else
         value = y * log(lam) - lam - log_gamma(y + 1.0_dp)
      end if
   end function poisson_logpmf

   pure elemental function binomial_logpmf(success, trials, probability) result(value)
      real(dp), intent(in) :: success !! Number of successful Bernoulli trials.
      real(dp), intent(in) :: trials !! Total number of Bernoulli trials.
      real(dp), intent(in) :: probability !! Success probability in the closed unit interval.
      real(dp) :: value
      real(dp) :: p
      if (success < 0.0_dp .or. trials < success) then
         value = -huge(1.0_dp)
         return
      end if
      p = clip_probability(probability)
      value = log_gamma(trials + 1.0_dp) - log_gamma(success + 1.0_dp) - &
              log_gamma(trials - success + 1.0_dp) + success * log(p) + &
              (trials - success) * log(1.0_dp - p)
   end function binomial_logpmf

   pure subroutine normalize_log_probabilities(logp, posterior, log_row_sum)
      real(dp), intent(in) :: logp(:,:) !! Row-wise component log weights before normalization, shape `(n, k)`.
      real(dp), intent(out) :: posterior(:,:) !! Row-wise normalized component probabilities, shape `(n, k)`.
      real(dp), intent(out) :: log_row_sum(:) !! Stable log of each row sum of exponentiated log weights, size `n`.
      integer :: i, k
      real(dp) :: m, s
      do i = 1, size(logp, 1)
         m = maxval(logp(i,:))
         if (m <= -0.5_dp * huge(1.0_dp)) then
            posterior(i,:) = 1.0_dp / real(size(logp, 2), dp)
            log_row_sum(i) = -huge(1.0_dp)
         else
            s = 0.0_dp
            do k = 1, size(logp, 2)
               posterior(i,k) = exp(logp(i,k) - m)
               s = s + posterior(i,k)
            end do
            if (s <= 0.0_dp) then
               posterior(i,:) = 1.0_dp / real(size(logp, 2), dp)
               log_row_sum(i) = m
            else
               posterior(i,:) = posterior(i,:) / s
               log_row_sum(i) = m + log(s)
            end if
         end if
      end do
   end subroutine normalize_log_probabilities

   pure subroutine weighted_mean_vector(x, weights, center, ok)
      real(dp), intent(in) :: x(:,:) !! Observations by variables, shape `(n, d)`.
      real(dp), intent(in) :: weights(:) !! Nonnegative case weights, size `n`.
      real(dp), intent(out) :: center(:) !! Weighted variable means, size `d`.
      logical, intent(out) :: ok !! True when a positive finite total weight permits estimation.
      real(dp) :: sw
      integer :: j
      sw = sum(weights)
      ok = sw > tiny(1.0_dp)
      if (.not. ok) then
         center = 0.0_dp
         return
      end if
      do j = 1, size(x, 2)
         center(j) = dot_product(weights, x(:,j)) / sw
      end do
   end subroutine weighted_mean_vector

   pure subroutine weighted_covariance(x, weights, center, covariance, ok)
      real(dp), intent(in) :: x(:,:) !! Observations by variables, shape `(n, d)`.
      real(dp), intent(in) :: weights(:) !! Nonnegative case weights, size `n`, interpreted as `cov.wt` frequency weights.
      real(dp), intent(out) :: center(:) !! Weighted means, size `d`.
      real(dp), intent(out) :: covariance(:,:) !! Unbiased `cov.wt`-style covariance matrix, shape `(d, d)`.
      logical, intent(out) :: ok !! True when the effective weight denominator is positive.
      real(dp) :: sw, denom, wi
      real(dp), allocatable :: wn(:)
      integer :: i, j, l
      call weighted_mean_vector(x, weights, center, ok)
      covariance = 0.0_dp
      if (.not. ok) return
      sw = sum(weights)
      allocate(wn(size(weights)))
      wn = weights / sw
      denom = 1.0_dp - sum(wn * wn)
      if (denom <= 100.0_dp * epsilon(1.0_dp)) then
         ok = .false.
         return
      end if
      do i = 1, size(x, 1)
         wi = wn(i) / denom
         do j = 1, size(x, 2)
            do l = 1, j
               covariance(j,l) = covariance(j,l) + wi * (x(i,j) - center(j)) * (x(i,l) - center(l))
            end do
         end do
      end do
      do j = 1, size(x, 2)
         do l = 1, j - 1
            covariance(l,j) = covariance(j,l)
         end do
      end do
   end subroutine weighted_covariance

   pure subroutine weighted_variance_unbiased(x, weights, mean_value, variance, ok)
      real(dp), intent(in) :: x(:) !! Scalar observations whose weighted variance is required.
      real(dp), intent(in) :: weights(:) !! Nonnegative case weights aligned with `x`.
      real(dp), intent(out) :: mean_value !! Weighted arithmetic mean.
      real(dp), intent(out) :: variance !! Unbiased `cov.wt`-style weighted variance.
      logical, intent(out) :: ok !! True when a positive effective denominator permits estimation.
      real(dp) :: sw, denom
      real(dp), allocatable :: wn(:)
      sw = sum(weights)
      ok = sw > tiny(1.0_dp)
      if (.not. ok) then
         mean_value = 0.0_dp
         variance = 0.0_dp
         return
      end if
      mean_value = dot_product(weights, x) / sw
      allocate(wn(size(weights)))
      wn = weights / sw
      denom = 1.0_dp - sum(wn * wn)
      if (denom <= 100.0_dp * epsilon(1.0_dp)) then
         variance = 0.0_dp
         ok = .false.
         return
      end if
      variance = sum(wn * (x - mean_value)**2) / denom
   end subroutine weighted_variance_unbiased

   pure subroutine symmetric_eigen_jacobi(a, values, vectors, info)
      real(dp), intent(in) :: a(:,:) !! Real symmetric matrix, shape `(n,n)`, whose eigenpairs are required.
      real(dp), intent(out) :: values(:) !! Eigenvalues sorted from largest to smallest, size `n`.
      real(dp), intent(out) :: vectors(:,:) !! Corresponding orthonormal eigenvectors stored by column, shape `(n,n)`.
      integer, intent(out) :: info !! Zero on convergence; nonzero for incompatible shapes or failure to diagonalize.
      real(dp), allocatable :: work(:,:)
      real(dp) :: app, aqq, apq, akp, akq, vip, viq, c, s, theta, offmax, scale, tmp
      integer :: i, iter, j, k, n, p, q, max_iter

      n = size(a,1)
      info = 0
      values = 0.0_dp
      vectors = 0.0_dp
      if (size(a,2) /= n .or. size(values) /= n .or. any(shape(vectors) /= [n,n])) then
         info = -1
         return
      end if
      if (n == 0) return
      allocate(work(n,n))
      work = 0.5_dp * (a + transpose(a))
      do i = 1, n
         vectors(i,i) = 1.0_dp
      end do
      max_iter = max(20, 80*n*n)
      do iter = 1, max_iter
         offmax = 0.0_dp
         p = 1
         q = 1
         do j = 2, n
            do i = 1, j - 1
               if (abs(work(i,j)) > offmax) then
                  offmax = abs(work(i,j))
                  p = i
                  q = j
               end if
            end do
         end do
         scale = max(1.0_dp, maxval(abs([(work(i,i), i=1,n)])))
         if (offmax <= 100.0_dp * epsilon(1.0_dp) * scale) exit
         app = work(p,p)
         aqq = work(q,q)
         apq = work(p,q)
         theta = 0.5_dp * atan2(2.0_dp * apq, app - aqq)
         c = cos(theta)
         s = sin(theta)
         do k = 1, n
            if (k == p .or. k == q) cycle
            akp = work(k,p)
            akq = work(k,q)
            work(k,p) = c * akp + s * akq
            work(p,k) = work(k,p)
            work(k,q) = -s * akp + c * akq
            work(q,k) = work(k,q)
         end do
         work(p,p) = c*c*app + 2.0_dp*s*c*apq + s*s*aqq
         work(q,q) = s*s*app - 2.0_dp*s*c*apq + c*c*aqq
         work(p,q) = 0.0_dp
         work(q,p) = 0.0_dp
         do k = 1, n
            vip = vectors(k,p)
            viq = vectors(k,q)
            vectors(k,p) = c * vip + s * viq
            vectors(k,q) = -s * vip + c * viq
         end do
      end do
      if (iter > max_iter) then
         info = 1
         return
      end if
      do i = 1, n
         values(i) = work(i,i)
      end do
      do i = 1, n - 1
         q = i
         do j = i + 1, n
            if (values(j) > values(q)) q = j
         end do
         if (q /= i) then
            tmp = values(i)
            values(i) = values(q)
            values(q) = tmp
            work(:,1) = vectors(:,i)
            vectors(:,i) = vectors(:,q)
            vectors(:,q) = work(:,1)
         end if
      end do
   end subroutine symmetric_eigen_jacobi

   pure subroutine solve_linear_system(a, b, x, info)
      real(dp), intent(in) :: a(:,:) !! Square coefficient matrix, shape `(n, n)`.
      real(dp), intent(in) :: b(:) !! Right-hand-side vector, size `n`.
      real(dp), intent(out) :: x(:) !! Solution vector, size `n`.
      integer, intent(out) :: info !! Zero on success; nonzero when dimensions or pivots are singular.
      real(dp), allocatable :: work(:,:), rhs(:), tmp_row(:)
      real(dp) :: factor, pivot_abs, rhs_tmp
      integer :: i, k, n, pivot
      n = size(a, 1)
      info = 0
      x = 0.0_dp
      if (size(a, 2) /= n .or. size(b) /= n .or. size(x) /= n) then
         info = -1
         return
      end if
      if (n == 0) return
      allocate(work(n,n), rhs(n), tmp_row(n))
      work = a
      rhs = b
      do k = 1, n - 1
         pivot = k
         pivot_abs = abs(work(k,k))
         do i = k + 1, n
            if (abs(work(i,k)) > pivot_abs) then
               pivot_abs = abs(work(i,k))
               pivot = i
            end if
         end do
         if (pivot_abs <= 100.0_dp * epsilon(1.0_dp) * max(1.0_dp, maxval(abs(work)))) then
            info = k
            return
         end if
         if (pivot /= k) then
            tmp_row = work(k,:)
            work(k,:) = work(pivot,:)
            work(pivot,:) = tmp_row
            rhs_tmp = rhs(k)
            rhs(k) = rhs(pivot)
            rhs(pivot) = rhs_tmp
         end if
         do i = k + 1, n
            factor = work(i,k) / work(k,k)
            work(i,k) = 0.0_dp
            work(i,k+1:n) = work(i,k+1:n) - factor * work(k,k+1:n)
            rhs(i) = rhs(i) - factor * rhs(k)
         end do
      end do
      if (abs(work(n,n)) <= 100.0_dp * epsilon(1.0_dp) * max(1.0_dp, maxval(abs(work)))) then
         info = n
         return
      end if
      do i = n, 1, -1
         if (i < n) then
            x(i) = (rhs(i) - dot_product(work(i,i+1:n), x(i+1:n))) / work(i,i)
         else
            x(i) = rhs(i) / work(i,i)
         end if
      end do
   end subroutine solve_linear_system

   pure subroutine weighted_least_squares(x, y, weights, beta, rank_value, residuals, info)
      real(dp), intent(in) :: x(:,:) !! Regression design matrix, shape `(n, p)`.
      real(dp), intent(in) :: y(:) !! Response vector, size `n`.
      real(dp), intent(in) :: weights(:) !! Nonnegative observation weights, size `n`.
      real(dp), intent(out) :: beta(:) !! Weighted least-squares coefficients, size `p`.
      integer, intent(out) :: rank_value !! Numerical rank estimated from the weighted QR factorization.
      real(dp), intent(out) :: residuals(:) !! Unweighted response residuals `y - X beta`, size `n`.
      integer, intent(out) :: info !! Zero on success; nonzero for dimension/rank failure.
      real(dp), allocatable :: a(:,:), q(:,:), r(:,:), v(:), qty(:), yw(:)
      real(dp) :: normv, threshold, scale
      integer :: i, j, n, p
      n = size(x, 1)
      p = size(x, 2)
      beta = 0.0_dp
      residuals = y
      rank_value = 0
      info = 0
      if (size(y) /= n .or. size(weights) /= n .or. size(beta) /= p .or. size(residuals) /= n) then
         info = -1
         return
      end if
      if (n == 0 .or. p == 0) then
         info = -2
         return
      end if
      allocate(a(n,p), q(n,p), r(p,p), v(n), qty(p), yw(n))
      a = x
      do i = 1, n
         scale = sqrt(max(0.0_dp, weights(i)))
         a(i,:) = a(i,:) * scale
         yw(i) = y(i) * scale
      end do
      q = 0.0_dp
      r = 0.0_dp
      threshold = sqrt(epsilon(1.0_dp)) * max(1.0_dp, maxval(abs(a)))
      do j = 1, p
         v = a(:,j)
         do i = 1, j - 1
            r(i,j) = dot_product(q(:,i), v)
            v = v - r(i,j) * q(:,i)
         end do
         do i = 1, j - 1
            scale = dot_product(q(:,i), v)
            r(i,j) = r(i,j) + scale
            v = v - scale * q(:,i)
         end do
         normv = sqrt(sum(v * v))
         r(j,j) = normv
         if (normv > threshold) then
            q(:,j) = v / normv
            rank_value = rank_value + 1
         else
            q(:,j) = 0.0_dp
         end if
      end do
      if (rank_value < p) then
         info = 1
         return
      end if
      do j = 1, p
         qty(j) = dot_product(q(:,j), yw)
      end do
      do i = p, 1, -1
         if (i < p) then
            beta(i) = (qty(i) - dot_product(r(i,i+1:p), beta(i+1:p))) / r(i,i)
         else
            beta(i) = qty(i) / r(i,i)
         end if
      end do
      residuals = y - matmul(x, beta)
   end subroutine weighted_least_squares

   pure subroutine cholesky_with_jitter(a, l, info)
      real(dp), intent(in) :: a(:,:) !! Symmetric matrix expected to be positive definite.
      real(dp), intent(out) :: l(:,:) !! Lower-triangular Cholesky factor when successful.
      integer, intent(out) :: info !! Zero on success; positive if regularization failed.
      real(dp), allocatable :: work(:,:)
      real(dp) :: jitter, s, base
      integer :: attempt, i, j, k, n
      n = size(a, 1)
      l = 0.0_dp
      info = 0
      if (size(a, 2) /= n .or. size(l, 1) /= n .or. size(l, 2) /= n) then
         info = -1
         return
      end if
      allocate(work(n,n))
      base = max(1.0_dp, maxval(abs(a)))
      jitter = 0.0_dp
      do attempt = 1, 8
         work = a
         if (jitter > 0.0_dp) then
            do i = 1, n
               work(i,i) = work(i,i) + jitter
            end do
         end if
         l = 0.0_dp
         info = 0
         do i = 1, n
            do j = 1, i
               s = work(i,j)
               do k = 1, j - 1
                  s = s - l(i,k) * l(j,k)
               end do
               if (i == j) then
                  if (s <= 0.0_dp) then
                     info = i
                     exit
                  end if
                  l(i,j) = sqrt(s)
               else
                  l(i,j) = s / l(j,j)
               end if
            end do
            if (info /= 0) exit
         end do
         if (info == 0) return
         if (attempt == 1) then
            jitter = 1.0e-10_dp * base
         else
            jitter = 10.0_dp * jitter
         end if
      end do
   end subroutine cholesky_with_jitter

   pure subroutine solve_lower(l, b, x)
      real(dp), intent(in) :: l(:,:) !! Nonsingular lower-triangular coefficient matrix.
      real(dp), intent(in) :: b(:) !! Right-hand side aligned with `l`.
      real(dp), intent(out) :: x(:) !! Forward-substitution solution.
      integer :: i
      do i = 1, size(b)
         if (i == 1) then
            x(i) = b(i) / l(i,i)
         else
            x(i) = (b(i) - dot_product(l(i,1:i-1), x(1:i-1))) / l(i,i)
         end if
      end do
   end subroutine solve_lower

   pure subroutine solve_upper_from_lower(l, b, x)
      real(dp), intent(in) :: l(:,:) !! Lower Cholesky factor whose transpose is solved against.
      real(dp), intent(in) :: b(:) !! Right-hand side for the transposed triangular system.
      real(dp), intent(out) :: x(:) !! Back-substitution solution.
      integer :: i, n
      n = size(b)
      do i = n, 1, -1
         if (i == n) then
            x(i) = b(i) / l(i,i)
         else
            x(i) = (b(i) - dot_product(l(i+1:n,i), x(i+1:n))) / l(i,i)
         end if
      end do
   end subroutine solve_upper_from_lower

   pure subroutine mvnormal_logpdf_rows(x, center, covariance, log_density, info)
      real(dp), intent(in) :: x(:,:) !! Observations by variables, shape `(n, d)`.
      real(dp), intent(in) :: center(:) !! Multivariate-normal center, size `d`.
      real(dp), intent(in) :: covariance(:,:) !! Positive-definite covariance, shape `(d, d)`.
      real(dp), intent(out) :: log_density(:) !! Log density for each observation, size `n`.
      integer, intent(out) :: info !! Zero on successful Cholesky factorization.
      real(dp), allocatable :: l(:,:), delta(:), z(:)
      real(dp) :: logdet, quad
      integer :: i, d
      d = size(center)
      if (size(x, 2) /= d .or. size(covariance, 1) /= d .or. size(covariance, 2) /= d .or. &
          size(log_density) /= size(x, 1)) then
         info = -1
         log_density = -huge(1.0_dp)
         return
      end if
      allocate(l(d,d), delta(d), z(d))
      call cholesky_with_jitter(covariance, l, info)
      if (info /= 0) then
         log_density = -huge(1.0_dp)
         return
      end if
      logdet = 2.0_dp * sum(log([(l(i,i), i = 1, d)]))
      do i = 1, size(x, 1)
         delta = x(i,:) - center
         call solve_lower(l, delta, z)
         quad = sum(z * z)
         log_density(i) = -0.5_dp * (real(d, dp) * log(2.0_dp * pi) + logdet + quad)
      end do
   end subroutine mvnormal_logpdf_rows

   pure subroutine inverse_logdet_spd(a, inverse, logdet, info)
      real(dp), intent(in) :: a(:,:) !! Positive-definite matrix to invert.
      real(dp), intent(out) :: inverse(:,:) !! Inverse matrix with the same square shape as `a`.
      real(dp), intent(out) :: logdet !! Natural logarithm of the positive determinant.
      integer, intent(out) :: info !! Zero on success, nonzero if factorization fails.
      real(dp), allocatable :: l(:,:), e(:), z(:), x(:)
      integer :: i, n
      n = size(a, 1)
      inverse = 0.0_dp
      logdet = 0.0_dp
      if (size(a, 2) /= n .or. size(inverse, 1) /= n .or. size(inverse, 2) /= n) then
         info = -1
         return
      end if
      allocate(l(n,n), e(n), z(n), x(n))
      call cholesky_with_jitter(a, l, info)
      if (info /= 0) return
      logdet = 2.0_dp * sum(log([(l(i,i), i = 1, n)]))
      do i = 1, n
         e = 0.0_dp
         e(i) = 1.0_dp
         call solve_lower(l, e, z)
         call solve_upper_from_lower(l, z, x)
         inverse(:,i) = x
      end do
   end subroutine inverse_logdet_spd

   pure elemental function digamma_approx(x) result(value)
      real(dp), intent(in) :: x !! Positive argument for the digamma approximation.
      real(dp) :: value
      real(dp) :: y, inv, inv2
      y = x
      if (y <= 0.0_dp) then
         value = huge(1.0_dp)
         return
      end if
      value = 0.0_dp
      do while (y < 8.0_dp)
         value = value - 1.0_dp / y
         y = y + 1.0_dp
      end do
      inv = 1.0_dp / y
      inv2 = inv * inv
      value = value + log(y) - 0.5_dp * inv - inv2 * (1.0_dp / 12.0_dp - &
              inv2 * (1.0_dp / 120.0_dp - inv2 * (1.0_dp / 252.0_dp)))
   end function digamma_approx

   pure elemental function trigamma_approx(x) result(value)
      real(dp), intent(in) :: x !! Positive argument for the trigamma approximation.
      real(dp) :: value
      real(dp) :: y, inv, inv2
      y = x
      if (y <= 0.0_dp) then
         value = huge(1.0_dp)
         return
      end if
      value = 0.0_dp
      do while (y < 8.0_dp)
         value = value + 1.0_dp / (y * y)
         y = y + 1.0_dp
      end do
      inv = 1.0_dp / y
      inv2 = inv * inv
      value = value + inv + 0.5_dp * inv2 + inv * inv2 / 6.0_dp - &
              inv * inv2 * inv2 / 30.0_dp + inv * inv2 * inv2 * inv2 / 42.0_dp
   end function trigamma_approx

end module flexmix_numeric
