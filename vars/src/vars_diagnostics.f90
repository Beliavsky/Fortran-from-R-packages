! SPDX-License-Identifier: GPL-2.0-or-later
! Computational translation of R package vars 1.6-1; see NOTICE.md and UPSTREAM.md.
module vars_diagnostics
   use r_kinds, only : dp
   use r_distributions, only : r_pchisq, r_pf
   use r_linalg, only : inverse_matrix, least_squares_svd, cholesky_factor
   use vars_types
   use vars_utils, only : covariance_matrix, determinant_logabs, trace_matrix
   use vars_utils, only : vech_lower, residual_standardize, center_columns
   implicit none
   private

   public :: arch_test_univariate, arch_test_multivariate
   public :: jarque_bera_univariate, jarque_bera_multivariate
   public :: portmanteau_tests, bg_serial_tests

contains

   subroutine arch_test_univariate(x, lags, result, info)
      real(dp), intent(in) :: x(:)
      integer, intent(in) :: lags
      type(vars_test_result), intent(out) :: result
      integer, intent(out) :: info
      real(dp), allocatable :: z(:), y(:), design(:, :), beta(:), resid(:)
      real(dp) :: mean_y, sst, sse, r2, mu, sd
      integer :: n, m, i, j, rank

      n = size(x)
      if (lags < 1 .or. n <= lags + 1) then
         info = vars_invalid_argument
         return
      end if
      mu = sum(x) / real(n, dp)
      sd = sqrt(sum((x - mu) ** 2) / real(n - 1, dp))
      if (sd <= 0.0_dp) then
         info = vars_singular
         return
      end if
      allocate(z(n))
      z = ((x - mu) / sd) ** 2
      m = n - lags
      allocate(y(m), design(m, lags + 1), beta(lags + 1), resid(m))
      y = z(lags + 1:n)
      design(:, 1) = 1.0_dp
      do j = 1, lags
         do i = 1, m
            design(i, j + 1) = z(lags + i - j)
         end do
      end do
      call least_squares_svd(design, y, beta, rank, info)
      if (info /= 0 .or. rank < lags + 1) then
         info = vars_singular
         return
      end if
      resid = y - matmul(design, beta)
      mean_y = sum(y) / real(m, dp)
      sst = sum((y - mean_y) ** 2)
      sse = sum(resid ** 2)
      if (sst > 0.0_dp) then
         r2 = max(0.0_dp, 1.0_dp - sse / sst)
      else
         r2 = 0.0_dp
      end if
      result%statistic = r2 * real(m, dp)
      result%df1 = real(lags, dp)
      result%df2 = 0.0_dp
      result%p_value = r_pchisq(result%statistic, result%df1, lower_tail = .false.)
      info = vars_success
   end subroutine arch_test_univariate

   subroutine arch_test_multivariate(x, lags, result, info)
      real(dp), intent(in) :: x(:, :)
      integer, intent(in) :: lags
      type(vars_test_result), intent(out) :: result
      integer, intent(out) :: info
      real(dp), allocatable :: z(:, :), vech_data(:, :), y(:, :), design(:, :), beta(:, :)
      real(dp), allocatable :: resid0(:, :), resid1(:, :), omega0(:, :), omega1(:, :), inv0(:, :)
      real(dp), allocatable :: temp(:, :)
      real(dp) :: r2m
      integer :: n, k, q, m, i, j, rank

      n = size(x, 1)
      k = size(x, 2)
      q = k * (k + 1) / 2
      if (lags < 1 .or. n <= lags + q) then
         info = vars_invalid_argument
         return
      end if
      allocate(z(n, k))
      call residual_standardize(x, z)
      allocate(vech_data(n, q), temp(k, k))
      do i = 1, n
         temp = 0.0_dp
         do j = 1, k
            temp(:, j) = z(i, :) * z(i, j)
         end do
         call vech_lower(temp, vech_data(i, :))
      end do
      m = n - lags
      allocate(y(m, q), design(m, 1 + lags * q), beta(1 + lags * q, q))
      allocate(resid0(m, q), resid1(m, q), omega0(q, q), omega1(q, q))
      y = vech_data(lags + 1:n, :)
      design(:, 1) = 1.0_dp
      do j = 1, lags
         design(:, 2 + (j - 1) * q:1 + j * q) = vech_data(lags + 1 - j:n - j, :)
      end do
      do j = 1, q
         resid0(:, j) = y(:, j) - sum(y(:, j)) / real(m, dp)
      end do
      call least_squares_svd(design, y, beta, rank, info)
      if (info /= 0 .or. rank < size(design, 2)) then
         info = vars_singular
         return
      end if
      resid1 = y - matmul(design, beta)
      call covariance_matrix(resid0, omega0)
      call covariance_matrix(resid1, omega1)
      call inverse_matrix(omega0, inv0, info)
      if (info /= 0) then
         info = vars_singular
         return
      end if
      r2m = 1.0_dp - 2.0_dp / real(k * (k + 1), dp) * trace_matrix(matmul(omega1, inv0))
      result%statistic = 0.5_dp * real(m * k * (k + 1), dp) * r2m
      result%df1 = real(lags * k * k * (k + 1) * (k + 1), dp) / 4.0_dp
      result%df2 = 0.0_dp
      result%p_value = r_pchisq(result%statistic, result%df1, lower_tail = .false.)
      info = vars_success
   end subroutine arch_test_multivariate

   subroutine jarque_bera_univariate(x, result, info)
      real(dp), intent(in) :: x(:)
      type(vars_test_result), intent(out) :: result
      integer, intent(out) :: info
      real(dp) :: m1, m2, m3, m4, b1, b2
      integer :: n

      n = size(x)
      if (n < 3) then
         info = vars_invalid_argument
         return
      end if
      m1 = sum(x) / real(n, dp)
      m2 = sum((x - m1) ** 2) / real(n, dp)
      if (m2 <= 0.0_dp) then
         info = vars_singular
         return
      end if
      m3 = sum((x - m1) ** 3) / real(n, dp)
      m4 = sum((x - m1) ** 4) / real(n, dp)
      b1 = (m3 / m2 ** 1.5_dp) ** 2
      b2 = m4 / m2 ** 2
      result%statistic = real(n, dp) * b1 / 6.0_dp + real(n, dp) * (b2 - 3.0_dp) ** 2 / 24.0_dp
      result%df1 = 2.0_dp
      result%df2 = 0.0_dp
      result%p_value = r_pchisq(result%statistic, 2.0_dp, lower_tail = .false.)
      info = vars_success
   end subroutine jarque_bera_univariate

   subroutine jarque_bera_multivariate(x, jb, skewness, kurtosis, info)
      real(dp), intent(in) :: x(:, :)
      type(vars_test_result), intent(out) :: jb, skewness, kurtosis
      integer, intent(out) :: info
      real(dp), allocatable :: centered(:, :), sigma(:, :), upper(:, :), inv_upper(:, :), z(:, :)
      real(dp), allocatable :: b1(:), b2(:)
      real(dp) :: s3, s4
      integer :: n, k, j

      n = size(x, 1)
      k = size(x, 2)
      if (n < 3 .or. k < 1) then
         info = vars_invalid_argument
         return
      end if
      allocate(centered(n, k), sigma(k, k), b1(k), b2(k))
      call center_columns(x, centered)
      sigma = matmul(transpose(centered), centered) / real(n, dp)
      call cholesky_factor(sigma, upper, info, upper = .true.)
      if (info /= 0) then
         info = vars_singular
         return
      end if
      call inverse_matrix(upper, inv_upper, info)
      if (info /= 0) return
      allocate(z(n, k))
      z = matmul(centered, inv_upper)
      do j = 1, k
         b1(j) = sum(z(:, j) ** 3) / real(n, dp)
         b2(j) = sum(z(:, j) ** 4) / real(n, dp)
      end do
      s3 = real(n, dp) * dot_product(b1, b1) / 6.0_dp
      s4 = real(n, dp) * dot_product(b2 - 3.0_dp, b2 - 3.0_dp) / 24.0_dp
      skewness%statistic = s3
      skewness%df1 = real(k, dp)
      skewness%p_value = r_pchisq(s3, real(k, dp), lower_tail = .false.)
      kurtosis%statistic = s4
      kurtosis%df1 = real(k, dp)
      kurtosis%p_value = r_pchisq(s4, real(k, dp), lower_tail = .false.)
      jb%statistic = s3 + s4
      jb%df1 = real(2 * k, dp)
      jb%p_value = r_pchisq(jb%statistic, jb%df1, lower_tail = .false.)
      info = vars_success
   end subroutine jarque_bera_multivariate

   subroutine portmanteau_tests(model, lags, asymptotic, adjusted, info, vec2var_adjustment)
      type(var_model), intent(in) :: model
      integer, intent(in) :: lags
      type(vars_test_result), intent(out) :: asymptotic, adjusted
      integer, intent(out) :: info
      logical, intent(in), optional :: vec2var_adjustment
      real(dp), allocatable :: c0(:, :), c0inv(:, :), ci(:, :)
      real(dp) :: tracesum, adjusted_sum, trace_i, df
      logical :: vec_adjust
      integer :: lag, k, n

      n = model%nobs
      k = model%k
      if (lags < 1 .or. lags >= n) then
         info = vars_invalid_argument
         return
      end if
      allocate(c0(k, k), ci(k, k))
      c0 = matmul(transpose(model%resid), model%resid) / real(n, dp)
      call inverse_matrix(c0, c0inv, info)
      if (info /= 0) then
         info = vars_singular
         return
      end if
      tracesum = 0.0_dp
      adjusted_sum = 0.0_dp
      do lag = 1, lags
         ci = matmul(transpose(model%resid(lag + 1:n, :)), model%resid(1:n - lag, :)) / real(n, dp)
         trace_i = trace_matrix(matmul(transpose(ci), matmul(c0inv, matmul(ci, c0inv))))
         tracesum = tracesum + trace_i
         adjusted_sum = adjusted_sum + trace_i / real(n - lag, dp)
      end do
      vec_adjust = .false.
      if (present(vec2var_adjustment)) vec_adjust = vec2var_adjustment
      df = real(k * k * lags - k * k * model%p, dp)
      if (vec_adjust) df = df + real(k, dp)
      asymptotic%statistic = real(n, dp) * tracesum
      asymptotic%df1 = df
      asymptotic%p_value = r_pchisq(asymptotic%statistic, df, lower_tail = .false.)
      adjusted%statistic = real(n * n, dp) * adjusted_sum
      adjusted%df1 = df
      adjusted%p_value = r_pchisq(adjusted%statistic, df, lower_tail = .false.)
      info = vars_success
   end subroutine portmanteau_tests

   subroutine bg_serial_tests(model, lags, lm_test, es_test, info)
      type(var_model), intent(in) :: model
      integer, intent(in) :: lags
      type(vars_test_result), intent(out) :: lm_test, es_test
      integer, intent(out) :: info
      real(dp), allocatable :: lagres(:, :), x0(:, :), x1(:, :), beta(:)
      real(dp), allocatable :: resid0(:, :), resid1(:, :), sigma0(:, :), sigma1(:, :), inv1(:, :)
      integer, allocatable :: idx(:)
      real(dp) :: logdet0, logdet1, r2r, m, q, nreg_r, big_n, rfactor, base
      integer :: eq, j, n, k, nactive, rank, sign0, sign1

      n = model%nobs
      k = model%k
      if (lags < 1 .or. lags >= n) then
         info = vars_invalid_argument
         return
      end if
      allocate(lagres(n, k * lags))
      lagres = 0.0_dp
      do j = 1, lags
         lagres(j + 1:n, (j - 1) * k + 1:j * k) = model%resid(1:n - j, :)
      end do
      allocate(resid0(n, k), resid1(n, k))
      do eq = 1, k
         nactive = count(model%active(eq, :))
         allocate(idx(nactive), x1(n, nactive), x0(n, nactive + k * lags))
         idx = pack([(j, j = 1, model%nreg)], model%active(eq, :))
         x1 = model%x(:, idx)
         x0(:, 1:nactive) = x1
         x0(:, nactive + 1:nactive + k * lags) = lagres
         allocate(beta(size(x1, 2)))
         call least_squares_svd(x1, model%resid(:, eq), beta, rank, info)
         if (info /= 0) return
         resid1(:, eq) = model%resid(:, eq) - matmul(x1, beta)
         deallocate(beta)
         allocate(beta(size(x0, 2)))
         call least_squares_svd(x0, model%resid(:, eq), beta, rank, info)
         if (info /= 0) return
         resid0(:, eq) = model%resid(:, eq) - matmul(x0, beta)
         deallocate(beta, idx, x1, x0)
      end do
      allocate(sigma0(k, k), sigma1(k, k))
      sigma0 = matmul(transpose(resid0), resid0) / real(n, dp)
      sigma1 = matmul(transpose(resid1), resid1) / real(n, dp)
      call inverse_matrix(sigma1, inv1, info)
      if (info /= 0) return
      lm_test%statistic = real(n, dp) * (real(k, dp) - trace_matrix(matmul(inv1, sigma0)))
      lm_test%df1 = real(lags * k * k, dp)
      lm_test%p_value = r_pchisq(lm_test%statistic, lm_test%df1, lower_tail = .false.)

      call determinant_logabs(sigma0, logdet0, sign0, info)
      if (info /= 0 .or. sign0 <= 0) return
      call determinant_logabs(sigma1, logdet1, sign1, info)
      if (info /= 0 .or. sign1 <= 0) return
      r2r = 1.0_dp - exp(logdet0 - logdet1)
      m = real(k * lags, dp)
      q = 0.5_dp * real(k, dp) * m - 1.0_dp
      nreg_r = real(model%nreg, dp)
      big_n = real(n, dp) - nreg_r - m - 0.5_dp * (real(k, dp) - m + 1.0_dp)
      rfactor = sqrt((real(k * k, dp) * m * m - 4.0_dp) / (real(k * k, dp) + m * m - 5.0_dp))
      base = max(tiny(1.0_dp), 1.0_dp - r2r)
      es_test%statistic = (1.0_dp - base ** (1.0_dp / rfactor)) / base ** (1.0_dp / rfactor) * &
         (big_n * rfactor - q) / (real(k, dp) * m)
      es_test%df1 = real(lags * k * k, dp)
      es_test%df2 = floor(big_n * rfactor - q)
      es_test%p_value = r_pf(es_test%statistic, es_test%df1, es_test%df2, lower_tail = .false.)
      info = vars_success
   end subroutine bg_serial_tests

end module vars_diagnostics
