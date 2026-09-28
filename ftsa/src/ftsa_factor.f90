! SPDX-License-Identifier: GPL-3.0-only
! Derived from R package ftsa 6.7 by Rob Hyndman and Han Lin Shang.
module ftsa_factor
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   use r_kinds, only : dp
   use r_linalg, only : symmetric_eigen
   use ftsa_types, only : hurst_result
   implicit none
   private

   public :: ftsa_er_gr, ftsa_h_estimate_rs

contains

   subroutine ftsa_er_gr(data, k_er, k_gr, eigenvalues, status)
      real(dp), intent(in) :: data(:, :) !! Observation-by-variable matrix used by upstream `ER_GR`.
      integer, intent(out) :: k_er !! Eigenvalue-ratio selected component count.
      integer, intent(out) :: k_gr !! Growth-ratio selected component count.
      real(dp), allocatable, intent(out) :: eigenvalues(:) !! Descending eigenvalues of `transpose(data)*data/(n*p)`.
      integer, intent(out) :: status !! Zero on success; nonzero when the selection range is too small or eigensolution fails.
      real(dp), allocatable :: covar(:, :), vectors(:, :), er(:), gr(:)
      real(dp) :: u1, u2, v0, v1, v2
      integer :: info, k, kmax, kmax_star, m, n, p

      n = size(data, 1)
      p = size(data, 2)
      m = min(n, p)
      k_er = 0
      k_gr = 0
      if (n < 2 .or. p < 2) then
         allocate(eigenvalues(0))
         status = 1
         return
      end if
      allocate(covar(p, p))
      covar = matmul(transpose(data), data) / real(n * p, dp)
      call symmetric_eigen(covar, eigenvalues, vectors, info, descending=.true.)
      if (info /= 0) then
         status = info
         return
      end if
      v0 = sum(eigenvalues(1:m))
      kmax_star = count(eigenvalues(1:m) >= v0 / real(m, dp))
      kmax = min(kmax_star, max(1, int(0.1_dp * real(m, dp))))
      if (kmax < 2) then
         status = 2
         return
      end if
      allocate(er(kmax - 1), gr(kmax - 1))
      do k = 1, kmax - 1
         if (abs(eigenvalues(k + 1)) <= tiny(1.0_dp)) then
            er(k) = huge(1.0_dp)
         else
            er(k) = eigenvalues(k) / eigenvalues(k + 1)
         end if
         v1 = sum(eigenvalues(k + 1:m))
         v2 = sum(eigenvalues(k + 2:m))
         if (v1 <= tiny(1.0_dp) .or. v2 <= tiny(1.0_dp)) then
            gr(k) = -huge(1.0_dp)
         else
            u1 = eigenvalues(k) / v1
            u2 = eigenvalues(k + 1) / v2
            if (log(1.0_dp + u2) <= tiny(1.0_dp)) then
               gr(k) = -huge(1.0_dp)
            else
               gr(k) = log(1.0_dp + u1) / log(1.0_dp + u2)
            end if
         end if
      end do
      k_er = maxloc(er, dim=1)
      k_gr = maxloc(gr, dim=1)
      status = 0
   end subroutine ftsa_er_gr

   subroutine ftsa_h_estimate_rs(fun_dat, result, band_const, average_scores, modified, status)
      real(dp), intent(in) :: fun_dat(:, :) !! Grid-by-time functional observations.
      type(hurst_result), intent(out) :: result !! Long-run covariance eigensystem and rescaled-range estimates.
      real(dp), intent(in), optional :: band_const !! Multiplier applied to the upstream lag span; default one.
      logical, intent(in), optional :: average_scores !! Average scores over 95%-variance components instead of first PC.
      logical, intent(in), optional :: modified !! Use Lo-style modified rescaled range when true.
      integer, intent(out) :: status !! Zero on success; nonzero for invalid dimensions or eigensolver failure.
      real(dp), allocatable :: centered(:, :), eigenvalues(:), vectors(:, :), gamma(:, :), score(:), partial(:)
      real(dp) :: band, mean_score, r_n, rho1, s_n, variance, weight_sum
      integer :: info, k, lag, m, n, ncomp, p, q, t
      logical :: average, use_modified

      p = size(fun_dat, 1)
      n = size(fun_dat, 2)
      band = 1.0_dp
      if (present(band_const)) band = band_const
      average = .false.
      if (present(average_scores)) average = average_scores
      use_modified = .false.
      if (present(modified)) use_modified = modified
      if (p < 1 .or. n < 3 .or. band <= 0.0_dp) then
         result%info = 1
         status = 1
         return
      end if
      m = max(1, int(real(min(n - 1, p), dp) * band))
      allocate(centered(p, n), result%covariance(p, p), gamma(p, p))
      do k = 1, p
         centered(k, :) = fun_dat(k, :) - sum(fun_dat(k, :)) / real(n, dp)
      end do
      result%covariance = 0.0_dp
      do lag = -m, m
         call lag_covariance(centered, lag, gamma)
         result%covariance = result%covariance + real(m - abs(lag), dp) * gamma
      end do
      call symmetric_eigen(result%covariance, eigenvalues, vectors, info, descending=.true.)
      if (info /= 0 .or. sum(max(eigenvalues, 0.0_dp)) <= tiny(1.0_dp)) then
         result%info = merge(info, 2, info /= 0)
         status = result%info
         return
      end if
      allocate(result%prop_eigenvalue(size(eigenvalues)))
      result%prop_eigenvalue = cumulative_positive_fraction(eigenvalues)
      ncomp = 1
      do k = 1, size(result%prop_eigenvalue)
         if (result%prop_eigenvalue(k) >= 0.95_dp) then
            ncomp = k
            exit
         end if
      end do
      result%ncomp = ncomp
      if (average .and. ncomp > 1) then
         allocate(result%eigenfunctions(p, ncomp), score(n))
         result%eigenfunctions = vectors(:, 1:ncomp)
         score = sum(matmul(transpose(result%eigenfunctions), fun_dat), dim=1) / real(ncomp, dp)
      else
         allocate(result%eigenfunctions(p, 1), score(n))
         result%eigenfunctions(:, 1) = vectors(:, 1)
         score = matmul(result%eigenfunctions(:, 1), fun_dat)
      end if
      allocate(result%score(n), partial(n))
      result%score = score
      mean_score = sum(score) / real(n, dp)
      partial(1) = score(1) - mean_score
      do t = 2, n
         partial(t) = partial(t - 1) + score(t) - mean_score
      end do
      r_n = maxval(partial) - minval(partial)
      variance = sum((score - mean_score)**2) / real(n, dp)
      if (variance <= tiny(1.0_dp) .or. r_n <= tiny(1.0_dp)) then
         result%h_est = 0.0_dp
         result%d_est = -0.5_dp
         result%alpha_est = 1.5_dp
         result%info = 2
         status = 2
         return
      end if
      if (.not. use_modified) then
         s_n = sqrt(variance)
         result%h_est = log(r_n / s_n) / log(real(n, dp))
      else
         rho1 = lag_correlation(score, 1)
         q = int((1.5_dp * real(n, dp))**(1.0_dp / 3.0_dp) * &
            abs(2.0_dp * rho1 / max(tiny(1.0_dp), 1.0_dp - rho1 * rho1))**(2.0_dp / 3.0_dp))
         q = min(max(0, q), n - 2)
         weight_sum = 0.0_dp
         do k = 1, q
            weight_sum = weight_sum + (1.0_dp - real(k, dp) / real(q + 1, dp)) * lag_covariance_scalar(score, k)
         end do
         s_n = sqrt(max(tiny(1.0_dp), variance + 2.0_dp * weight_sum))
         result%h_est = r_n / (s_n * sqrt(real(n, dp)))
      end if
      result%d_est = result%h_est - 0.5_dp
      result%alpha_est = 1.5_dp - result%h_est
      result%info = 0
      status = 0
   end subroutine ftsa_h_estimate_rs

   pure subroutine lag_covariance(centered, lag, covariance)
      real(dp), intent(in) :: centered(:, :) !! Grid-by-time centered functional observations.
      integer, intent(in) :: lag !! Signed time lag.
      real(dp), intent(out) :: covariance(:, :) !! Lag covariance divided by full series length.
      integer :: h, n, p, t

      p = size(centered, 1)
      n = size(centered, 2)
      h = abs(lag)
      covariance = 0.0_dp
      if (lag >= 0) then
         do t = 1, n - h
            covariance = covariance + spread(centered(:, t), 2, p) * spread(centered(:, t + h), 1, p)
         end do
      else
         do t = 1, n - h
            covariance = covariance + spread(centered(:, t + h), 2, p) * spread(centered(:, t), 1, p)
         end do
      end if
      covariance = covariance / real(n, dp)
   end subroutine lag_covariance

   pure function cumulative_positive_fraction(values) result(fraction)
      real(dp), intent(in) :: values(:) !! Descending eigenvalues.
      real(dp) :: fraction(size(values))
      real(dp) :: total
      integer :: i

      total = sum(max(values, 0.0_dp))
      fraction = 0.0_dp
      if (total <= tiny(1.0_dp)) return
      fraction(1) = max(values(1), 0.0_dp) / total
      do i = 2, size(values)
         fraction(i) = fraction(i - 1) + max(values(i), 0.0_dp) / total
      end do
   end function cumulative_positive_fraction

   pure real(dp) function lag_covariance_scalar(x, lag) result(value)
      real(dp), intent(in) :: x(:) !! Scalar time series.
      integer, intent(in) :: lag !! Positive lag less than series length.
      real(dp) :: mean_x
      integer :: n

      n = size(x)
      mean_x = sum(x) / real(n, dp)
      value = sum((x(1:n-lag) - mean_x) * (x(1+lag:n) - mean_x)) / real(n, dp)
   end function lag_covariance_scalar

   pure real(dp) function lag_correlation(x, lag) result(value)
      real(dp), intent(in) :: x(:) !! Scalar time series.
      integer, intent(in) :: lag !! Positive lag less than series length.
      real(dp) :: variance

      variance = lag_covariance_scalar(x, 0)
      if (variance <= tiny(1.0_dp)) then
         value = 0.0_dp
      else
         value = lag_covariance_scalar(x, lag) / variance
      end if
   end function lag_correlation

end module ftsa_factor
