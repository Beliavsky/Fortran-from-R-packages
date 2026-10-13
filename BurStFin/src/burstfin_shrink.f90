! Public-domain Burns Statistics Ltd. algorithms, translated to modern Fortran.
module burstfin_shrink
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_is_nan, ieee_quiet_nan, ieee_value
   use r_kinds, only: dp
   use r_descriptive, only: r_covariance
   use r_quantiles, only: r_quantile_type7
   use r_linalg, only: symmetric_eigen
   implicit none
   private
   public :: var_shrink_eqcor
contains
   subroutine var_shrink_eqcor(x, covariance, shrinkage, status, weights, shrink, center, center_values, &
         vol_shrink, sd_min, quan_sd, tol, compatible)
      real(dp), intent(in) :: x(:, :) !! Returns by observation and asset; NaNs are pairwise missing, infinities are invalid.
      real(dp), allocatable, intent(out) :: covariance(:, :) !! Equal-correlation-target shrinkage covariance in original units.
      real(dp), intent(out) :: shrinkage !! Applied intensity in [0,1]; initialized to NaN and meaningful only when status is zero.
      integer, intent(out) :: status !! 0 success, 1 invalid inputs, 2 insufficient data, 3 degenerate variance, 4 solver failure.
      real(dp), intent(in), optional :: weights(:) !! Nonnegative row weights; scalar means equal, default linear 0.5 to 1.5.
      real(dp), intent(in), optional :: shrink !! Explicit intensity clipped to [0,1]; absent invokes upstream automatic estimate.
      logical, intent(in), optional :: center !! Use upstream weighted centers when true (default), zero centers otherwise.
      real(dp), intent(in), optional :: center_values(:) !! Explicit finite asset centers, overriding logical centering.
      real(dp), intent(in), optional :: vol_shrink !! Volatility shrinkage toward the mean, capped at one; default zero.
      integer, intent(in), optional :: sd_min !! Minimum observed rows for a reliable volatility, capped at row count; default 20.
      real(dp), intent(in), optional :: quan_sd !! Type-7 volatility replacement quantile in [0,1]; default 0.9.
      real(dp), intent(in), optional :: tol !! Eigenvalue floor as a fraction of the largest eigenvalue; default 0.0001.
      logical, intent(in), optional :: compatible !! Multiply sample covariance by (n-1)/n for Ledoit-Wolf comparison; default false.
      real(dp), allocatable :: z(:, :), w(:), wt(:), y(:, :), sample(:, :), prior(:, :), sd(:), centers(:)
      real(dp), allocatable :: asset_scale(:)
      real(dp), allocatable :: cross(:, :), pi_mat(:, :), theta1(:, :), theta2(:, :), values(:), vectors(:, :)
      logical, allocatable :: keep(:), observed(:, :), good(:), both(:)
      integer, allocatable :: rows(:), nobs(:)
      real(dp) :: vs, qs, tolerance, scale, maximum, meanvol, meancor, gamma, pi_hat, rho, threshold
      integer :: n, p, nr, minimum, i, j, k, pairs, info
      logical :: do_center, compat, had_missing
      status = 1
      shrinkage = ieee_value(0.0_dp, ieee_quiet_nan)
      n = size(x, 1)
      p = size(x, 2)
      if (n < 2 .or. p < 2) return
      if (any(.not. ieee_is_finite(x) .and. .not. ieee_is_nan(x))) return
      vs = 0.0_dp
      qs = 0.9_dp
      tolerance = 0.0001_dp
      minimum = 20
      do_center = .true.
      compat = .false.
      if (present(vol_shrink)) vs = vol_shrink
      if (present(quan_sd)) qs = quan_sd
      if (present(tol)) tolerance = tol
      if (present(sd_min)) minimum = sd_min
      if (present(center)) do_center = center
      if (present(compatible)) compat = compatible
      if (.not. all(ieee_is_finite([vs, qs, tolerance]))) return
      if (minimum < 2 .or. qs < 0.0_dp .or. qs > 1.0_dp) return
      if (present(shrink)) then
         if (.not. ieee_is_finite(shrink)) return
      end if
      if (present(center_values)) then
         if (size(center_values) /= p) return
         if (.not. all(ieee_is_finite(center_values))) return
      end if
      allocate(wt(n))
      if (present(weights)) then
         if (.not. all(ieee_is_finite(weights))) return
         if (any(weights < 0.0_dp)) return
         if (size(weights) == 1) then
            if (weights(1) <= 0.0_dp) return
            wt = 1.0_dp
         else if (size(weights) == n) then
            wt = weights
         else
            return
         end if
      else
         do i = 1, n
            wt(i) = 0.5_dp + real(i - 1, dp) / real(n - 1, dp)
         end do
      end if
      keep = wt > 0.0_dp
      nr = count(keep)
      status = 2
      if (nr < 2) return
      rows = pack([(i, i = 1, n)], keep)
      w = wt(rows)
      w = w / maxval(w)
      w = w / (sum(w) / real(nr, dp))
      z = x(rows, :)
      observed = .not. ieee_is_nan(z)
      had_missing = any(.not. observed)
      nobs = count(observed, dim = 1)
      good = nobs >= min(minimum, nr)
      if (count(good) < 2) return
      where (.not. observed) z = 0.0_dp
      scale = maxval(abs(z))
      if (present(center_values)) scale = max(scale, maxval(abs(center_values)))
      status = 3
      if (scale <= 0.0_dp) return
      allocate(asset_scale(p))
      asset_scale = scale
      if (present(shrink)) then
         ! Explicit shrinkage needs no fourth moments: per-asset scaling preserves mixed dynamic ranges.
         do j = 1, p
            asset_scale(j) = maxval(abs(z(:, j)))
            if (present(center_values)) asset_scale(j) = max(asset_scale(j), abs(center_values(j)))
            if (asset_scale(j) == 0.0_dp) asset_scale(j) = 1.0_dp
         end do
      end if
      do j = 1, p
         z(:, j) = z(:, j) / asset_scale(j)
      end do
      allocate(centers(p))
      centers = 0.0_dp
      if (do_center) then
         do j = 1, p
            if (nobs(j) > 0) centers(j) = sum(z(:, j) * w, mask = observed(:, j)) / real(nobs(j), dp)
         end do
      end if
      if (present(center_values)) centers = center_values / asset_scale
      do j = 1, p
         where (observed(:, j)) z(:, j) = z(:, j) - centers(j)
      end do
      allocate(y(nr, p))
      do j = 1, p
         y(:, j) = z(:, j) * sqrt(w)
      end do
      allocate(sample(p, p), prior(p, p), sd(p))
      do j = 1, p
         do i = 1, j
            both = observed(:, i) .and. observed(:, j)
            sample(i, j) = r_covariance(pack(y(:, i), both), pack(y(:, j), both))
            sample(j, i) = sample(i, j)
         end do
      end do
      if (compat) sample = sample * real(nr - 1, dp) / real(nr, dp)
      do j = 1, p
         sd(j) = sqrt(max(0.0_dp, sample(j, j)))
         if (ieee_is_nan(sample(j, j))) sd(j) = ieee_value(0.0_dp, ieee_quiet_nan)
      end do
      if (any(pack(sd, good) <= 0.0_dp)) return
      meancor = 0.0_dp
      pairs = 0
      do j = 1, p
         if (.not. good(j)) cycle
         do i = j + 1, p
            if (.not. good(i)) cycle
            if (.not. ieee_is_finite(sample(i, j))) return
            ! Correlation target uses the original sample diagonal, before volatility shrinkage.
            meancor = meancor + sample(i, j) / sqrt(sample(i, i)) / sqrt(sample(j, j))
            pairs = pairs + 1
         end do
      end do
      if (pairs < 1) return
      meancor = meancor / real(pairs, dp)
      if (present(shrink)) then
         sd = sd * asset_scale
         do j = 1, p
            do i = 1, p
               sample(i, j) = (sample(i, j) * asset_scale(i)) * asset_scale(j)
            end do
         end do
         scale = 1.0_dp
      end if
      maximum = r_quantile_type7(pack(sd, good), qs)
      where (.not. good) sd = maximum
      if (vs > 0.0_dp) then
         vs = min(vs, 1.0_dp)
         meanvol = sum(sd) / real(p, dp)
         sd = vs * meanvol + (1.0_dp - vs) * sd
      end if
      do j = 1, p
         do i = 1, p
            prior(i, j) = meancor * sd(i) * sd(j)
         end do
         prior(j, j) = sd(j) ** 2
      end do
      where (ieee_is_nan(sample)) sample = prior
      if (present(shrink)) then
         shrinkage = min(1.0_dp, max(0.0_dp, shrink))
      else
         gamma = sum((sample - prior) ** 2)
         if (gamma <= tiny(1.0_dp)) then
            ! The estimator is independent of intensity when its target equals the sample.
            shrinkage = 0.0_dp
         else
            allocate(pi_mat(p, p), theta1(p, p), theta2(p, p), cross(p, p))
            pi_mat = 0.0_dp
            theta1 = 0.0_dp
            theta2 = 0.0_dp
            do k = 1, nr
               do j = 1, p
                  do i = 1, p
                     cross(i, j) = z(k, i) * z(k, j) - sample(i, j)
                  end do
               end do
               pi_mat = pi_mat + w(k) * cross ** 2
               do j = 1, p
                  do i = 1, p
                     theta1(i, j) = theta1(i, j) + w(k) * cross(i, i) * cross(i, j)
                     theta2(i, j) = theta2(i, j) + w(k) * cross(j, j) * cross(i, j)
                  end do
               end do
            end do
            pi_hat = sum(pi_mat) / real(nr, dp)
            rho = 0.0_dp
            do j = 1, p
               rho = rho + pi_mat(j, j) / real(nr, dp)
               do i = 1, p
                  if (i == j) cycle
                  rho = rho + meancor * 0.5_dp * &
                     (theta1(i, j) * sd(j) / sd(i) + theta2(i, j) * sd(i) / sd(j)) / real(nr, dp)
               end do
            end do
            shrinkage = min(1.0_dp, max(0.0_dp, (pi_hat - rho) / gamma / real(nr, dp)))
         end if
      end if
      covariance = shrinkage * prior + (1.0_dp - shrinkage) * sample
      if (tolerance > 0.0_dp .or. (had_missing .and. present(shrink))) then
         call symmetric_eigen(covariance, values, vectors, info)
         status = 4
         if (info /= 0) return
         threshold = tolerance * maxval(values)
         if (minval(values) < threshold) then
            values = max(values, threshold)
            covariance = 0.0_dp
            do k = 1, p
               do j = 1, p
                  covariance(:, j) = covariance(:, j) + vectors(:, k) * values(k) * vectors(j, k)
               end do
            end do
         end if
      end if
      covariance = (covariance * scale) * scale
      status = 3
      if (.not. all(ieee_is_finite(covariance))) return
      status = 0
   end subroutine var_shrink_eqcor
end module burstfin_shrink
