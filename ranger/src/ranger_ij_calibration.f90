! SPDX-License-Identifier: GPL-3.0-only
! Modern Fortran translation of ranger 0.18.0 computational code.
! The empirical-Bayes calibration follows R/infinitesimalJackknife.R,
! adapted from randomForestCI by Stefan Wager. See NOTICE.md.
module ranger_ij_calibration
   use r_kinds, only : dp, i64
   use ranger_rng, only : ranger_rng_state, sample_indices
   implicit none
   private

   integer, parameter :: eb_nbin = 1000
   real(dp), parameter :: eb_unif_fraction = 0.1_dp
   real(dp), parameter :: sqrt_two_pi = 2.506628274631000502415765284811_dp

   public :: calibrate_ij_variances

contains

   subroutine calibrate_ij_variances(pred, inbag, variance_full, variance_calibrated, seed, status)
      real(dp), intent(in) :: pred(:,:)
      integer, intent(in) :: inbag(:,:)
      real(dp), intent(in) :: variance_full(:)
      real(dp), intent(out) :: variance_calibrated(size(variance_full))
      integer(i64), intent(in), optional :: seed
      integer, intent(out), optional :: status
      type(ranger_rng_state) :: rng
      integer, allocatable :: chosen(:)
      real(dp), allocatable :: variance_sub(:)
      integer :: ntree, nsample, rc
      real(dp) :: sigma2_ss, delta, sigma2
      integer(i64) :: local_seed

      if (present(status)) status = 0
      ! rInfJack disables calibration entirely for <= 20 prediction
      ! points, returning the raw (possibly negative) IJ variances unchanged.
      variance_calibrated = variance_full
      if (size(variance_full) <= 20) return
      ntree = size(pred, 2)
      if (size(inbag, 2) /= ntree .or. ntree <= 1) then
         if (present(status)) status = 1
         return
      end if

      nsample = (ntree + 1) / 2
      allocate(chosen(nsample), variance_sub(size(variance_full)))
      local_seed = 1976_i64
      if (present(seed)) local_seed = seed
      call rng%seed(local_seed)
      call sample_indices(rng, ntree, nsample, .false., chosen, status=rc)
      if (rc /= 0) then
         if (present(status)) status = 2
         return
      end if
      call ij_variance_subset(pred, inbag, chosen, variance_sub)
      sigma2_ss = sum((variance_sub - variance_full) ** 2) / real(size(variance_full), dp)
      delta = real(nsample, dp) / real(ntree, dp)
      sigma2 = (delta ** 2 + (1.0_dp - delta) ** 2) / (2.0_dp * (1.0_dp - delta) ** 2) * sigma2_ss
      call calibrate_eb(variance_full, sigma2, variance_calibrated, rc)
      if (rc /= 0) then
         variance_calibrated = variance_full
         if (present(status)) status = 3
      end if
   end subroutine calibrate_ij_variances

   subroutine ij_variance_subset(pred, inbag, tree_index, variance)
      real(dp), intent(in) :: pred(:,:)
      integer, intent(in) :: inbag(:,:), tree_index(:)
      real(dp), intent(out) :: variance(size(pred, 1))
      real(dp), allocatable :: centered(:,:), navg(:), covariance(:,:)
      real(dp) :: mean_pred, nvar, bootvar, sample_fraction, inflation
      integer :: ntrain, ntest, bcount, i, j, b, tree
      logical :: no_replacement

      ntrain = size(inbag, 1)
      ntest = size(pred, 1)
      bcount = size(tree_index)
      allocate(centered(ntest, bcount), navg(ntrain), covariance(ntrain, ntest))
      do j = 1, ntest
         mean_pred = sum(pred(j, tree_index)) / real(bcount, dp)
         do b = 1, bcount
            centered(j, b) = pred(j, tree_index(b)) - mean_pred
         end do
      end do
      navg = 0.0_dp
      do b = 1, bcount
         tree = tree_index(b)
         navg = navg + real(inbag(:, tree), dp)
      end do
      navg = navg / real(bcount, dp)
      covariance = 0.0_dp
      do i = 1, ntrain
         do j = 1, ntest
            do b = 1, bcount
               tree = tree_index(b)
               covariance(i, j) = covariance(i, j) + &
                  (real(inbag(i, tree), dp) - navg(i)) * centered(j, b)
            end do
            covariance(i, j) = covariance(i, j) / real(bcount, dp)
         end do
      end do
      variance = sum(covariance * covariance, dim=1)

      nvar = 0.0_dp
      do i = 1, ntrain
         do b = 1, bcount
            tree = tree_index(b)
            nvar = nvar + real(inbag(i, tree), dp) ** 2 / real(bcount, dp)
         end do
         nvar = nvar - navg(i) ** 2
      end do
      nvar = nvar / real(ntrain, dp)
      do j = 1, ntest
         bootvar = sum(centered(j, :) ** 2) / real(bcount, dp)
         variance(j) = variance(j) - real(ntrain, dp) * nvar * bootvar / real(bcount, dp)
      end do

      no_replacement = maxval(inbag) == 1
      if (no_replacement) then
         sample_fraction = sum(real(inbag, dp)) / real(size(inbag, 1) * size(inbag, 2), dp)
         if (sample_fraction < 1.0_dp) then
            inflation = 1.0_dp / (1.0_dp - sample_fraction) ** 2
            variance = inflation * variance
         end if
      end if
   end subroutine ij_variance_subset

   subroutine calibrate_eb(vars, sigma2, calibrated, status)
      real(dp), intent(in) :: vars(:), sigma2
      real(dp), intent(out) :: calibrated(size(vars))
      integer, intent(out) :: status
      real(dp), allocatable :: gx(:), gg(:), qx(:), qy(:)
      real(dp) :: sigma
      integer :: i, nq

      status = 0
      if (sigma2 <= 0.0_dp .or. exact_equal(minval(vars), maxval(vars))) then
         calibrated = max(vars, 0.0_dp)
         return
      end if
      sigma = sqrt(sigma2)
      call gfit(vars, sigma, gx, gg, status)
      if (status /= 0) return

      ! Upstream ranger 0.18.0 tests length(vars >= 200), which is nonzero
      ! for every nonempty vector. Preserve that behavior exactly: calibration
      ! always uses the 2%-spaced quantile interpolation path.
      allocate(qx(51), qy(51))
      do i = 1, 51
         qx(i) = quantile_type7(vars, real(i - 1, dp) / 50.0_dp)
      end do
      call unique_sorted(qx, nq)
      do i = 1, nq
         qy(i) = gbayes(qx(i), gx, gg, sigma)
      end do
      do i = 1, size(vars)
         calibrated(i) = linear_interp(qx(1:nq), qy(1:nq), vars(i))
      end do
      calibrated = max(calibrated, 0.0_dp)
   end subroutine calibrate_eb

   subroutine gfit(x, sigma, xvals, g, status)
      real(dp), intent(in) :: x(:), sigma
      real(dp), allocatable, intent(out) :: xvals(:), g(:)
      integer, intent(out) :: status
      real(dp), allocatable :: noise(:), noise_rotate(:)
      real(dp) :: sx, lo, hi, binw, eta(2)
      integer :: i, zero_idx

      status = 0
      if (sigma <= 0.0_dp .or. size(x) < 2) then
         status = 1
         allocate(xvals(0), g(0))
         return
      end if
      sx = sample_sd(x)
      lo = min(minval(x) - 2.0_dp * sx, 0.0_dp)
      hi = max(maxval(x) + 2.0_dp * sx, sx)
      if (hi <= lo) then
         status = 2
         allocate(xvals(0), g(0))
         return
      end if
      allocate(xvals(eb_nbin), g(eb_nbin), noise(eb_nbin), noise_rotate(eb_nbin))
      do i = 1, eb_nbin
         xvals(i) = lo + real(i - 1, dp) * (hi - lo) / real(eb_nbin - 1, dp)
      end do
      binw = xvals(2) - xvals(1)
      zero_idx = maxloc(merge([(i, i = 1, eb_nbin)], 0, xvals <= 0.0_dp), dim=1)
      do i = 1, eb_nbin
         noise(i) = normal_density(xvals(i) / sigma) * binw / sigma
      end do
      noise_rotate(1:eb_nbin - zero_idx + 1) = noise(zero_idx:eb_nbin)
      if (zero_idx > 1) noise_rotate(eb_nbin - zero_idx + 2:eb_nbin) = noise(1:zero_idx - 1)

      eta = -1.0_dp
      call minimize_gfit(x, xvals, noise_rotate, eta)
      call make_prior(xvals, eta, g, status)
   end subroutine gfit

   subroutine minimize_gfit(x, xvals, noise_rotate, eta)
      real(dp), intent(in) :: x(:), xvals(:), noise_rotate(:)
      real(dp), intent(inout) :: eta(2)
      real(dp) :: f, fnew, grad(2), hess(2,2), step(2), candidate(2), scale, det
      real(dp) :: h1, h2, fpp, fpm, fmp, fmm, fp, fm
      integer :: iter, ls

      f = gfit_objective(eta, x, xvals, noise_rotate)
      do iter = 1, 100
         h1 = 1.0e-4_dp * max(1.0_dp, abs(eta(1)))
         h2 = 1.0e-4_dp * max(1.0_dp, abs(eta(2)))
         fp = gfit_objective(eta + [h1, 0.0_dp], x, xvals, noise_rotate)
         fm = gfit_objective(eta - [h1, 0.0_dp], x, xvals, noise_rotate)
         grad(1) = (fp - fm) / (2.0_dp * h1)
         hess(1,1) = (fp - 2.0_dp * f + fm) / h1 ** 2
         fp = gfit_objective(eta + [0.0_dp, h2], x, xvals, noise_rotate)
         fm = gfit_objective(eta - [0.0_dp, h2], x, xvals, noise_rotate)
         grad(2) = (fp - fm) / (2.0_dp * h2)
         hess(2,2) = (fp - 2.0_dp * f + fm) / h2 ** 2
         fpp = gfit_objective(eta + [h1, h2], x, xvals, noise_rotate)
         fpm = gfit_objective(eta + [h1, -h2], x, xvals, noise_rotate)
         fmp = gfit_objective(eta + [-h1, h2], x, xvals, noise_rotate)
         fmm = gfit_objective(eta - [h1, h2], x, xvals, noise_rotate)
         hess(1,2) = (fpp - fpm - fmp + fmm) / (4.0_dp * h1 * h2)
         hess(2,1) = hess(1,2)
         if (maxval(abs(grad)) <= 1.0e-6_dp * max(1.0_dp, abs(f))) exit

         det = hess(1,1) * hess(2,2) - hess(1,2) * hess(2,1)
         if (det > 1.0e-12_dp .and. hess(1,1) > 0.0_dp .and. hess(2,2) > 0.0_dp) then
            step(1) = -(hess(2,2) * grad(1) - hess(1,2) * grad(2)) / det
            step(2) = -(-hess(2,1) * grad(1) + hess(1,1) * grad(2)) / det
         else
            step = -grad / max(1.0_dp, sqrt(sum(grad ** 2)))
         end if
         if (sqrt(sum(step ** 2)) > 1000.0_dp) step = step * (1000.0_dp / sqrt(sum(step ** 2)))
         scale = 1.0_dp
         do ls = 1, 40
            candidate = eta + scale * step
            fnew = gfit_objective(candidate, x, xvals, noise_rotate)
            if (fnew < f) exit
            scale = 0.5_dp * scale
         end do
         if (fnew >= f) exit
         if (maxval(abs(candidate - eta) / max(1.0_dp, abs(eta))) <= 1.0e-6_dp) then
            eta = candidate
            exit
         end if
         eta = candidate
         f = fnew
      end do
   end subroutine minimize_gfit

   real(dp) function gfit_objective(eta, x, xvals, noise_rotate) result(value)
      real(dp), intent(in) :: eta(2), x(:), xvals(:), noise_rotate(:)
      real(dp), allocatable :: g(:), f(:), nll(:)
      integer :: rc, i

      allocate(g(size(xvals)), f(size(xvals)), nll(size(xvals)))
      call make_prior(xvals, eta, g, rc)
      if (rc /= 0) then
         value = 1000.0_dp * (real(size(x), dp) + sum(eta ** 2))
         return
      end if
      call circular_convolve(g, noise_rotate, f)
      do i = 1, size(f)
         nll(i) = -log(max(f(i), 1.0e-7_dp))
      end do
      value = 0.0_dp
      do i = 1, size(x)
         value = value + linear_interp(xvals, nll, x(i))
      end do
   end function gfit_objective

   subroutine make_prior(xvals, eta, g, status)
      real(dp), intent(in) :: xvals(:), eta(2)
      real(dp), intent(out) :: g(size(xvals))
      integer, intent(out) :: status
      real(dp) :: exponent, total
      integer :: i, nnonneg

      status = 0
      g = 0.0_dp
      nnonneg = count(xvals >= 0.0_dp)
      do i = 1, size(xvals)
         if (xvals(i) < 0.0_dp) cycle
         exponent = eta(1) * xvals(i) + eta(2) * xvals(i) ** 2
         if (exponent >= log(huge(1.0_dp)) - 2.0_dp) then
            status = 1
            return
         end if
         g(i) = exp(exponent)
      end do
      total = sum(g)
      if (total <= 100.0_dp * epsilon(1.0_dp)) then
         status = 2
         return
      end if
      g = (1.0_dp - eb_unif_fraction) * g / total
      where (xvals >= 0.0_dp)
         g = g + eb_unif_fraction / real(nnonneg, dp)
      end where
   end subroutine make_prior

   subroutine circular_convolve(x, y, result)
      real(dp), intent(in) :: x(:), y(:)
      real(dp), intent(out) :: result(size(x))
      integer :: n, k, i, j

      n = size(x)
      result = 0.0_dp
      do k = 1, n
         do i = 1, n
            j = k - n + i
            if (j < 1) j = j + n
            result(k) = result(k) + x(j) * y(i)
         end do
      end do
   end subroutine circular_convolve

   real(dp) function gbayes(x0, gx, gg, sigma) result(value)
      real(dp), intent(in) :: x0, gx(:), gg(:), sigma
      real(dp) :: w, total, weighted
      integer :: i

      total = 0.0_dp
      weighted = 0.0_dp
      do i = 1, size(gx)
         w = normal_density((gx(i) - x0) / sigma) * gg(i)
         total = total + w
         weighted = weighted + w * gx(i)
      end do
      if (total > 0.0_dp) then
         value = weighted / total
      else
         value = max(x0, 0.0_dp)
      end if
   end function gbayes

   pure real(dp) function normal_density(z) result(value)
      real(dp), intent(in) :: z
      value = exp(-0.5_dp * z ** 2) / sqrt_two_pi
   end function normal_density

   pure real(dp) function sample_sd(x) result(value)
      real(dp), intent(in) :: x(:)
      real(dp) :: mean_x
      if (size(x) <= 1) then
         value = 0.0_dp
      else
         mean_x = sum(x) / real(size(x), dp)
         value = sqrt(sum((x - mean_x) ** 2) / real(size(x) - 1, dp))
      end if
   end function sample_sd

   real(dp) function quantile_type7(x, probability) result(value)
      real(dp), intent(in) :: x(:), probability
      real(dp), allocatable :: work(:)
      real(dp) :: h, frac, key
      integer :: i, j, lo

      allocate(work(size(x)))
      work = x
      do i = 2, size(work)
         key = work(i)
         j = i - 1
         do while (j >= 1)
            if (work(j) <= key) exit
            work(j + 1) = work(j)
            j = j - 1
         end do
         work(j + 1) = key
      end do
      h = 1.0_dp + real(size(work) - 1, dp) * min(1.0_dp, max(0.0_dp, probability))
      lo = int(floor(h))
      if (lo >= size(work)) then
         value = work(size(work))
      else
         frac = h - real(lo, dp)
         value = (1.0_dp - frac) * work(lo) + frac * work(lo + 1)
      end if
   end function quantile_type7

   subroutine unique_sorted(x, n)
      real(dp), intent(inout) :: x(:)
      integer, intent(out) :: n
      real(dp) :: key
      integer :: i, j

      do i = 2, size(x)
         key = x(i)
         j = i - 1
         do while (j >= 1)
            if (x(j) <= key) exit
            x(j + 1) = x(j)
            j = j - 1
         end do
         x(j + 1) = key
      end do
      n = 1
      do i = 2, size(x)
         if (.not. exact_equal(x(i), x(n))) then
            n = n + 1
            x(n) = x(i)
         end if
      end do
   end subroutine unique_sorted

   pure real(dp) function linear_interp(x, y, x0) result(value)
      real(dp), intent(in) :: x(:), y(:), x0
      integer :: lo, hi, mid
      real(dp) :: weight

      if (x0 <= x(1)) then
         value = y(1)
         return
      end if
      if (x0 >= x(size(x))) then
         value = y(size(y))
         return
      end if
      lo = 1
      hi = size(x)
      do while (hi - lo > 1)
         mid = (lo + hi) / 2
         if (x0 >= x(mid)) then
            lo = mid
         else
            hi = mid
         end if
      end do
      if (exact_equal(x(hi), x(lo))) then
         value = y(lo)
      else
         weight = (x0 - x(lo)) / (x(hi) - x(lo))
         value = (1.0_dp - weight) * y(lo) + weight * y(hi)
      end if
   end function linear_interp

   pure logical function exact_equal(a, b) result(equal)
      real(dp), intent(in) :: a, b
      equal = abs(a - b) <= 0.0_dp
   end function exact_equal

end module ranger_ij_calibration
