! SPDX-License-Identifier: AGPL-3.0-only
! Derived from glmmTMB 1.1.14 computational sources; see NOTICE.md.
module glmmtmb_covariance
   use glmmtmb_kinds, only: dp
   use glmmtmb_codes
   use glmmtmb_math, only: invlogit, logaddexp, logcosh_safe
   use tmb_distributions, only: dnorm
   use tmb_numerics, only: mvnorm_nll, n01_nll, unstructured_corr
   use, intrinsic :: ieee_arithmetic, only: ieee_is_nan, ieee_quiet_nan, ieee_value
   implicit none
   private
   public :: covariance_term_nll, matern_corr, reduced_rank_loadings, reduced_rank_transform
contains
   pure real(dp) function scaled_mvn_nll(x, corr, sd) result(ans)
      real(dp), intent(in) :: x(:) !! Zero-mean random-effect vector on its data scale.
      real(dp), intent(in) :: corr(:, :) !! Positive-definite correlation matrix matching x.
      real(dp), intent(in) :: sd(:) !! Positive marginal standard deviations matching x.
      real(dp) :: sigma(size(x), size(x))
      integer :: i, j, n
      n = size(x)
      if (size(corr, 1) /= n .or. size(corr, 2) /= n .or. size(sd) /= n .or. any(sd <= 0.0_dp)) then
         ans = ieee_value(ans, ieee_quiet_nan)
         return
      end if
      do j = 1, n
         do i = 1, n
            sigma(i, j) = corr(i, j) * sd(i) * sd(j)
         end do
      end do
      ans = mvnorm_nll(x, sigma)
   end function scaled_mvn_nll

   pure real(dp) function log_bessel_k_integral(nu, x) result(ans)
      real(dp), intent(in) :: nu !! Positive order of the modified Bessel K function.
      real(dp), intent(in) :: x !! Strictly positive Bessel argument.
      integer, parameter :: nq = 4096
      real(dp) :: h, t, logterm, logsum, tmax
      integer :: i
      if (nu <= 0.0_dp .or. x <= 0.0_dp) then
         ans = ieee_value(ans, ieee_quiet_nan)
         return
      end if
      tmax = max(12.0_dp, min(30.0_dp, log(1.0_dp + 2.0_dp * nu / x) + 6.0_dp))
      h = tmax / real(nq, dp)
      logsum = -huge(1.0_dp)
      do i = 1, nq
         t = (real(i, dp) - 0.5_dp) * h
         logterm = -x * cosh(t) + logcosh_safe(nu * t)
         logsum = logaddexp(logsum, logterm)
      end do
      ans = logsum + log(h)
   end function log_bessel_k_integral

   pure real(dp) function matern_corr(distance, range, smoothness) result(ans)
      real(dp), intent(in) :: distance !! Nonnegative pairwise distance between two random-effect locations.
      real(dp), intent(in) :: range !! Strictly positive Matérn range parameter phi.
      real(dp), intent(in) :: smoothness !! Strictly positive Matérn smoothness parameter kappa.
      real(dp) :: x, logk, logcorr
      if (distance < 0.0_dp .or. range <= 0.0_dp .or. smoothness <= 0.0_dp) then
         ans = ieee_value(ans, ieee_quiet_nan)
      else if (distance == 0.0_dp) then
         ans = 1.0_dp
      else
         x = distance / range
         logk = log_bessel_k_integral(smoothness, x)
         logcorr = -log_gamma(smoothness) - (smoothness - 1.0_dp) * log(2.0_dp)
         logcorr = logcorr + smoothness * log(x) + logk
         ans = min(1.0_dp, exp(logcorr))
      end if
   end function matern_corr

   pure subroutine reduced_rank_loadings(theta, p, lambda, rank, status)
      real(dp), intent(in) :: theta(:) !! Reduced-rank loading parameters in glmmTMB's diagonal-then-lower ordering.
      integer, intent(in) :: p !! Number of rows in the factor-loading matrix.
      real(dp), allocatable, intent(out) :: lambda(:, :) !! Constructed p-by-rank lower-trapezoidal factor-loading matrix.
      integer, intent(out) :: rank !! Rank inferred from p and the number of parameters.
      integer, intent(out) :: status !! Zero on success, nonzero when the parameter count does not imply a valid rank.
      integer :: i, j, k, nt
      real(dp) :: disc
      nt = size(theta)
      disc = real((2 * p + 1) * (2 * p + 1) - 8 * nt, dp)
      if (p <= 0 .or. disc < 0.0_dp) then
         allocate(lambda(0, 0))
         rank = 0
         status = 1
         return
      end if
      rank = nint(0.5_dp * (real(2 * p + 1, dp) - sqrt(disc)))
      if (rank < 1 .or. rank > p .or. nt /= rank * (2 * p - rank + 1) / 2) then
         allocate(lambda(0, 0))
         status = 2
         return
      end if
      allocate(lambda(p, rank))
      lambda = 0.0_dp
      do j = 1, rank
         lambda(j, j) = theta(j)
      end do
      k = rank
      do j = 1, rank
         do i = j + 1, p
            k = k + 1
            lambda(i, j) = theta(k)
         end do
      end do
      status = 0
   end subroutine reduced_rank_loadings

   pure subroutine reduced_rank_transform(u_spherical, lambda, b_data, status)
      real(dp), intent(in) :: u_spherical(:, :) !! Spherical random effects, with at least rank rows and replicates in columns.
      real(dp), intent(in) :: lambda(:, :) !! p-by-rank reduced-rank factor-loading matrix.
      real(dp), intent(out) :: b_data(:, :) !! Data-scale random effects lambda times the first rank spherical rows.
      integer, intent(out) :: status !! Zero on success, nonzero for incompatible matrix dimensions.
      integer :: rank
      rank = size(lambda, 2)
      if (size(u_spherical, 1) < rank .or. size(b_data, 1) /= size(lambda, 1) .or. &
          size(b_data, 2) /= size(u_spherical, 2)) then
         b_data = 0.0_dp
         status = 1
         return
      end if
      b_data = matmul(lambda, u_spherical(1:rank, :))
      status = 0
   end subroutine reduced_rank_transform

   pure subroutine covariance_term_nll(u, theta, code, nll, corr, sd, status, times, dist, factor_loadings)
      real(dp), intent(in) :: u(:, :) !! Random-effect matrix, rows form one covariance block and columns are block replicates.
      real(dp), intent(in) :: theta(:) !! Unconstrained covariance parameters in the ordering used by glmmTMB.
      integer, intent(in) :: code !! glmmTMB covariance-structure code from glmmtmb_codes.
      real(dp), intent(out) :: nll !! Negative log density of all columns of u under the requested covariance structure.
      real(dp), allocatable, intent(out) :: corr(:, :) !! Implied correlation matrix, or an empty matrix when not applicable.
      real(dp), allocatable, intent(out) :: sd(:) !! Implied marginal standard deviations, or an empty vector when not applicable.
      integer, intent(out) :: status !! Zero on success, nonzero for incompatible dimensions, parameters, or covariance matrices.
      real(dp), intent(in), optional :: times(:) !! Sorted time coordinates required for the Ornstein-Uhlenbeck structure.
      real(dp), intent(in), optional :: dist(:, :) !! Pairwise distance matrix required by spatial covariance structures.
      real(dp), allocatable, intent(out), optional :: factor_loadings(:, :) !! Reduced-rank loading matrix.
      real(dp), allocatable :: work_corr(:, :), work_sd(:), lambda(:, :)
      real(dp) :: a, rho, phi, raw, dval
      integer :: i, j, k, n, info, rank
      n = size(u, 1)
      nll = 0.0_dp
      status = 0
      allocate(corr(0, 0), sd(0))
      if (present(factor_loadings)) allocate(factor_loadings(0, 0))
      if (n <= 0) then
         status = 1
         return
      end if

      select case (code)
      case (diag_covstruct)
         if (size(theta) /= n) then
            status = 2
            return
         end if
         deallocate(sd)
         allocate(sd(n))
         sd = exp(theta)
         do j = 1, size(u, 2)
            do i = 1, n
               nll = nll - dnorm(u(i, j), 0.0_dp, sd(i), .true.)
            end do
         end do

      case (homdiag_covstruct)
         if (size(theta) /= 1) then
            status = 2
            return
         end if
         deallocate(sd)
         allocate(sd(n))
         sd = exp(theta(1))
         do j = 1, size(u, 2)
            do i = 1, n
               nll = nll - dnorm(u(i, j), 0.0_dp, sd(i), .true.)
            end do
         end do

      case (us_covstruct, equalto_covstruct, propto_covstruct)
         k = n * (n - 1) / 2
         if (code == propto_covstruct) then
            if (size(theta) /= n + k + 1) then
               status = 2
               return
            end if
         else
            if (size(theta) /= n + k) then
               status = 2
               return
            end if
         end if
         allocate(work_corr(n, n), work_sd(n))
         call unstructured_corr(theta(n + 1:n + k), work_corr, info)
         if (info /= 0) then
            status = 3
            return
         end if
         if (code == propto_covstruct) then
            work_sd = exp(theta(1:n) + 0.5_dp * theta(size(theta)))
         else
            work_sd = exp(theta(1:n))
         end if
         do j = 1, size(u, 2)
            nll = nll + scaled_mvn_nll(u(:, j), work_corr, work_sd)
         end do
         call move_alloc(work_corr, corr)
         call move_alloc(work_sd, sd)

      case (cs_covstruct, homcs_covstruct)
         if (n < 2) then
            status = 2
            return
         end if
         if (code == cs_covstruct) then
            if (size(theta) /= n + 1) then
               status = 2
               return
            end if
            allocate(work_sd(n))
            work_sd = exp(theta(1:n))
            raw = theta(n + 1)
         else
            if (size(theta) /= 2) then
               status = 2
               return
            end if
            allocate(work_sd(n))
            work_sd = exp(theta(1))
            raw = theta(2)
         end if
         a = 1.0_dp / real(n - 1, dp)
         rho = invlogit(raw) * (1.0_dp + a) - a
         allocate(work_corr(n, n))
         work_corr = rho
         do i = 1, n
            work_corr(i, i) = 1.0_dp
         end do
         do j = 1, size(u, 2)
            nll = nll + scaled_mvn_nll(u(:, j), work_corr, work_sd)
         end do
         call move_alloc(work_corr, corr)
         call move_alloc(work_sd, sd)

      case (toep_covstruct, homtoep_covstruct)
         if (code == toep_covstruct) then
            if (size(theta) /= 2 * n - 1) then
               status = 2
               return
            end if
            allocate(work_sd(n))
            work_sd = exp(theta(1:n))
            k = n
         else
            if (size(theta) /= n) then
               status = 2
               return
            end if
            allocate(work_sd(n))
            work_sd = exp(theta(1))
            k = 1
         end if
         allocate(work_corr(n, n))
         do j = 1, n
            do i = 1, n
               if (i == j) then
                  work_corr(i, j) = 1.0_dp
               else
                  raw = theta(k + abs(i - j))
                  work_corr(i, j) = raw / sqrt(1.0_dp + raw * raw)
               end if
            end do
         end do
         do j = 1, size(u, 2)
            nll = nll + scaled_mvn_nll(u(:, j), work_corr, work_sd)
         end do
         call move_alloc(work_corr, corr)
         call move_alloc(work_sd, sd)

      case (ar1_covstruct, hetar1_covstruct)
         if (code == ar1_covstruct) then
            if (size(theta) /= 2) then
               status = 2
               return
            end if
            allocate(work_sd(n))
            work_sd = exp(theta(1))
            raw = theta(2)
         else
            if (size(theta) /= n + 1) then
               status = 2
               return
            end if
            allocate(work_sd(n))
            work_sd = exp(theta(1:n))
            raw = theta(n + 1)
         end if
         phi = raw / sqrt(1.0_dp + raw * raw)
         allocate(work_corr(n, n))
         do j = 1, n
            do i = 1, n
               work_corr(i, j) = phi**abs(i - j)
            end do
         end do
         do j = 1, size(u, 2)
            nll = nll + scaled_mvn_nll(u(:, j), work_corr, work_sd)
         end do
         call move_alloc(work_corr, corr)
         call move_alloc(work_sd, sd)

      case (ou_covstruct)
         if (.not. present(times) .or. size(theta) /= 2 .or. size(times) /= n) then
            status = 2
            return
         end if
         allocate(work_sd(n), work_corr(n, n))
         work_sd = exp(theta(1))
         do j = 1, n
            do i = 1, n
               work_corr(i, j) = exp(-exp(theta(2)) * abs(times(i) - times(j)))
            end do
         end do
         do j = 1, size(u, 2)
            nll = nll + scaled_mvn_nll(u(:, j), work_corr, work_sd)
         end do
         call move_alloc(work_corr, corr)
         call move_alloc(work_sd, sd)

      case (exp_covstruct, gau_covstruct, mat_covstruct)
         if (.not. present(dist)) then
            status = 2
            return
         end if
         if (size(dist, 1) /= n .or. size(dist, 2) /= n) then
            status = 2
            return
         end if
         if ((code == mat_covstruct .and. size(theta) /= 3) .or. &
             (code /= mat_covstruct .and. size(theta) /= 2)) then
            status = 2
            return
         end if
         allocate(work_sd(n), work_corr(n, n))
         work_sd = exp(theta(1))
         do j = 1, n
            do i = 1, n
               dval = dist(i, j)
               if (i == j) then
                  work_corr(i, j) = 1.0_dp
               else if (code == exp_covstruct) then
                  work_corr(i, j) = exp(-dval * exp(-theta(2)))
               else if (code == gau_covstruct) then
                  work_corr(i, j) = exp(-dval * dval * exp(-2.0_dp * theta(2)))
               else
                  work_corr(i, j) = matern_corr(dval, exp(theta(2)), exp(theta(3)))
               end if
            end do
         end do
         do j = 1, size(u, 2)
            nll = nll + scaled_mvn_nll(u(:, j), work_corr, work_sd)
         end do
         call move_alloc(work_corr, corr)
         call move_alloc(work_sd, sd)

      case (rr_covstruct)
         call reduced_rank_loadings(theta, n, lambda, rank, info)
         if (info /= 0) then
            status = 2
            return
         end if
         do j = 1, size(u, 2)
            nll = nll + sum(n01_nll(u(:, j)))
         end do
         if (present(factor_loadings)) then
            deallocate(factor_loadings)
            allocate(factor_loadings(size(lambda, 1), size(lambda, 2)))
            factor_loadings = lambda
         end if
         allocate(work_corr(n, n), work_sd(n))
         work_corr = matmul(lambda, transpose(lambda))
         do i = 1, n
            work_sd(i) = sqrt(work_corr(i, i))
         end do
         do j = 1, n
            do i = 1, n
               if (work_sd(i) > 0.0_dp .and. work_sd(j) > 0.0_dp) then
                  work_corr(i, j) = work_corr(i, j) / (work_sd(i) * work_sd(j))
               end if
            end do
         end do
         call move_alloc(work_corr, corr)
         call move_alloc(work_sd, sd)

      case default
         status = 9
      end select

      if (status == 0 .and. ieee_is_nan(nll)) status = 4
   end subroutine covariance_term_nll
end module glmmtmb_covariance
