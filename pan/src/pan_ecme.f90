! SPDX-License-Identifier: GPL-3.0-only
! Derived from computational code in R package pan 2.0.
! Upstream authorship/maintenance: Joseph L. Schafer and Jing Hua Zhao.
! Maximum-likelihood mixed-model kernel corresponding to ecme() in R package pan.
module pan_ecme
   use pan_kinds, only : dp
   use pan_linalg, only : is_spd, spd_inverse, spd_logdet, spd_solve_vec, symmetrize
   use pan_types, only : PAN_ERR_ARGUMENT, PAN_ERR_DIMENSION, PAN_ERR_LINALG, PAN_OK, ecme_result
   implicit none
   private

   real(dp), parameter :: log_two_pi = log(2.0_dp * acos(-1.0_dp))

   public :: ecme_fit

contains

   subroutine ecme_fit(y, subj, occ, pred, xcol, result, zcol, vmax, beta_start, psi_start, sigma2_start, maxits, tol)
      real(dp), intent(in) :: y(:) !! Complete univariate response vector stacked by subject or cluster.
      integer, intent(in) :: subj(:) !! Sorted subject/cluster labels, one for each response.
      integer, intent(in) :: occ(:) !! One-based occasion labels selecting rows and columns from vmax.
      real(dp), intent(in) :: pred(:, :) !! Complete predictor matrix with one row per response.
      integer, intent(in) :: xcol(:) !! One-based predictor columns used for fixed effects.
      type(ecme_result), intent(out) :: result !! ML estimates, covariance, log likelihood, and empirical Bayes effects.
      integer, intent(in), optional :: zcol(:) !! One-based random-effect predictor columns; absent means no random effects.
      real(dp), intent(in), optional :: vmax(:, :) !! Known maximum residual correlation matrix; identity if absent.
      real(dp), intent(in), optional :: beta_start(:) !! Optional starting fixed-effect coefficients.
      real(dp), intent(in), optional :: psi_start(:, :) !! Optional starting random-effect covariance matrix.
      real(dp), intent(in), optional :: sigma2_start !! Optional positive starting residual scale.
      integer, intent(in), optional :: maxits !! Maximum EM/ECME-target iterations; default 1000.
      real(dp), intent(in), optional :: tol !! Relative convergence tolerance; default 1.0e-4.

      integer :: max_iterations
      integer :: q
      integer :: stat
      real(dp) :: tolerance
      integer, allocatable :: first(:)
      integer, allocatable :: last(:)

      call clear_ecme_result(result)

      max_iterations = 1000
      if (present(maxits)) max_iterations = maxits
      tolerance = 1.0e-4_dp
      if (present(tol)) tolerance = tol

      call validate_ecme(y, subj, occ, pred, xcol, vmax, max_iterations, tolerance, stat)
      if (stat /= PAN_OK) then
         call set_ecme_error(result, stat, "invalid ecme data, design, occasion, or control arguments")
         return
      end if

      call cluster_bounds_ecme(subj, first, last, stat)
      if (stat /= PAN_OK) then
         call set_ecme_error(result, stat, "subject labels must be sorted into contiguous blocks")
         return
      end if

      q = 0
      if (present(zcol)) q = size(zcol)

      if (q == 0) then
         call fit_gls_no_random(y, occ, pred, xcol, first, last, vmax, result, stat)
      else
         if (any(zcol < 1) .or. any(zcol > size(pred, 2))) then
            call set_ecme_error(result, PAN_ERR_ARGUMENT, "zcol contains an invalid predictor column")
            return
         end if
         call fit_em_random(y, occ, pred, xcol, zcol, first, last, vmax, beta_start, psi_start, &
            sigma2_start, max_iterations, tolerance, result, stat)
      end if

      if (stat /= PAN_OK .and. result%status == PAN_OK) then
         call set_ecme_error(result, stat, "numerical failure in ecme mixed-model fitting")
      end if
   end subroutine ecme_fit

   subroutine fit_gls_no_random(y, occ, pred, xcol, first, last, vmax, result, status)
      real(dp), intent(in) :: y(:) !! Complete response vector.
      integer, intent(in) :: occ(:) !! One-based occasion labels.
      real(dp), intent(in) :: pred(:, :) !! Predictor matrix.
      integer, intent(in) :: xcol(:) !! One-based fixed-effect predictor columns.
      integer, intent(in) :: first(:) !! First response row for each cluster.
      integer, intent(in) :: last(:) !! Last response row for each cluster.
      real(dp), intent(in), optional :: vmax(:, :) !! Known maximum residual correlation matrix; identity if absent.
      type(ecme_result), intent(inout) :: result !! Result object receiving exact GLS estimates.
      integer, intent(out) :: status !! PAN_OK on success or a numerical error code.

      integer :: ni
      integer :: s
      integer :: stat
      real(dp) :: qform
      real(dp), allocatable :: beta(:)
      real(dp), allocatable :: delta(:)
      real(dp), allocatable :: info(:, :)
      real(dp), allocatable :: info_inv(:, :)
      real(dp), allocatable :: rhs(:)
      real(dp), allocatable :: vi(:, :)
      real(dp), allocatable :: vi_inv(:, :)
      real(dp), allocatable :: x(:, :)
      real(dp), allocatable :: yi(:)

      allocate(info(size(xcol), size(xcol)), info_inv(size(xcol), size(xcol)))
      allocate(rhs(size(xcol)), beta(size(xcol)))
      info = 0.0_dp
      rhs = 0.0_dp

      do s = 1, size(first)
         ni = last(s) - first(s) + 1
         allocate(vi(ni, ni), vi_inv(ni, ni), x(ni, size(xcol)), yi(ni))
         call cluster_vi(occ(first(s):last(s)), vmax, vi)
         call spd_inverse(vi, vi_inv, stat)
         if (stat /= 0) then
            status = PAN_ERR_LINALG
            return
         end if

         x = pred(first(s):last(s), xcol)
         yi = y(first(s):last(s))
         info = info + matmul(transpose(x), matmul(vi_inv, x))
         rhs = rhs + matmul(transpose(x), matmul(vi_inv, yi))
         deallocate(vi, vi_inv, x, yi)
      end do

      call spd_inverse(info, info_inv, stat)
      if (stat /= 0) then
         status = PAN_ERR_LINALG
         return
      end if
      beta = matmul(info_inv, rhs)

      qform = 0.0_dp
      do s = 1, size(first)
         ni = last(s) - first(s) + 1
         allocate(vi(ni, ni), vi_inv(ni, ni), delta(ni))
         call cluster_vi(occ(first(s):last(s)), vmax, vi)
         call spd_inverse(vi, vi_inv, stat)
         if (stat /= 0) then
            status = PAN_ERR_LINALG
            return
         end if
         delta = y(first(s):last(s)) - matmul(pred(first(s):last(s), xcol), beta)
         qform = qform + dot_product(delta, matmul(vi_inv, delta))
         deallocate(vi, vi_inv, delta)
      end do

      result%beta = beta
      result%sigma2 = qform / real(size(y), dp)
      allocate(result%psi(0, 0), result%bhat(0, size(first)), result%cov_b(0, 0, size(first)))
      result%cov_beta = result%sigma2 * info_inv
      allocate(result%loglik(1))
      result%loglik(1) = marginal_loglik_no_random(y, occ, pred, xcol, first, last, vmax, beta, result%sigma2, stat)
      result%iter = 1
      result%converged = .true.
      result%status = PAN_OK
      result%message = "ok"
      status = PAN_OK
   end subroutine fit_gls_no_random

   subroutine fit_em_random(y, occ, pred, xcol, zcol, first, last, vmax, beta_start, psi_start, sigma2_start, &
         maxits, tol, result, status)
      real(dp), intent(in) :: y(:) !! Complete response vector.
      integer, intent(in) :: occ(:) !! One-based occasion labels.
      real(dp), intent(in) :: pred(:, :) !! Predictor matrix.
      integer, intent(in) :: xcol(:) !! One-based fixed-effect predictor columns.
      integer, intent(in) :: zcol(:) !! One-based random-effect predictor columns.
      integer, intent(in) :: first(:) !! First response row for each cluster.
      integer, intent(in) :: last(:) !! Last response row for each cluster.
      real(dp), intent(in), optional :: vmax(:, :) !! Known maximum residual correlation matrix; identity if absent.
      real(dp), intent(in), optional :: beta_start(:) !! Optional starting fixed-effect coefficients.
      real(dp), intent(in), optional :: psi_start(:, :) !! Optional starting random-effect covariance.
      real(dp), intent(in), optional :: sigma2_start !! Optional starting positive residual scale.
      integer, intent(in) :: maxits !! Maximum number of EM iterations.
      real(dp), intent(in) :: tol !! Relative parameter convergence tolerance.
      type(ecme_result), intent(inout) :: result !! Result object receiving fitted parameters and diagnostics.
      integer, intent(out) :: status !! PAN_OK on success or a numerical error code.

      integer :: it
      integer :: m
      integer :: p
      integer :: q
      integer :: stat
      real(dp) :: diff
      real(dp) :: old_sigma2
      real(dp), allocatable :: beta(:)
      real(dp), allocatable :: bhat(:, :)
      real(dp), allocatable :: covb(:, :, :)
      real(dp), allocatable :: ll(:)
      real(dp), allocatable :: old_beta(:)
      real(dp), allocatable :: old_psi(:, :)
      real(dp), allocatable :: psi(:, :)
      real(dp) :: sigma2

      p = size(xcol)
      q = size(zcol)
      m = size(first)
      allocate(beta(p), psi(q, q), bhat(q, m), covb(q, q, m))
      allocate(old_beta(p), old_psi(q, q), ll(maxits))

      call starting_values_random(y, occ, pred, xcol, zcol, first, last, vmax, beta, psi, sigma2, status)
      if (status /= PAN_OK) return

      if (present(beta_start)) then
         if (size(beta_start) /= p) then
            status = PAN_ERR_DIMENSION
            return
         end if
         beta = beta_start
      end if
      if (present(psi_start)) then
         if (size(psi_start, 1) /= q .or. size(psi_start, 2) /= q .or. .not. is_spd(psi_start)) then
            status = PAN_ERR_ARGUMENT
            return
         end if
         psi = psi_start
      end if
      if (present(sigma2_start)) then
         if (sigma2_start <= 0.0_dp) then
            status = PAN_ERR_ARGUMENT
            return
         end if
         sigma2 = sigma2_start
      end if

      result%converged = .false.
      do it = 1, maxits
         old_beta = beta
         old_psi = psi
         old_sigma2 = sigma2

         call em_iteration(y, occ, pred, xcol, zcol, first, last, vmax, beta, psi, sigma2, bhat, covb, status)
         if (status /= PAN_OK) return

         ll(it) = marginal_loglik_random(y, occ, pred, xcol, zcol, first, last, vmax, beta, psi, sigma2, stat)
         if (stat /= PAN_OK) then
            status = stat
            return
         end if

         diff = max_relative_change(beta, old_beta)
         diff = max(diff, max_relative_change_matrix(psi, old_psi))
         diff = max(diff, abs(sigma2 - old_sigma2) / max(abs(old_sigma2), 1.0e-12_dp))

         if (diff <= tol) then
            result%converged = .true.
            exit
         end if
      end do

      result%iter = min(it, maxits)
      result%beta = beta
      result%sigma2 = sigma2
      result%psi = psi
      result%bhat = bhat
      result%cov_b = covb
      result%loglik = ll(1:result%iter)

      call final_cov_beta(y, occ, pred, xcol, zcol, first, last, vmax, psi, sigma2, result%cov_beta, status)
      if (status /= PAN_OK) return

      call e_step_random(y, occ, pred, xcol, zcol, first, last, vmax, beta, psi, sigma2, result%bhat, result%cov_b, status)
      if (status /= PAN_OK) return

      result%status = PAN_OK
      result%message = "ok"
      status = PAN_OK
   end subroutine fit_em_random

   subroutine em_iteration(y, occ, pred, xcol, zcol, first, last, vmax, beta, psi, sigma2, bhat, covb, status)
      real(dp), intent(in) :: y(:) !! Complete response vector.
      integer, intent(in) :: occ(:) !! One-based occasion labels.
      real(dp), intent(in) :: pred(:, :) !! Predictor matrix.
      integer, intent(in) :: xcol(:) !! One-based fixed-effect predictor columns.
      integer, intent(in) :: zcol(:) !! One-based random-effect predictor columns.
      integer, intent(in) :: first(:) !! First response row for each cluster.
      integer, intent(in) :: last(:) !! Last response row for each cluster.
      real(dp), intent(in), optional :: vmax(:, :) !! Known maximum residual correlation matrix; identity if absent.
      real(dp), intent(inout) :: beta(:) !! Fixed-effect coefficients updated by the EM cycle.
      real(dp), intent(inout) :: psi(:, :) !! Random-effect covariance updated by the EM cycle.
      real(dp), intent(inout) :: sigma2 !! Positive residual scale updated by the EM cycle.
      real(dp), intent(out) :: bhat(:, :) !! Posterior random-effect means under pre-update parameters.
      real(dp), intent(out) :: covb(:, :, :) !! Posterior random-effect covariances under pre-update parameters.
      integer, intent(out) :: status !! PAN_OK on success or a numerical error code.

      integer :: ni
      integer :: s
      integer :: stat
      real(dp) :: residual_ss
      real(dp), allocatable :: info(:, :)
      real(dp), allocatable :: info_inv(:, :)
      real(dp), allocatable :: rhs(:)
      real(dp), allocatable :: vi(:, :)
      real(dp), allocatable :: vi_inv(:, :)
      real(dp), allocatable :: x(:, :)
      real(dp), allocatable :: yadj(:)
      real(dp), allocatable :: z(:, :)
      real(dp), allocatable :: zcz(:, :)

      call e_step_random(y, occ, pred, xcol, zcol, first, last, vmax, beta, psi, sigma2, bhat, covb, status)
      if (status /= PAN_OK) return

      allocate(info(size(xcol), size(xcol)), info_inv(size(xcol), size(xcol)), rhs(size(xcol)))
      info = 0.0_dp
      rhs = 0.0_dp

      do s = 1, size(first)
         ni = last(s) - first(s) + 1
         allocate(vi(ni, ni), vi_inv(ni, ni), x(ni, size(xcol)), z(ni, size(zcol)), yadj(ni))
         call cluster_vi(occ(first(s):last(s)), vmax, vi)
         call spd_inverse(vi, vi_inv, stat)
         if (stat /= 0) then
            status = PAN_ERR_LINALG
            return
         end if

         x = pred(first(s):last(s), xcol)
         z = pred(first(s):last(s), zcol)
         yadj = y(first(s):last(s)) - matmul(z, bhat(:, s))
         info = info + matmul(transpose(x), matmul(vi_inv, x))
         rhs = rhs + matmul(transpose(x), matmul(vi_inv, yadj))
         deallocate(vi, vi_inv, x, z, yadj)
      end do

      call spd_inverse(info, info_inv, stat)
      if (stat /= 0) then
         status = PAN_ERR_LINALG
         return
      end if
      beta = matmul(info_inv, rhs)

      psi = 0.0_dp
      do s = 1, size(first)
         psi = psi + covb(:, :, s) + outer_product_ecme(bhat(:, s), bhat(:, s))
      end do
      psi = psi / real(size(first), dp)
      call symmetrize(psi)

      residual_ss = 0.0_dp
      do s = 1, size(first)
         ni = last(s) - first(s) + 1
         allocate(vi(ni, ni), vi_inv(ni, ni), z(ni, size(zcol)), yadj(ni), zcz(ni, ni))
         call cluster_vi(occ(first(s):last(s)), vmax, vi)
         call spd_inverse(vi, vi_inv, stat)
         if (stat /= 0) then
            status = PAN_ERR_LINALG
            return
         end if

         z = pred(first(s):last(s), zcol)
         yadj = y(first(s):last(s)) - matmul(pred(first(s):last(s), xcol), beta) - matmul(z, bhat(:, s))
         zcz = matmul(matmul(z, covb(:, :, s)), transpose(z))
         residual_ss = residual_ss + dot_product(yadj, matmul(vi_inv, yadj)) + trace_product(vi_inv, zcz)
         deallocate(vi, vi_inv, z, yadj, zcz)
      end do
      sigma2 = residual_ss / real(size(y), dp)

      if (sigma2 <= 0.0_dp .or. .not. is_spd(psi)) then
         status = PAN_ERR_LINALG
      else
         status = PAN_OK
      end if
   end subroutine em_iteration

   subroutine e_step_random(y, occ, pred, xcol, zcol, first, last, vmax, beta, psi, sigma2, bhat, covb, status)
      real(dp), intent(in) :: y(:) !! Complete response vector.
      integer, intent(in) :: occ(:) !! One-based occasion labels.
      real(dp), intent(in) :: pred(:, :) !! Predictor matrix.
      integer, intent(in) :: xcol(:) !! One-based fixed-effect predictor columns.
      integer, intent(in) :: zcol(:) !! One-based random-effect predictor columns.
      integer, intent(in) :: first(:) !! First response row for each cluster.
      integer, intent(in) :: last(:) !! Last response row for each cluster.
      real(dp), intent(in), optional :: vmax(:, :) !! Known maximum residual correlation matrix; identity if absent.
      real(dp), intent(in) :: beta(:) !! Current fixed-effect coefficients.
      real(dp), intent(in) :: psi(:, :) !! Current random-effect covariance matrix.
      real(dp), intent(in) :: sigma2 !! Current positive residual scale.
      real(dp), intent(out) :: bhat(:, :) !! Conditional random-effect means by cluster.
      real(dp), intent(out) :: covb(:, :, :) !! Conditional random-effect covariance matrices by cluster.
      integer, intent(out) :: status !! PAN_OK on success or a numerical error code.

      integer :: ni
      integer :: s
      integer :: stat
      real(dp), allocatable :: delta(:)
      real(dp), allocatable :: temp(:, :)
      real(dp), allocatable :: v(:, :)
      real(dp), allocatable :: vi(:, :)
      real(dp), allocatable :: vinv(:, :)
      real(dp), allocatable :: z(:, :)

      do s = 1, size(first)
         ni = last(s) - first(s) + 1
         allocate(vi(ni, ni), v(ni, ni), vinv(ni, ni), z(ni, size(zcol)), delta(ni))
         allocate(temp(size(zcol), ni))

         call cluster_vi(occ(first(s):last(s)), vmax, vi)
         z = pred(first(s):last(s), zcol)
         v = sigma2 * vi + matmul(matmul(z, psi), transpose(z))
         call symmetrize(v)
         call spd_inverse(v, vinv, stat)
         if (stat /= 0) then
            status = PAN_ERR_LINALG
            return
         end if

         delta = y(first(s):last(s)) - matmul(pred(first(s):last(s), xcol), beta)
         temp = matmul(matmul(psi, transpose(z)), vinv)
         bhat(:, s) = matmul(temp, delta)
         covb(:, :, s) = psi - matmul(matmul(temp, z), psi)
         call symmetrize(covb(:, :, s))

         if (.not. is_spd(covb(:, :, s))) call stabilize_cov_ecme(covb(:, :, s))
         if (.not. is_spd(covb(:, :, s))) then
            status = PAN_ERR_LINALG
            return
         end if

         deallocate(vi, v, vinv, z, delta, temp)
      end do

      status = PAN_OK
   end subroutine e_step_random

   subroutine starting_values_random(y, occ, pred, xcol, zcol, first, last, vmax, beta, psi, sigma2, status)
      real(dp), intent(in) :: y(:) !! Complete response vector.
      integer, intent(in) :: occ(:) !! One-based occasion labels.
      real(dp), intent(in) :: pred(:, :) !! Predictor matrix.
      integer, intent(in) :: xcol(:) !! One-based fixed-effect predictor columns.
      integer, intent(in) :: zcol(:) !! One-based random-effect predictor columns.
      integer, intent(in) :: first(:) !! First response row for each cluster.
      integer, intent(in) :: last(:) !! Last response row for each cluster.
      real(dp), intent(in), optional :: vmax(:, :) !! Known maximum residual correlation matrix; identity if absent.
      real(dp), intent(out) :: beta(:) !! Initial fixed-effect coefficients.
      real(dp), intent(out) :: psi(:, :) !! Initial random-effect covariance matrix.
      real(dp), intent(out) :: sigma2 !! Initial positive residual scale.
      integer, intent(out) :: status !! PAN_OK on success or a numerical error code.

      type(ecme_result) :: gls
      integer :: i
      integer :: stat

      call fit_gls_no_random(y, occ, pred, xcol, first, last, vmax, gls, stat)
      if (stat /= PAN_OK) then
         status = stat
         return
      end if

      beta = gls%beta
      sigma2 = max(gls%sigma2, 1.0e-8_dp)
      psi = 0.0_dp
      do i = 1, size(zcol)
         psi(i, i) = max(0.25_dp * sigma2, 1.0e-6_dp)
      end do
      status = PAN_OK
   end subroutine starting_values_random

   subroutine final_cov_beta(y, occ, pred, xcol, zcol, first, last, vmax, psi, sigma2, cov_beta, status)
      real(dp), intent(in) :: y(:) !! Response vector; used only for dimension consistency.
      integer, intent(in) :: occ(:) !! One-based occasion labels.
      real(dp), intent(in) :: pred(:, :) !! Predictor matrix.
      integer, intent(in) :: xcol(:) !! One-based fixed-effect predictor columns.
      integer, intent(in) :: zcol(:) !! One-based random-effect predictor columns.
      integer, intent(in) :: first(:) !! First response row for each cluster.
      integer, intent(in) :: last(:) !! Last response row for each cluster.
      real(dp), intent(in), optional :: vmax(:, :) !! Known maximum residual correlation matrix; identity if absent.
      real(dp), intent(in) :: psi(:, :) !! Final random-effect covariance matrix.
      real(dp), intent(in) :: sigma2 !! Final residual scale.
      real(dp), allocatable, intent(out) :: cov_beta(:, :) !! Inverse marginal fixed-effect information matrix.
      integer, intent(out) :: status !! PAN_OK on success or PAN_ERR_LINALG if information is singular.

      integer :: ni
      integer :: s
      integer :: stat
      real(dp), allocatable :: info(:, :)
      real(dp), allocatable :: v(:, :)
      real(dp), allocatable :: vi(:, :)
      real(dp), allocatable :: vinv(:, :)
      real(dp), allocatable :: x(:, :)
      real(dp), allocatable :: z(:, :)

      if (size(y) /= size(pred, 1)) then
         status = PAN_ERR_DIMENSION
         return
      end if

      allocate(info(size(xcol), size(xcol)))
      info = 0.0_dp

      do s = 1, size(first)
         ni = last(s) - first(s) + 1
         allocate(vi(ni, ni), v(ni, ni), vinv(ni, ni), x(ni, size(xcol)), z(ni, size(zcol)))
         call cluster_vi(occ(first(s):last(s)), vmax, vi)
         x = pred(first(s):last(s), xcol)
         z = pred(first(s):last(s), zcol)
         v = sigma2 * vi + matmul(matmul(z, psi), transpose(z))
         call symmetrize(v)
         call spd_inverse(v, vinv, stat)
         if (stat /= 0) then
            status = PAN_ERR_LINALG
            return
         end if
         info = info + matmul(transpose(x), matmul(vinv, x))
         deallocate(vi, v, vinv, x, z)
      end do

      allocate(cov_beta(size(xcol), size(xcol)))
      call spd_inverse(info, cov_beta, stat)
      if (stat /= 0) then
         status = PAN_ERR_LINALG
      else
         status = PAN_OK
      end if
   end subroutine final_cov_beta

   function marginal_loglik_random(y, occ, pred, xcol, zcol, first, last, vmax, beta, psi, sigma2, status) result(ll)
      real(dp), intent(in) :: y(:) !! Complete response vector.
      integer, intent(in) :: occ(:) !! One-based occasion labels.
      real(dp), intent(in) :: pred(:, :) !! Predictor matrix.
      integer, intent(in) :: xcol(:) !! One-based fixed-effect predictor columns.
      integer, intent(in) :: zcol(:) !! One-based random-effect predictor columns.
      integer, intent(in) :: first(:) !! First response row for each cluster.
      integer, intent(in) :: last(:) !! Last response row for each cluster.
      real(dp), intent(in), optional :: vmax(:, :) !! Known maximum residual correlation matrix; identity if absent.
      real(dp), intent(in) :: beta(:) !! Fixed-effect coefficients at which the likelihood is evaluated.
      real(dp), intent(in) :: psi(:, :) !! Random-effect covariance at which the likelihood is evaluated.
      real(dp), intent(in) :: sigma2 !! Positive residual scale at which the likelihood is evaluated.
      integer, intent(out) :: status !! PAN_OK on success or PAN_ERR_LINALG if a cluster covariance is invalid.
      real(dp) :: ll

      integer :: ni
      integer :: s
      integer :: stat
      real(dp) :: logdet
      real(dp), allocatable :: delta(:)
      real(dp), allocatable :: sol(:)
      real(dp), allocatable :: v(:, :)
      real(dp), allocatable :: vi(:, :)
      real(dp), allocatable :: z(:, :)

      ll = 0.0_dp
      do s = 1, size(first)
         ni = last(s) - first(s) + 1
         allocate(vi(ni, ni), v(ni, ni), z(ni, size(zcol)), delta(ni), sol(ni))
         call cluster_vi(occ(first(s):last(s)), vmax, vi)
         z = pred(first(s):last(s), zcol)
         v = sigma2 * vi + matmul(matmul(z, psi), transpose(z))
         call symmetrize(v)

         call spd_logdet(v, logdet, stat)
         if (stat /= 0) then
            status = PAN_ERR_LINALG
            return
         end if
         delta = y(first(s):last(s)) - matmul(pred(first(s):last(s), xcol), beta)
         call spd_solve_vec(v, delta, sol, stat)
         if (stat /= 0) then
            status = PAN_ERR_LINALG
            return
         end if

         ll = ll - 0.5_dp * (real(ni, dp) * log_two_pi + logdet + dot_product(delta, sol))
         deallocate(vi, v, z, delta, sol)
      end do
      status = PAN_OK
   end function marginal_loglik_random

   function marginal_loglik_no_random(y, occ, pred, xcol, first, last, vmax, beta, sigma2, status) result(ll)
      real(dp), intent(in) :: y(:) !! Complete response vector.
      integer, intent(in) :: occ(:) !! One-based occasion labels.
      real(dp), intent(in) :: pred(:, :) !! Predictor matrix.
      integer, intent(in) :: xcol(:) !! One-based fixed-effect predictor columns.
      integer, intent(in) :: first(:) !! First response row for each cluster.
      integer, intent(in) :: last(:) !! Last response row for each cluster.
      real(dp), intent(in), optional :: vmax(:, :) !! Known maximum residual correlation matrix; identity if absent.
      real(dp), intent(in) :: beta(:) !! Fixed-effect coefficients at which the likelihood is evaluated.
      real(dp), intent(in) :: sigma2 !! Positive residual scale at which the likelihood is evaluated.
      integer, intent(out) :: status !! PAN_OK on success or PAN_ERR_LINALG if a cluster covariance is invalid.
      real(dp) :: ll

      integer :: ni
      integer :: s
      integer :: stat
      real(dp) :: logdet
      real(dp), allocatable :: delta(:)
      real(dp), allocatable :: sol(:)
      real(dp), allocatable :: v(:, :)
      real(dp), allocatable :: vi(:, :)

      ll = 0.0_dp
      do s = 1, size(first)
         ni = last(s) - first(s) + 1
         allocate(vi(ni, ni), v(ni, ni), delta(ni), sol(ni))
         call cluster_vi(occ(first(s):last(s)), vmax, vi)
         v = sigma2 * vi

         call spd_logdet(v, logdet, stat)
         if (stat /= 0) then
            status = PAN_ERR_LINALG
            return
         end if
         delta = y(first(s):last(s)) - matmul(pred(first(s):last(s), xcol), beta)
         call spd_solve_vec(v, delta, sol, stat)
         if (stat /= 0) then
            status = PAN_ERR_LINALG
            return
         end if

         ll = ll - 0.5_dp * (real(ni, dp) * log_two_pi + logdet + dot_product(delta, sol))
         deallocate(vi, v, delta, sol)
      end do
      status = PAN_OK
   end function marginal_loglik_no_random

   pure subroutine cluster_vi(occ_cluster, vmax, vi)
      integer, intent(in) :: occ_cluster(:) !! One-based occasion labels for one cluster.
      real(dp), intent(in), optional :: vmax(:, :) !! Maximum known residual correlation matrix; identity if absent.
      real(dp), intent(out) :: vi(:, :) !! Cluster residual correlation submatrix.

      integer :: i
      integer :: j

      if (present(vmax)) then
         do j = 1, size(occ_cluster)
            do i = 1, size(occ_cluster)
               vi(i, j) = vmax(occ_cluster(i), occ_cluster(j))
            end do
         end do
      else
         vi = 0.0_dp
         do i = 1, size(occ_cluster)
            vi(i, i) = 1.0_dp
         end do
      end if
   end subroutine cluster_vi

   pure subroutine validate_ecme(y, subj, occ, pred, xcol, vmax, maxits, tol, status)
      real(dp), intent(in) :: y(:) !! Response vector whose dimensions are validated.
      integer, intent(in) :: subj(:) !! Subject labels whose dimensions and ordering are validated.
      integer, intent(in) :: occ(:) !! Occasion labels validated against vmax or their observed maximum.
      real(dp), intent(in) :: pred(:, :) !! Predictor matrix whose dimensions are validated.
      integer, intent(in) :: xcol(:) !! One-based fixed-effect predictor columns.
      real(dp), intent(in), optional :: vmax(:, :) !! Optional known residual correlation matrix.
      integer, intent(in) :: maxits !! Maximum iteration count; must be positive.
      real(dp), intent(in) :: tol !! Positive convergence tolerance.
      integer, intent(out) :: status !! PAN_OK if all structural argument checks pass.

      integer :: i
      integer :: nmax

      status = PAN_OK
      if (size(y) <= 0 .or. size(xcol) <= 0) then
         status = PAN_ERR_DIMENSION
         return
      end if
      if (size(subj) /= size(y) .or. size(occ) /= size(y) .or. size(pred, 1) /= size(y)) then
         status = PAN_ERR_DIMENSION
         return
      end if
      if (any(xcol < 1) .or. any(xcol > size(pred, 2)) .or. any(occ < 1)) then
         status = PAN_ERR_ARGUMENT
         return
      end if
      if (maxits <= 0 .or. tol <= 0.0_dp) then
         status = PAN_ERR_ARGUMENT
         return
      end if
      do i = 2, size(subj)
         if (subj(i) < subj(i - 1)) then
            status = PAN_ERR_ARGUMENT
            return
         end if
      end do

      nmax = maxval(occ)
      if (present(vmax)) then
         if (size(vmax, 1) < nmax .or. size(vmax, 2) < nmax) then
            status = PAN_ERR_DIMENSION
            return
         end if
         if (size(vmax, 1) /= size(vmax, 2) .or. .not. is_spd(vmax)) then
            status = PAN_ERR_ARGUMENT
            return
         end if
      end if
   end subroutine validate_ecme

   pure subroutine cluster_bounds_ecme(subj, first, last, status)
      integer, intent(in) :: subj(:) !! Sorted subject labels grouped contiguously.
      integer, allocatable, intent(out) :: first(:) !! First row of each subject block.
      integer, allocatable, intent(out) :: last(:) !! Last row of each subject block.
      integer, intent(out) :: status !! PAN_OK on success or PAN_ERR_ARGUMENT for decreasing labels.

      integer :: i
      integer :: m
      integer :: s

      do i = 2, size(subj)
         if (subj(i) < subj(i - 1)) then
            allocate(first(0), last(0))
            status = PAN_ERR_ARGUMENT
            return
         end if
      end do

      m = 1
      do i = 2, size(subj)
         if (subj(i) /= subj(i - 1)) m = m + 1
      end do
      allocate(first(m), last(m))

      first(1) = 1
      s = 1
      do i = 2, size(subj)
         if (subj(i) /= subj(i - 1)) then
            last(s) = i - 1
            s = s + 1
            first(s) = i
         end if
      end do
      last(m) = size(subj)
      status = PAN_OK
   end subroutine cluster_bounds_ecme

   pure function outer_product_ecme(x, y) result(a)
      real(dp), intent(in) :: x(:) !! Left vector in the dyadic product.
      real(dp), intent(in) :: y(:) !! Right vector in the dyadic product.
      real(dp) :: a(size(x), size(y))

      integer :: i
      integer :: j

      do j = 1, size(y)
         do i = 1, size(x)
            a(i, j) = x(i) * y(j)
         end do
      end do
   end function outer_product_ecme

   pure function trace_product(a, b) result(value)
      real(dp), intent(in) :: a(:, :) !! First square matrix in trace(A B).
      real(dp), intent(in) :: b(:, :) !! Second square matrix in trace(A B).
      real(dp) :: value

      integer :: i
      integer :: j

      value = 0.0_dp
      do i = 1, size(a, 1)
         do j = 1, size(a, 2)
            value = value + a(i, j) * b(j, i)
         end do
      end do
   end function trace_product

   pure function max_relative_change(x, old) result(value)
      real(dp), intent(in) :: x(:) !! New parameter vector.
      real(dp), intent(in) :: old(:) !! Previous parameter vector.
      real(dp) :: value

      value = maxval(abs(x - old) / max(abs(old), 1.0e-12_dp))
   end function max_relative_change

   pure function max_relative_change_matrix(x, old) result(value)
      real(dp), intent(in) :: x(:, :) !! New parameter matrix.
      real(dp), intent(in) :: old(:, :) !! Previous parameter matrix.
      real(dp) :: value

      value = maxval(abs(x - old) / max(abs(old), 1.0e-12_dp))
   end function max_relative_change_matrix

   pure subroutine stabilize_cov_ecme(a)
      real(dp), intent(inout) :: a(:, :) !! Nearly positive-definite covariance receiving small diagonal jitter.

      integer :: i
      integer :: k
      real(dp) :: scale

      scale = max(1.0_dp, maxval(abs(a)))
      do k = 1, 8
         do i = 1, size(a, 1)
            a(i, i) = a(i, i) + scale * 10.0_dp**(-12 + k)
         end do
         if (is_spd(a)) return
      end do
   end subroutine stabilize_cov_ecme

   subroutine clear_ecme_result(result)
      type(ecme_result), intent(out) :: result !! Result object reset before fitting.

      result%status = PAN_OK
      result%message = ""
      result%iter = 0
      result%converged = .false.
      result%sigma2 = 0.0_dp
   end subroutine clear_ecme_result

   subroutine set_ecme_error(result, status, message)
      type(ecme_result), intent(inout) :: result !! Result object receiving the error status and diagnostic.
      integer, intent(in) :: status !! Nonzero package status code.
      character(len=*), intent(in) :: message !! Human-readable diagnostic explaining the failure.

      result%status = status
      result%message = message
   end subroutine set_ecme_error

end module pan_ecme
