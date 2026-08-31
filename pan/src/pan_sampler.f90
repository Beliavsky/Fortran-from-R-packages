! SPDX-License-Identifier: GPL-3.0-only
! Derived from computational code in R package pan 2.0.
! Upstream authorship/maintenance: Joseph L. Schafer and Jing Hua Zhao.
! Gibbs samplers corresponding to pan() and pan.bd() from R package pan 2.0.
module pan_sampler
   use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
   use pan_kinds, only : dp
   use pan_linalg, only : invwishart_draw, is_spd, matrix_normal_draw, mvn_draw, spd_inverse, spd_solve_vec, symmetrize
   use pan_rng, only : rng_seed, rng_state
   use pan_types, only : PAN_ERR_ARGUMENT, PAN_ERR_DIMENSION, PAN_ERR_LINALG, PAN_OK, pan_bd_prior, pan_bd_result, &
      pan_bd_state, pan_prior, pan_result, pan_state
   implicit none
   private

   public :: pan_mcmc
   public :: pan_bd_mcmc

contains

   subroutine pan_mcmc(y, subj, pred, xcol, zcol, prior, seed, iter, result, start)
      real(dp), intent(in) :: y(:, :) !! Response matrix; IEEE NaNs mark missing response components.
      integer, intent(in) :: subj(:) !! Sorted subject/cluster labels, one label for each row of y.
      real(dp), intent(in) :: pred(:, :) !! Complete predictor matrix with the same number of rows as y.
      integer, intent(in) :: xcol(:) !! One-based predictor columns used for fixed effects.
      integer, intent(in) :: zcol(:) !! One-based predictor columns used for random effects.
      type(pan_prior), intent(in) :: prior !! Inverse-Wishart hyperparameters a/Binv and c/Dinv.
      integer, intent(in) :: seed !! Positive deterministic random-number seed.
      integer, intent(in) :: iter !! Number of Gibbs cycles to retain; must be positive.
      type(pan_result), intent(out) :: result !! Simulated parameter chains, final imputation, and restart state.
      type(pan_state), intent(in), optional :: start !! Optional complete restart state from an earlier run.

      integer :: m
      integer :: n
      integer :: p
      integer :: q
      integer :: r
      integer :: stat
      integer, allocatable :: first(:)
      integer, allocatable :: last(:)
      logical, allocatable :: observed(:, :)
      logical, allocatable :: usable(:)
      real(dp), allocatable :: b(:, :, :)
      real(dp), allocatable :: beta(:, :)
      real(dp), allocatable :: sigma(:, :)
      real(dp), allocatable :: psi(:, :)
      real(dp), allocatable :: ywork(:, :)
      type(rng_state) :: rng

      call clear_pan_result(result)
      n = size(y, 1)
      r = size(y, 2)
      p = size(xcol)
      q = size(zcol)

      call validate_common(y, subj, pred, xcol, zcol, iter, stat)
      if (stat /= PAN_OK) then
         call set_pan_error(result, stat, "invalid data, design-column, subject-order, or iteration dimensions")
         return
      end if

      if (.not. allocated(prior%binv) .or. .not. allocated(prior%dinv)) then
         call set_pan_error(result, PAN_ERR_ARGUMENT, "prior Binv and Dinv must be allocated")
         return
      end if
      if (size(prior%binv, 1) /= r .or. size(prior%binv, 2) /= r) then
         call set_pan_error(result, PAN_ERR_DIMENSION, "prior Binv must be r by r")
         return
      end if
      if (size(prior%dinv, 1) /= q * r .or. size(prior%dinv, 2) /= q * r) then
         call set_pan_error(result, PAN_ERR_DIMENSION, "prior Dinv must be (q*r) by (q*r)")
         return
      end if
      if (prior%a < real(r, dp) .or. prior%c < real(q * r, dp)) then
         call set_pan_error(result, PAN_ERR_ARGUMENT, "prior degrees of freedom are below pan's documented minima")
         return
      end if
      if (.not. is_spd(prior%binv) .or. .not. is_spd(prior%dinv)) then
         call set_pan_error(result, PAN_ERR_ARGUMENT, "prior scale matrices must be positive definite")
         return
      end if

      call cluster_bounds(subj, first, last, m, stat)
      if (stat /= PAN_OK) then
         call set_pan_error(result, stat, "subject labels must occur in contiguous sorted blocks")
         return
      end if

      allocate(observed(n, r), usable(n))
      observed = .not. ieee_is_nan(y)
      usable = any(observed, dim=2)

      if (count(usable) <= p) then
         call set_pan_error(result, PAN_ERR_ARGUMENT, "too few nonempty response rows for the fixed-effect design")
         return
      end if

      allocate(result%beta(p, r, iter), result%sigma(r, r, iter), result%psi(q * r, q * r, iter))
      allocate(result%y(n, r))
      allocate(beta(p, r), sigma(r, r), psi(q * r, q * r), ywork(n, r), b(q, r, m))
      b = 0.0_dp

      call rng_seed(rng, seed)

      if (present(start)) then
         call load_start_full(start, y, observed, p, q, r, beta, sigma, psi, ywork, stat)
         if (stat /= PAN_OK) then
            call set_pan_error(result, stat, "restart state has incompatible dimensions or invalid covariance matrices")
            return
         end if
      else
         call initialize_full(y, observed, usable, pred, xcol, prior, beta, sigma, psi, ywork, stat)
         if (stat /= PAN_OK) then
            call set_pan_error(result, stat, "could not construct valid initial values")
            return
         end if
      end if

      call run_full_chain(y, observed, usable, pred, xcol, zcol, first, last, prior, rng, &
         beta, sigma, psi, b, ywork, result, stat)
      if (stat /= PAN_OK) then
         call set_pan_error(result, stat, "numerical failure during Gibbs sampling")
         return
      end if

      result%y = ywork
      call copy_full_state(beta, sigma, psi, ywork, result%last)
      result%status = PAN_OK
      result%message = "ok"
   end subroutine pan_mcmc

   subroutine pan_bd_mcmc(y, subj, pred, xcol, zcol, prior, seed, iter, result, start)
      real(dp), intent(in) :: y(:, :) !! Response matrix; IEEE NaNs mark missing response components.
      integer, intent(in) :: subj(:) !! Sorted subject/cluster labels, one label for each row of y.
      real(dp), intent(in) :: pred(:, :) !! Complete predictor matrix with the same number of rows as y.
      integer, intent(in) :: xcol(:) !! One-based predictor columns used for fixed effects.
      integer, intent(in) :: zcol(:) !! One-based predictor columns used for random effects.
      type(pan_bd_prior), intent(in) :: prior !! Block-diagonal inverse-Wishart hyperparameters.
      integer, intent(in) :: seed !! Positive deterministic random-number seed.
      integer, intent(in) :: iter !! Number of Gibbs cycles to retain; must be positive.
      type(pan_bd_result), intent(out) :: result !! Simulated chains, final imputation, and block-diagonal restart state.
      type(pan_bd_state), intent(in), optional :: start !! Optional complete restart state from an earlier pan.bd run.

      integer :: j
      integer :: m
      integer :: n
      integer :: p
      integer :: q
      integer :: r
      integer :: stat
      integer, allocatable :: first(:)
      integer, allocatable :: last(:)
      logical, allocatable :: observed(:, :)
      logical, allocatable :: usable(:)
      real(dp), allocatable :: b(:, :, :)
      real(dp), allocatable :: beta(:, :)
      real(dp), allocatable :: sigma(:, :)
      real(dp), allocatable :: psi(:, :, :)
      real(dp), allocatable :: ywork(:, :)
      type(rng_state) :: rng

      call clear_pan_bd_result(result)
      n = size(y, 1)
      r = size(y, 2)
      p = size(xcol)
      q = size(zcol)

      call validate_common(y, subj, pred, xcol, zcol, iter, stat)
      if (stat /= PAN_OK) then
         call set_pan_bd_error(result, stat, "invalid data, design-column, subject-order, or iteration dimensions")
         return
      end if

      if (.not. allocated(prior%binv) .or. .not. allocated(prior%c) .or. .not. allocated(prior%dinv)) then
         call set_pan_bd_error(result, PAN_ERR_ARGUMENT, "block prior Binv, c, and Dinv must be allocated")
         return
      end if
      if (size(prior%binv, 1) /= r .or. size(prior%binv, 2) /= r) then
         call set_pan_bd_error(result, PAN_ERR_DIMENSION, "prior Binv must be r by r")
         return
      end if
      if (size(prior%c) /= r) then
         call set_pan_bd_error(result, PAN_ERR_DIMENSION, "block prior c must have length r")
         return
      end if
      if (size(prior%dinv, 1) /= q .or. size(prior%dinv, 2) /= q .or. size(prior%dinv, 3) /= r) then
         call set_pan_bd_error(result, PAN_ERR_DIMENSION, "block prior Dinv must be q by q by r")
         return
      end if
      if (prior%a < real(r, dp) .or. any(prior%c < real(q, dp))) then
         call set_pan_bd_error(result, PAN_ERR_ARGUMENT, "prior degrees of freedom are below pan.bd's documented minima")
         return
      end if
      if (.not. is_spd(prior%binv)) then
         call set_pan_bd_error(result, PAN_ERR_ARGUMENT, "prior Binv must be positive definite")
         return
      end if
      do j = 1, r
         if (.not. is_spd(prior%dinv(:, :, j))) then
            call set_pan_bd_error(result, PAN_ERR_ARGUMENT, "every prior Dinv block must be positive definite")
            return
         end if
      end do

      call cluster_bounds(subj, first, last, m, stat)
      if (stat /= PAN_OK) then
         call set_pan_bd_error(result, stat, "subject labels must occur in contiguous sorted blocks")
         return
      end if

      allocate(observed(n, r), usable(n))
      observed = .not. ieee_is_nan(y)
      usable = any(observed, dim=2)

      if (count(usable) <= p) then
         call set_pan_bd_error(result, PAN_ERR_ARGUMENT, "too few nonempty response rows for the fixed-effect design")
         return
      end if

      allocate(result%beta(p, r, iter), result%sigma(r, r, iter), result%psi(q, q, r, iter))
      allocate(result%y(n, r))
      allocate(beta(p, r), sigma(r, r), psi(q, q, r), ywork(n, r), b(q, r, m))
      b = 0.0_dp

      call rng_seed(rng, seed)

      if (present(start)) then
         call load_start_bd(start, y, observed, p, q, r, beta, sigma, psi, ywork, stat)
         if (stat /= PAN_OK) then
            call set_pan_bd_error(result, stat, "restart state has incompatible dimensions or invalid covariance matrices")
            return
         end if
      else
         call initialize_bd(y, observed, usable, pred, xcol, prior, beta, sigma, psi, ywork, stat)
         if (stat /= PAN_OK) then
            call set_pan_bd_error(result, stat, "could not construct valid initial values")
            return
         end if
      end if

      call run_bd_chain(y, observed, usable, pred, xcol, zcol, first, last, prior, rng, &
         beta, sigma, psi, b, ywork, result, stat)
      if (stat /= PAN_OK) then
         call set_pan_bd_error(result, stat, "numerical failure during block-diagonal Gibbs sampling")
         return
      end if

      result%y = ywork
      call copy_bd_state(beta, sigma, psi, ywork, result%last)
      result%status = PAN_OK
      result%message = "ok"
   end subroutine pan_bd_mcmc

   subroutine run_full_chain(y_original, observed, usable, pred, xcol, zcol, first, last, prior, rng, &
         beta, sigma, psi, b, ywork, result, status)
      real(dp), intent(in) :: y_original(:, :) !! Original response matrix used to preserve observed entries exactly.
      logical, intent(in) :: observed(:, :) !! True where the original response component was observed.
      logical, intent(in) :: usable(:) !! True for rows containing at least one originally observed response.
      real(dp), intent(in) :: pred(:, :) !! Complete predictor matrix.
      integer, intent(in) :: xcol(:) !! One-based fixed-effect predictor columns.
      integer, intent(in) :: zcol(:) !! One-based random-effect predictor columns.
      integer, intent(in) :: first(:) !! First row index of each contiguous subject block.
      integer, intent(in) :: last(:) !! Last row index of each contiguous subject block.
      type(pan_prior), intent(in) :: prior !! Full random-effect inverse-Wishart prior.
      type(rng_state), intent(inout) :: rng !! Mutable deterministic random-number state.
      real(dp), intent(inout) :: beta(:, :) !! Current fixed-effect coefficient matrix.
      real(dp), intent(inout) :: sigma(:, :) !! Current residual covariance matrix.
      real(dp), intent(inout) :: psi(:, :) !! Current full random-effect covariance matrix.
      real(dp), intent(inout) :: b(:, :, :) !! Current subject-specific random effects.
      real(dp), intent(inout) :: ywork(:, :) !! Current completed response matrix.
      type(pan_result), intent(inout) :: result !! Output chains populated one cycle at a time.
      integer, intent(out) :: status !! PAN_OK on success or a numerical error code.

      integer :: it
      real(dp), allocatable :: beta_mean(:, :)
      real(dp), allocatable :: xtx_inv(:, :)

      status = PAN_OK
      allocate(beta_mean(size(beta, 1), size(beta, 2)), xtx_inv(size(beta, 1), size(beta, 1)))

      do it = 1, size(result%beta, 3)
         call draw_random_effects(ywork, usable, pred, xcol, zcol, first, last, beta, sigma, psi, rng, b, status)
         if (status /= PAN_OK) return

         call draw_psi_full(b, prior, rng, psi, status)
         if (status /= PAN_OK) return

         call fixed_effect_mean(ywork, usable, pred, xcol, zcol, first, last, b, beta_mean, xtx_inv, status)
         if (status /= PAN_OK) return

         call draw_sigma_beta(ywork, usable, pred, xcol, zcol, first, last, b, prior%a, prior%binv, &
            beta_mean, xtx_inv, rng, sigma, beta, status)
         if (status /= PAN_OK) return

         call impute_missing(y_original, observed, pred, xcol, zcol, first, last, beta, b, sigma, rng, ywork, status)
         if (status /= PAN_OK) return

         result%beta(:, :, it) = beta
         result%sigma(:, :, it) = sigma
         result%psi(:, :, it) = psi
      end do
   end subroutine run_full_chain

   subroutine run_bd_chain(y_original, observed, usable, pred, xcol, zcol, first, last, prior, rng, &
         beta, sigma, psi, b, ywork, result, status)
      real(dp), intent(in) :: y_original(:, :) !! Original response matrix used to preserve observed entries exactly.
      logical, intent(in) :: observed(:, :) !! True where the original response component was observed.
      logical, intent(in) :: usable(:) !! True for rows containing at least one originally observed response.
      real(dp), intent(in) :: pred(:, :) !! Complete predictor matrix.
      integer, intent(in) :: xcol(:) !! One-based fixed-effect predictor columns.
      integer, intent(in) :: zcol(:) !! One-based random-effect predictor columns.
      integer, intent(in) :: first(:) !! First row index of each contiguous subject block.
      integer, intent(in) :: last(:) !! Last row index of each contiguous subject block.
      type(pan_bd_prior), intent(in) :: prior !! Block-diagonal random-effect inverse-Wishart prior.
      type(rng_state), intent(inout) :: rng !! Mutable deterministic random-number state.
      real(dp), intent(inout) :: beta(:, :) !! Current fixed-effect coefficient matrix.
      real(dp), intent(inout) :: sigma(:, :) !! Current residual covariance matrix.
      real(dp), intent(inout) :: psi(:, :, :) !! Current per-response random-effect covariance blocks.
      real(dp), intent(inout) :: b(:, :, :) !! Current subject-specific random effects.
      real(dp), intent(inout) :: ywork(:, :) !! Current completed response matrix.
      type(pan_bd_result), intent(inout) :: result !! Output chains populated one cycle at a time.
      integer, intent(out) :: status !! PAN_OK on success or a numerical error code.

      integer :: it
      real(dp), allocatable :: beta_mean(:, :)
      real(dp), allocatable :: psi_full(:, :)
      real(dp), allocatable :: xtx_inv(:, :)

      status = PAN_OK
      allocate(beta_mean(size(beta, 1), size(beta, 2)), xtx_inv(size(beta, 1), size(beta, 1)))
      allocate(psi_full(size(zcol) * size(beta, 2), size(zcol) * size(beta, 2)))

      do it = 1, size(result%beta, 3)
         call expand_bd_psi(psi, psi_full)
         call draw_random_effects(ywork, usable, pred, xcol, zcol, first, last, beta, sigma, psi_full, rng, b, status)
         if (status /= PAN_OK) return

         call draw_psi_bd(b, prior, rng, psi, status)
         if (status /= PAN_OK) return

         call fixed_effect_mean(ywork, usable, pred, xcol, zcol, first, last, b, beta_mean, xtx_inv, status)
         if (status /= PAN_OK) return

         call draw_sigma_beta(ywork, usable, pred, xcol, zcol, first, last, b, prior%a, prior%binv, &
            beta_mean, xtx_inv, rng, sigma, beta, status)
         if (status /= PAN_OK) return

         call impute_missing(y_original, observed, pred, xcol, zcol, first, last, beta, b, sigma, rng, ywork, status)
         if (status /= PAN_OK) return

         result%beta(:, :, it) = beta
         result%sigma(:, :, it) = sigma
         result%psi(:, :, :, it) = psi
      end do
   end subroutine run_bd_chain

   subroutine draw_random_effects(ywork, usable, pred, xcol, zcol, first, last, beta, sigma, psi, rng, b, status)
      real(dp), intent(in) :: ywork(:, :) !! Current completed response matrix.
      logical, intent(in) :: usable(:) !! True for rows used in parameter updates.
      real(dp), intent(in) :: pred(:, :) !! Complete predictor matrix.
      integer, intent(in) :: xcol(:) !! One-based fixed-effect predictor columns.
      integer, intent(in) :: zcol(:) !! One-based random-effect predictor columns.
      integer, intent(in) :: first(:) !! First row index of each subject.
      integer, intent(in) :: last(:) !! Last row index of each subject.
      real(dp), intent(in) :: beta(:, :) !! Current fixed-effect coefficient matrix.
      real(dp), intent(in) :: sigma(:, :) !! Current residual covariance matrix.
      real(dp), intent(in) :: psi(:, :) !! Current full covariance of vectorized random effects.
      type(rng_state), intent(inout) :: rng !! Mutable random-number state.
      real(dp), intent(out) :: b(:, :, :) !! Newly drawn q by r random-effect matrix for each subject.
      integer, intent(out) :: status !! PAN_OK on success or PAN_ERR_LINALG on factorization failure.

      integer :: a
      integer :: i
      integer :: ia
      integer :: j
      integer :: k
      integer :: q
      integer :: r
      integer :: s
      integer :: stat
      integer :: row
      real(dp), allocatable :: cov(:, :)
      real(dp), allocatable :: draw(:)
      real(dp), allocatable :: h(:)
      real(dp), allocatable :: mean(:)
      real(dp), allocatable :: precision(:, :)
      real(dp), allocatable :: psi_inv(:, :)
      real(dp), allocatable :: residual(:)
      real(dp), allocatable :: sigma_inv(:, :)
      real(dp), allocatable :: ztz(:, :)

      q = size(zcol)
      r = size(beta, 2)
      allocate(cov(q * r, q * r), draw(q * r), h(q * r), mean(q * r))
      allocate(precision(q * r, q * r), psi_inv(q * r, q * r), sigma_inv(r, r))
      allocate(residual(r), ztz(q, q))

      call spd_inverse(psi, psi_inv, stat)
      if (stat /= 0) then
         status = PAN_ERR_LINALG
         return
      end if
      call spd_inverse(sigma, sigma_inv, stat)
      if (stat /= 0) then
         status = PAN_ERR_LINALG
         return
      end if

      do s = 1, size(first)
         ztz = 0.0_dp
         h = 0.0_dp

         do row = first(s), last(s)
            if (.not. usable(row)) cycle

            residual = ywork(row, :) - matmul(pred(row, xcol), beta)

            do i = 1, q
               do j = 1, q
                  ztz(i, j) = ztz(i, j) + pred(row, zcol(i)) * pred(row, zcol(j))
               end do
            end do

            do a = 1, r
               do i = 1, q
                  ia = (a - 1) * q + i
                  do k = 1, r
                     h(ia) = h(ia) + pred(row, zcol(i)) * sigma_inv(a, k) * residual(k)
                  end do
               end do
            end do
         end do

         precision = psi_inv
         do a = 1, r
            do k = 1, r
               do i = 1, q
                  do j = 1, q
                     ia = (a - 1) * q + i
                     precision(ia, (k - 1) * q + j) = precision(ia, (k - 1) * q + j) + &
                        sigma_inv(a, k) * ztz(i, j)
                  end do
               end do
            end do
         end do
         call symmetrize(precision)

         call spd_inverse(precision, cov, stat)
         if (stat /= 0) then
            status = PAN_ERR_LINALG
            return
         end if
         mean = matmul(cov, h)

         call mvn_draw(rng, mean, cov, draw, stat)
         if (stat /= 0) then
            status = PAN_ERR_LINALG
            return
         end if

         do a = 1, r
            do i = 1, q
               b(i, a, s) = draw((a - 1) * q + i)
            end do
         end do
      end do

      status = PAN_OK
   end subroutine draw_random_effects

   subroutine draw_psi_full(b, prior, rng, psi, status)
      real(dp), intent(in) :: b(:, :, :) !! q by r by m subject-specific random-effect draws.
      type(pan_prior), intent(in) :: prior !! Full inverse-Wishart prior for vectorized random effects.
      type(rng_state), intent(inout) :: rng !! Mutable random-number state.
      real(dp), intent(out) :: psi(:, :) !! Newly drawn full covariance of vectorized random effects.
      integer, intent(out) :: status !! PAN_OK on success or PAN_ERR_LINALG if the draw fails.

      integer :: a
      integer :: i
      integer :: q
      integer :: r
      integer :: s
      integer :: stat
      real(dp), allocatable :: scale(:, :)
      real(dp), allocatable :: v(:)

      q = size(b, 1)
      r = size(b, 2)
      allocate(scale(q * r, q * r), v(q * r))
      scale = prior%dinv

      do s = 1, size(b, 3)
         do a = 1, r
            do i = 1, q
               v((a - 1) * q + i) = b(i, a, s)
            end do
         end do
         scale = scale + outer_product(v, v)
      end do

      call invwishart_draw(rng, scale, prior%c + real(size(b, 3), dp), psi, stat)
      if (stat /= 0) then
         status = PAN_ERR_LINALG
      else
         status = PAN_OK
      end if
   end subroutine draw_psi_full

   subroutine draw_psi_bd(b, prior, rng, psi, status)
      real(dp), intent(in) :: b(:, :, :) !! q by r by m subject-specific random-effect draws.
      type(pan_bd_prior), intent(in) :: prior !! Independent inverse-Wishart priors for each response block.
      type(rng_state), intent(inout) :: rng !! Mutable random-number state.
      real(dp), intent(out) :: psi(:, :, :) !! Newly drawn q by q covariance block for each response.
      integer, intent(out) :: status !! PAN_OK on success or PAN_ERR_LINALG if any block draw fails.

      integer :: a
      integer :: s
      integer :: stat
      real(dp), allocatable :: scale(:, :)

      allocate(scale(size(b, 1), size(b, 1)))

      do a = 1, size(b, 2)
         scale = prior%dinv(:, :, a)
         do s = 1, size(b, 3)
            scale = scale + outer_product(b(:, a, s), b(:, a, s))
         end do

         call invwishart_draw(rng, scale, prior%c(a) + real(size(b, 3), dp), psi(:, :, a), stat)
         if (stat /= 0) then
            status = PAN_ERR_LINALG
            return
         end if
      end do

      status = PAN_OK
   end subroutine draw_psi_bd

   subroutine fixed_effect_mean(ywork, usable, pred, xcol, zcol, first, last, b, beta_mean, xtx_inv, status)
      real(dp), intent(in) :: ywork(:, :) !! Current completed response matrix.
      logical, intent(in) :: usable(:) !! True for rows included in parameter updates.
      real(dp), intent(in) :: pred(:, :) !! Complete predictor matrix.
      integer, intent(in) :: xcol(:) !! One-based fixed-effect predictor columns.
      integer, intent(in) :: zcol(:) !! One-based random-effect predictor columns.
      integer, intent(in) :: first(:) !! First row index of each subject.
      integer, intent(in) :: last(:) !! Last row index of each subject.
      real(dp), intent(in) :: b(:, :, :) !! Current q by r by m random-effect draws.
      real(dp), intent(out) :: beta_mean(:, :) !! Conditional least-squares mean of the fixed effects.
      real(dp), intent(out) :: xtx_inv(:, :) !! Inverse X^T X over rows with at least one observed response.
      integer, intent(out) :: status !! PAN_OK on success or PAN_ERR_LINALG if X^T X is singular.

      integer :: row
      integer :: s
      integer :: stat
      real(dp), allocatable :: xtx(:, :)
      real(dp), allocatable :: xty(:, :)
      real(dp), allocatable :: yadj(:)

      allocate(xtx(size(xcol), size(xcol)), xty(size(xcol), size(ywork, 2)), yadj(size(ywork, 2)))
      xtx = 0.0_dp
      xty = 0.0_dp

      do s = 1, size(first)
         do row = first(s), last(s)
            if (.not. usable(row)) cycle
            yadj = ywork(row, :) - matmul(pred(row, zcol), b(:, :, s))
            xtx = xtx + outer_product(pred(row, xcol), pred(row, xcol))
            xty = xty + outer_product(pred(row, xcol), yadj)
         end do
      end do

      call spd_inverse(xtx, xtx_inv, stat)
      if (stat /= 0) then
         status = PAN_ERR_LINALG
         return
      end if

      beta_mean = matmul(xtx_inv, xty)
      status = PAN_OK
   end subroutine fixed_effect_mean

   subroutine draw_sigma_beta(ywork, usable, pred, xcol, zcol, first, last, b, a_prior, binv, &
         beta_mean, xtx_inv, rng, sigma, beta, status)
      real(dp), intent(in) :: ywork(:, :) !! Current completed response matrix.
      logical, intent(in) :: usable(:) !! True for rows included in parameter updates.
      real(dp), intent(in) :: pred(:, :) !! Complete predictor matrix.
      integer, intent(in) :: xcol(:) !! One-based fixed-effect predictor columns.
      integer, intent(in) :: zcol(:) !! One-based random-effect predictor columns.
      integer, intent(in) :: first(:) !! First row index of each subject.
      integer, intent(in) :: last(:) !! Last row index of each subject.
      real(dp), intent(in) :: b(:, :, :) !! Current q by r by m random-effect draws.
      real(dp), intent(in) :: a_prior !! Residual inverse-Wishart prior degrees of freedom.
      real(dp), intent(in) :: binv(:, :) !! Residual inverse-Wishart prior scale matrix.
      real(dp), intent(in) :: beta_mean(:, :) !! Conditional least-squares mean of fixed effects.
      real(dp), intent(in) :: xtx_inv(:, :) !! Inverse fixed-effect cross-product matrix.
      type(rng_state), intent(inout) :: rng !! Mutable random-number state.
      real(dp), intent(out) :: sigma(:, :) !! Newly drawn residual covariance matrix.
      real(dp), intent(out) :: beta(:, :) !! Newly drawn fixed-effect coefficient matrix.
      integer, intent(out) :: status !! PAN_OK on success or PAN_ERR_LINALG if a covariance draw fails.

      integer :: row
      integer :: s
      integer :: stat
      real(dp), allocatable :: e(:)
      real(dp), allocatable :: scale(:, :)

      allocate(e(size(ywork, 2)), scale(size(ywork, 2), size(ywork, 2)))
      scale = binv

      do s = 1, size(first)
         do row = first(s), last(s)
            if (.not. usable(row)) cycle
            e = ywork(row, :) - matmul(pred(row, xcol), beta_mean) - matmul(pred(row, zcol), b(:, :, s))
            scale = scale + outer_product(e, e)
         end do
      end do

      call invwishart_draw(rng, scale, a_prior + real(count(usable) - size(xcol), dp), sigma, stat)
      if (stat /= 0) then
         status = PAN_ERR_LINALG
         return
      end if

      call matrix_normal_draw(rng, beta_mean, xtx_inv, sigma, beta, stat)
      if (stat /= 0) then
         status = PAN_ERR_LINALG
         return
      end if

      status = PAN_OK
   end subroutine draw_sigma_beta

   subroutine impute_missing(y_original, observed, pred, xcol, zcol, first, last, beta, b, sigma, rng, ywork, status)
      real(dp), intent(in) :: y_original(:, :) !! Original response matrix containing observed values and NaNs.
      logical, intent(in) :: observed(:, :) !! True where y_original is observed and must remain unchanged.
      real(dp), intent(in) :: pred(:, :) !! Complete predictor matrix.
      integer, intent(in) :: xcol(:) !! One-based fixed-effect predictor columns.
      integer, intent(in) :: zcol(:) !! One-based random-effect predictor columns.
      integer, intent(in) :: first(:) !! First row index of each subject.
      integer, intent(in) :: last(:) !! Last row index of each subject.
      real(dp), intent(in) :: beta(:, :) !! Current fixed-effect coefficient matrix.
      real(dp), intent(in) :: b(:, :, :) !! Current q by r by m random-effect draws.
      real(dp), intent(in) :: sigma(:, :) !! Current residual covariance matrix.
      type(rng_state), intent(inout) :: rng !! Mutable random-number state.
      real(dp), intent(inout) :: ywork(:, :) !! Completed response matrix updated only at missing entries.
      integer, intent(out) :: status !! PAN_OK on success or PAN_ERR_LINALG if conditional covariance calculations fail.

      integer :: i
      integer :: nm
      integer :: no
      integer :: row
      integer :: s
      integer :: stat
      integer, allocatable :: miss(:)
      integer, allocatable :: obs(:)
      real(dp), allocatable :: cond_cov(:, :)
      real(dp), allocatable :: cond_mean(:)
      real(dp), allocatable :: draw(:)
      real(dp), allocatable :: mean(:)
      real(dp), allocatable :: residual_obs(:)
      real(dp), allocatable :: soo(:, :)
      real(dp), allocatable :: som(:, :)
      real(dp), allocatable :: sol(:)

      do s = 1, size(first)
         do row = first(s), last(s)
            mean = matmul(pred(row, xcol), beta) + matmul(pred(row, zcol), b(:, :, s))
            no = count(observed(row, :))
            nm = size(observed, 2) - no

            if (nm == 0) then
               ywork(row, :) = y_original(row, :)
               cycle
            end if

            allocate(obs(no), miss(nm))
            call logical_indices(observed(row, :), obs, miss)

            if (no == 0) then
               allocate(draw(size(mean)))
               call mvn_draw(rng, mean, sigma, draw, stat)
               if (stat /= 0) then
                  status = PAN_ERR_LINALG
                  return
               end if
               ywork(row, :) = draw
               deallocate(obs, miss, draw, mean)
               cycle
            end if

            allocate(soo(no, no), som(no, nm), residual_obs(no), sol(no))
            allocate(cond_mean(nm), cond_cov(nm, nm), draw(nm))

            soo = sigma(obs, obs)
            som = sigma(obs, miss)
            residual_obs = y_original(row, obs) - mean(obs)

            call spd_solve_vec(soo, residual_obs, sol, stat)
            if (stat /= 0) then
               status = PAN_ERR_LINALG
               return
            end if
            cond_mean = mean(miss) + matmul(transpose(som), sol)

            call conditional_covariance(sigma, obs, miss, cond_cov, stat)
            if (stat /= PAN_OK) then
               status = stat
               return
            end if

            call mvn_draw(rng, cond_mean, cond_cov, draw, stat)
            if (stat /= 0) then
               status = PAN_ERR_LINALG
               return
            end if

            do i = 1, no
               ywork(row, obs(i)) = y_original(row, obs(i))
            end do
            do i = 1, nm
               ywork(row, miss(i)) = draw(i)
            end do

            deallocate(obs, miss, soo, som, residual_obs, sol, cond_mean, cond_cov, draw, mean)
         end do
      end do

      status = PAN_OK
   end subroutine impute_missing

   pure subroutine conditional_covariance(sigma, obs, miss, cov, status)
      real(dp), intent(in) :: sigma(:, :) !! Full residual covariance matrix.
      integer, intent(in) :: obs(:) !! Indices of observed response components.
      integer, intent(in) :: miss(:) !! Indices of missing response components.
      real(dp), intent(out) :: cov(:, :) !! Conditional covariance of missing components given observed components.
      integer, intent(out) :: status !! PAN_OK on success or PAN_ERR_LINALG if the observed block is singular.

      integer :: stat
      real(dp), allocatable :: inv_oo_om(:, :)
      real(dp), allocatable :: soo(:, :)
      real(dp), allocatable :: som(:, :)

      allocate(soo(size(obs), size(obs)), som(size(obs), size(miss)))
      allocate(inv_oo_om(size(obs), size(miss)))
      soo = sigma(obs, obs)
      som = sigma(obs, miss)

      call solve_spd_columns(soo, som, inv_oo_om, stat)
      if (stat /= 0) then
         status = PAN_ERR_LINALG
         return
      end if

      cov = sigma(miss, miss) - matmul(transpose(som), inv_oo_om)
      call symmetrize(cov)

      if (.not. is_spd(cov)) then
         call stabilize_covariance(cov)
      end if
      if (.not. is_spd(cov)) then
         status = PAN_ERR_LINALG
      else
         status = PAN_OK
      end if
   end subroutine conditional_covariance

   pure subroutine solve_spd_columns(a, b, x, info)
      real(dp), intent(in) :: a(:, :) !! Symmetric positive-definite coefficient matrix.
      real(dp), intent(in) :: b(:, :) !! Matrix of right-hand-side columns.
      real(dp), intent(out) :: x(:, :) !! Solution matrix A^{-1} B.
      integer, intent(out) :: info !! Zero on success or the failing factorization code.

      integer :: j
      integer :: stat

      info = 0
      do j = 1, size(b, 2)
         call spd_solve_vec(a, b(:, j), x(:, j), stat)
         if (stat /= 0) then
            info = stat
            return
         end if
      end do
   end subroutine solve_spd_columns

   pure subroutine stabilize_covariance(a)
      real(dp), intent(inout) :: a(:, :) !! Nearly positive-definite covariance matrix receiving diagonal jitter.

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
   end subroutine stabilize_covariance

   pure subroutine initialize_full(y, observed, usable, pred, xcol, prior, beta, sigma, psi, ywork, status)
      real(dp), intent(in) :: y(:, :) !! Original response matrix with NaNs at missing components.
      logical, intent(in) :: observed(:, :) !! Original response-observation mask.
      logical, intent(in) :: usable(:) !! True for rows containing at least one observed response.
      real(dp), intent(in) :: pred(:, :) !! Complete predictor matrix.
      integer, intent(in) :: xcol(:) !! One-based fixed-effect predictor columns.
      type(pan_prior), intent(in) :: prior !! Full covariance prior used for initial Psi.
      real(dp), intent(out) :: beta(:, :) !! Initial fixed-effect coefficient matrix.
      real(dp), intent(out) :: sigma(:, :) !! Initial residual covariance matrix.
      real(dp), intent(out) :: psi(:, :) !! Initial full random-effect covariance matrix.
      real(dp), intent(out) :: ywork(:, :) !! Mean-imputed initial complete response matrix.
      integer, intent(out) :: status !! PAN_OK on success or an initialization failure code.

      call mean_impute(y, observed, usable, ywork, status)
      if (status /= PAN_OK) return

      call initial_ols(ywork, usable, pred, xcol, beta, sigma, status)
      if (status /= PAN_OK) return

      psi = prior%dinv / prior%c
      call symmetrize(psi)
      if (.not. is_spd(psi)) status = PAN_ERR_LINALG
   end subroutine initialize_full

   pure subroutine initialize_bd(y, observed, usable, pred, xcol, prior, beta, sigma, psi, ywork, status)
      real(dp), intent(in) :: y(:, :) !! Original response matrix with NaNs at missing components.
      logical, intent(in) :: observed(:, :) !! Original response-observation mask.
      logical, intent(in) :: usable(:) !! True for rows containing at least one observed response.
      real(dp), intent(in) :: pred(:, :) !! Complete predictor matrix.
      integer, intent(in) :: xcol(:) !! One-based fixed-effect predictor columns.
      type(pan_bd_prior), intent(in) :: prior !! Block covariance prior used for initial Psi blocks.
      real(dp), intent(out) :: beta(:, :) !! Initial fixed-effect coefficient matrix.
      real(dp), intent(out) :: sigma(:, :) !! Initial residual covariance matrix.
      real(dp), intent(out) :: psi(:, :, :) !! Initial q by q random-effect covariance blocks.
      real(dp), intent(out) :: ywork(:, :) !! Mean-imputed initial complete response matrix.
      integer, intent(out) :: status !! PAN_OK on success or an initialization failure code.

      integer :: j

      call mean_impute(y, observed, usable, ywork, status)
      if (status /= PAN_OK) return

      call initial_ols(ywork, usable, pred, xcol, beta, sigma, status)
      if (status /= PAN_OK) return

      do j = 1, size(psi, 3)
         psi(:, :, j) = prior%dinv(:, :, j) / prior%c(j)
         call symmetrize(psi(:, :, j))
         if (.not. is_spd(psi(:, :, j))) then
            status = PAN_ERR_LINALG
            return
         end if
      end do
   end subroutine initialize_bd

   pure subroutine mean_impute(y, observed, usable, ywork, status)
      real(dp), intent(in) :: y(:, :) !! Original response matrix with NaNs at missing components.
      logical, intent(in) :: observed(:, :) !! Original response-observation mask.
      logical, intent(in) :: usable(:) !! True for rows containing at least one observed response.
      real(dp), intent(out) :: ywork(:, :) !! Response matrix with column means substituted for missing components.
      integer, intent(out) :: status !! PAN_OK on success or PAN_ERR_ARGUMENT if a response column is entirely missing.

      integer :: i
      integer :: j
      integer :: nobs
      real(dp) :: mean_value

      ywork = 0.0_dp
      do j = 1, size(y, 2)
         nobs = 0
         mean_value = 0.0_dp
         do i = 1, size(y, 1)
            if (usable(i) .and. observed(i, j)) then
               mean_value = mean_value + y(i, j)
               nobs = nobs + 1
            end if
         end do
         if (nobs == 0) then
            status = PAN_ERR_ARGUMENT
            return
         end if
         mean_value = mean_value / real(nobs, dp)

         do i = 1, size(y, 1)
            if (observed(i, j)) then
               ywork(i, j) = y(i, j)
            else
               ywork(i, j) = mean_value
            end if
         end do
      end do

      status = PAN_OK
   end subroutine mean_impute

   pure subroutine initial_ols(ywork, usable, pred, xcol, beta, sigma, status)
      real(dp), intent(in) :: ywork(:, :) !! Mean-imputed complete response matrix.
      logical, intent(in) :: usable(:) !! True for rows containing at least one originally observed response.
      real(dp), intent(in) :: pred(:, :) !! Complete predictor matrix.
      integer, intent(in) :: xcol(:) !! One-based fixed-effect predictor columns.
      real(dp), intent(out) :: beta(:, :) !! Ordinary least-squares starting coefficients.
      real(dp), intent(out) :: sigma(:, :) !! Residual cross-product starting covariance.
      integer, intent(out) :: status !! PAN_OK on success or PAN_ERR_LINALG if X^T X is singular.

      integer :: i
      integer :: stat
      real(dp), allocatable :: e(:)
      real(dp), allocatable :: inv_xtx(:, :)
      real(dp), allocatable :: xtx(:, :)
      real(dp), allocatable :: xty(:, :)

      allocate(xtx(size(xcol), size(xcol)), xty(size(xcol), size(ywork, 2)))
      allocate(inv_xtx(size(xcol), size(xcol)), e(size(ywork, 2)))
      xtx = 0.0_dp
      xty = 0.0_dp

      do i = 1, size(ywork, 1)
         if (.not. usable(i)) cycle
         xtx = xtx + outer_product(pred(i, xcol), pred(i, xcol))
         xty = xty + outer_product(pred(i, xcol), ywork(i, :))
      end do

      call spd_inverse(xtx, inv_xtx, stat)
      if (stat /= 0) then
         status = PAN_ERR_LINALG
         return
      end if
      beta = matmul(inv_xtx, xty)

      sigma = 0.0_dp
      do i = 1, size(ywork, 1)
         if (.not. usable(i)) cycle
         e = ywork(i, :) - matmul(pred(i, xcol), beta)
         sigma = sigma + outer_product(e, e)
      end do
      sigma = sigma / real(count(usable), dp)
      call symmetrize(sigma)
      if (.not. is_spd(sigma)) call stabilize_covariance(sigma)

      if (.not. is_spd(sigma)) then
         status = PAN_ERR_LINALG
      else
         status = PAN_OK
      end if
   end subroutine initial_ols

   pure subroutine load_start_full(start, y, observed, p, q, r, beta, sigma, psi, ywork, status)
      type(pan_state), intent(in) :: start !! User-supplied restart state.
      real(dp), intent(in) :: y(:, :) !! Original response matrix used to verify observed restart entries.
      logical, intent(in) :: observed(:, :) !! Original response-observation mask.
      integer, intent(in) :: p !! Number of fixed-effect predictors.
      integer, intent(in) :: q !! Number of random-effect predictors.
      integer, intent(in) :: r !! Number of response variables.
      real(dp), intent(out) :: beta(:, :) !! Loaded fixed-effect coefficient matrix.
      real(dp), intent(out) :: sigma(:, :) !! Loaded residual covariance matrix.
      real(dp), intent(out) :: psi(:, :) !! Loaded full random-effect covariance matrix.
      real(dp), intent(out) :: ywork(:, :) !! Loaded completed response matrix with observed entries restored.
      integer, intent(out) :: status !! PAN_OK on success or an argument/dimension error.

      integer :: i
      integer :: j

      if (.not. allocated(start%beta) .or. .not. allocated(start%sigma) .or. &
          .not. allocated(start%psi) .or. .not. allocated(start%y)) then
         status = PAN_ERR_ARGUMENT
         return
      end if
      if (any(shape(start%beta) /= [p, r]) .or. any(shape(start%sigma) /= [r, r]) .or. &
          any(shape(start%psi) /= [q * r, q * r]) .or. any(shape(start%y) /= shape(y))) then
         status = PAN_ERR_DIMENSION
         return
      end if
      if (.not. is_spd(start%sigma) .or. .not. is_spd(start%psi)) then
         status = PAN_ERR_ARGUMENT
         return
      end if

      beta = start%beta
      sigma = start%sigma
      psi = start%psi
      ywork = start%y
      do j = 1, r
         do i = 1, size(y, 1)
            if (observed(i, j)) ywork(i, j) = y(i, j)
         end do
      end do
      status = PAN_OK
   end subroutine load_start_full

   pure subroutine load_start_bd(start, y, observed, p, q, r, beta, sigma, psi, ywork, status)
      type(pan_bd_state), intent(in) :: start !! User-supplied block-diagonal restart state.
      real(dp), intent(in) :: y(:, :) !! Original response matrix used to verify observed restart entries.
      logical, intent(in) :: observed(:, :) !! Original response-observation mask.
      integer, intent(in) :: p !! Number of fixed-effect predictors.
      integer, intent(in) :: q !! Number of random-effect predictors.
      integer, intent(in) :: r !! Number of response variables.
      real(dp), intent(out) :: beta(:, :) !! Loaded fixed-effect coefficient matrix.
      real(dp), intent(out) :: sigma(:, :) !! Loaded residual covariance matrix.
      real(dp), intent(out) :: psi(:, :, :) !! Loaded q by q covariance block for each response.
      real(dp), intent(out) :: ywork(:, :) !! Loaded completed response matrix with observed entries restored.
      integer, intent(out) :: status !! PAN_OK on success or an argument/dimension error.

      integer :: a
      integer :: i
      integer :: j

      if (.not. allocated(start%beta) .or. .not. allocated(start%sigma) .or. &
          .not. allocated(start%psi) .or. .not. allocated(start%y)) then
         status = PAN_ERR_ARGUMENT
         return
      end if
      if (any(shape(start%beta) /= [p, r]) .or. any(shape(start%sigma) /= [r, r]) .or. &
          any(shape(start%psi) /= [q, q, r]) .or. any(shape(start%y) /= shape(y))) then
         status = PAN_ERR_DIMENSION
         return
      end if
      if (.not. is_spd(start%sigma)) then
         status = PAN_ERR_ARGUMENT
         return
      end if
      do a = 1, r
         if (.not. is_spd(start%psi(:, :, a))) then
            status = PAN_ERR_ARGUMENT
            return
         end if
      end do

      beta = start%beta
      sigma = start%sigma
      psi = start%psi
      ywork = start%y
      do j = 1, r
         do i = 1, size(y, 1)
            if (observed(i, j)) ywork(i, j) = y(i, j)
         end do
      end do
      status = PAN_OK
   end subroutine load_start_bd

   pure subroutine expand_bd_psi(psi, full)
      real(dp), intent(in) :: psi(:, :, :) !! q by q covariance block for each response.
      real(dp), intent(out) :: full(:, :) !! Block-diagonal full covariance in response-major vectorization order.

      integer :: a
      integer :: q
      integer :: lo
      integer :: hi

      q = size(psi, 1)
      full = 0.0_dp
      do a = 1, size(psi, 3)
         lo = (a - 1) * q + 1
         hi = a * q
         full(lo:hi, lo:hi) = psi(:, :, a)
      end do
   end subroutine expand_bd_psi

   pure function outer_product(x, y) result(a)
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
   end function outer_product

   pure subroutine logical_indices(mask, true_idx, false_idx)
      logical, intent(in) :: mask(:) !! Logical mask to split into true and false index arrays.
      integer, intent(out) :: true_idx(:) !! One-based positions where mask is true.
      integer, intent(out) :: false_idx(:) !! One-based positions where mask is false.

      integer :: i
      integer :: nf
      integer :: nt

      nt = 0
      nf = 0
      do i = 1, size(mask)
         if (mask(i)) then
            nt = nt + 1
            true_idx(nt) = i
         else
            nf = nf + 1
            false_idx(nf) = i
         end if
      end do
   end subroutine logical_indices

   pure subroutine validate_common(y, subj, pred, xcol, zcol, iter, status)
      real(dp), intent(in) :: y(:, :) !! Response matrix whose dimensions are validated.
      integer, intent(in) :: subj(:) !! Subject labels whose ordering and length are validated.
      real(dp), intent(in) :: pred(:, :) !! Predictor matrix whose dimensions and finiteness are validated.
      integer, intent(in) :: xcol(:) !! One-based fixed-effect predictor columns.
      integer, intent(in) :: zcol(:) !! One-based random-effect predictor columns.
      integer, intent(in) :: iter !! Requested number of Gibbs iterations.
      integer, intent(out) :: status !! PAN_OK if all common structural checks pass.

      integer :: i

      status = PAN_OK
      if (size(y, 1) <= 0 .or. size(y, 2) <= 0 .or. size(xcol) <= 0 .or. size(zcol) <= 0) then
         status = PAN_ERR_DIMENSION
         return
      end if
      if (size(subj) /= size(y, 1) .or. size(pred, 1) /= size(y, 1)) then
         status = PAN_ERR_DIMENSION
         return
      end if
      if (iter <= 0) then
         status = PAN_ERR_ARGUMENT
         return
      end if
      if (any(xcol < 1) .or. any(xcol > size(pred, 2)) .or. any(zcol < 1) .or. any(zcol > size(pred, 2))) then
         status = PAN_ERR_ARGUMENT
         return
      end if
      if (any(ieee_is_nan(pred))) then
         status = PAN_ERR_ARGUMENT
         return
      end if
      do i = 2, size(subj)
         if (subj(i) < subj(i - 1)) then
            status = PAN_ERR_ARGUMENT
            return
         end if
      end do
   end subroutine validate_common

   pure subroutine cluster_bounds(subj, first, last, m, status)
      integer, intent(in) :: subj(:) !! Sorted subject labels grouped contiguously.
      integer, allocatable, intent(out) :: first(:) !! First row of each subject block.
      integer, allocatable, intent(out) :: last(:) !! Last row of each subject block.
      integer, intent(out) :: m !! Number of distinct contiguous subject blocks.
      integer, intent(out) :: status !! PAN_OK on success or PAN_ERR_ARGUMENT for decreasing labels.

      integer :: i
      integer :: s

      status = PAN_OK
      if (size(subj) == 0) then
         allocate(first(0), last(0))
         m = 0
         return
      end if

      do i = 2, size(subj)
         if (subj(i) < subj(i - 1)) then
            status = PAN_ERR_ARGUMENT
            allocate(first(0), last(0))
            m = 0
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
   end subroutine cluster_bounds

   subroutine copy_full_state(beta, sigma, psi, y, state)
      real(dp), intent(in) :: beta(:, :) !! Final fixed-effect coefficient matrix.
      real(dp), intent(in) :: sigma(:, :) !! Final residual covariance matrix.
      real(dp), intent(in) :: psi(:, :) !! Final full random-effect covariance matrix.
      real(dp), intent(in) :: y(:, :) !! Final completed response matrix.
      type(pan_state), intent(out) :: state !! Allocatable restart state receiving independent copies.

      state%beta = beta
      state%sigma = sigma
      state%psi = psi
      state%y = y
   end subroutine copy_full_state

   subroutine copy_bd_state(beta, sigma, psi, y, state)
      real(dp), intent(in) :: beta(:, :) !! Final fixed-effect coefficient matrix.
      real(dp), intent(in) :: sigma(:, :) !! Final residual covariance matrix.
      real(dp), intent(in) :: psi(:, :, :) !! Final block-diagonal random-effect covariance blocks.
      real(dp), intent(in) :: y(:, :) !! Final completed response matrix.
      type(pan_bd_state), intent(out) :: state !! Allocatable restart state receiving independent copies.

      state%beta = beta
      state%sigma = sigma
      state%psi = psi
      state%y = y
   end subroutine copy_bd_state

   subroutine clear_pan_result(result)
      type(pan_result), intent(out) :: result !! Result object reset to a predictable empty successful state.

      result%status = PAN_OK
      result%message = ""
   end subroutine clear_pan_result

   subroutine clear_pan_bd_result(result)
      type(pan_bd_result), intent(out) :: result !! Block-diagonal result object reset to a predictable empty state.

      result%status = PAN_OK
      result%message = ""
   end subroutine clear_pan_bd_result

   subroutine set_pan_error(result, status, message)
      type(pan_result), intent(inout) :: result !! Result object receiving the error status and message.
      integer, intent(in) :: status !! Nonzero package status code.
      character(len=*), intent(in) :: message !! Human-readable diagnostic explaining the failure.

      result%status = status
      result%message = message
   end subroutine set_pan_error

   subroutine set_pan_bd_error(result, status, message)
      type(pan_bd_result), intent(inout) :: result !! Block result receiving the error status and message.
      integer, intent(in) :: status !! Nonzero package status code.
      character(len=*), intent(in) :: message !! Human-readable diagnostic explaining the failure.

      result%status = status
      result%message = message
   end subroutine set_pan_bd_error

end module pan_sampler
