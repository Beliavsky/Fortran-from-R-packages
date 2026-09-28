module psych_factor_extended
   use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
   use r_kinds, only : dp
   use r_linalg, only : symmetric_eigen
   use psych_correlations, only : pseudoinverse, squared_multiple_correlations
   use psych_distributions, only : chi_square_cdf
   use psych_factor, only : factor_residuals, procrustes_rotate
   implicit none
   private
   public :: minimum_residual_factor, weighted_least_squares_factor
   public :: factor_fit_statistics, factor_regression_partition, factor_extension_loadings
   public :: equamax_rotate, varimin_rotate, varimin_objective_gradient
   public :: bimin_objective_gradient, target_objective_gradient
   public :: bifactor_rotate, biquartimin_rotate, target_orthogonal_rotate
   public :: sort_factor_loadings, factor_to_irt_parameters, eigenvalue_factor_count

contains

   subroutine minimum_residual_factor(r, n_factors, loadings, uniqueness, objective, info, max_iter, tolerance)
      real(dp), intent(in) :: r(:, :) !! Correlation matrix to approximate by a common-factor model.
      integer, intent(in) :: n_factors !! Number of common factors.
      real(dp), allocatable, intent(out) :: loadings(:, :) !! MINRES item-by-factor loading matrix.
      real(dp), allocatable, intent(out) :: uniqueness(:) !! Optimized uniqueness estimates constrained to `[0.005,1]`.
      real(dp), intent(out) :: objective !! Sum of squared off-diagonal residual correlations.
      integer, intent(out) :: info !! Zero when all eigensolutions succeed.
      integer, intent(in), optional :: max_iter !! Maximum coordinate-search sweeps; default 80.
      real(dp), intent(in), optional :: tolerance !! Minimum uniqueness-search step; default `1e-5`.
      real(dp), allocatable :: smc(:)
      integer :: p, smc_info

      p = size(r, 1)
      allocate(uniqueness(p), smc(p))
      call squared_multiple_correlations(r, smc, smc_info)
      if (smc_info == 0) then
         uniqueness = min(1.0_dp, max(0.005_dp, 1.0_dp - smc))
      else
         uniqueness = 0.5_dp
      end if
      call optimize_uniqueness(r, n_factors, uniqueness, .false., loadings, objective, info, max_iter, tolerance)
   end subroutine minimum_residual_factor

   subroutine weighted_least_squares_factor(r, n_factors, loadings, uniqueness, objective, info, max_iter, tolerance)
      real(dp), intent(in) :: r(:, :) !! Correlation matrix to approximate by weighted least squares.
      integer, intent(in) :: n_factors !! Number of common factors.
      real(dp), allocatable, intent(out) :: loadings(:, :) !! WLS item-by-factor loading matrix.
      real(dp), allocatable, intent(out) :: uniqueness(:) !! Optimized uniqueness estimates constrained to `[0.005,1]`.
      real(dp), intent(out) :: objective !! Weighted sum of squared off-diagonal residual correlations.
      integer, intent(out) :: info !! Zero when inversion and eigensolutions succeed.
      integer, intent(in), optional :: max_iter !! Maximum coordinate-search sweeps; default 80.
      real(dp), intent(in), optional :: tolerance !! Minimum uniqueness-search step; default `1e-5`.
      real(dp), allocatable :: smc(:)
      integer :: p, smc_info

      p = size(r, 1)
      allocate(uniqueness(p), smc(p))
      call squared_multiple_correlations(r, smc, smc_info)
      if (smc_info == 0) then
         uniqueness = min(1.0_dp, max(0.005_dp, 1.0_dp - smc))
      else
         uniqueness = 0.5_dp
      end if
      call optimize_uniqueness(r, n_factors, uniqueness, .true., loadings, objective, info, max_iter, tolerance)
   end subroutine weighted_least_squares_factor

   subroutine factor_fit_statistics(r, loadings, n_obs, residual_rms, offdiag_fit, chisq, df, p_value)
      real(dp), intent(in) :: r(:, :) !! Observed correlation matrix.
      real(dp), intent(in) :: loadings(:, :) !! Item-by-factor loading matrix.
      integer, intent(in) :: n_obs !! Sample size for the approximate chi-square statistic.
      real(dp), intent(out) :: residual_rms !! Root mean squared off-diagonal residual correlation.
      real(dp), intent(out) :: offdiag_fit !! One minus residual off-diagonal energy divided by observed off-diagonal energy.
      real(dp), intent(out) :: chisq !! Approximate residual chi-square statistic.
      real(dp), intent(out) :: df !! Factor-model residual degrees of freedom.
      real(dp), intent(out) :: p_value !! Upper-tail chi-square probability.
      real(dp), allocatable :: residual(:, :)
      real(dp) :: observed_ss, residual_ss
      integer :: i, j, p, k, n_off

      p = size(r, 1)
      k = size(loadings, 2)
      allocate(residual(p, p))
      call factor_residuals(r, loadings, residual)
      observed_ss = 0.0_dp
      residual_ss = 0.0_dp
      n_off = 0
      do j = 1, p
         do i = j + 1, p
            observed_ss = observed_ss + r(i, j) * r(i, j)
            residual_ss = residual_ss + residual(i, j) * residual(i, j)
            n_off = n_off + 1
         end do
      end do
      if (n_off > 0) then
         residual_rms = sqrt(residual_ss / real(n_off, dp))
      else
         residual_rms = 0.0_dp
      end if
      if (observed_ss > tiny(1.0_dp)) then
         offdiag_fit = 1.0_dp - residual_ss / observed_ss
      else
         offdiag_fit = 1.0_dp
      end if
      df = real(max(0, ((p - k) * (p - k) - (p + k)) / 2), dp)
      chisq = real(max(0, n_obs - 1), dp) * residual_ss
      if (df > 0.0_dp) then
         p_value = 1.0_dp - chi_square_cdf(max(0.0_dp, chisq), df)
      else
         p_value = 1.0_dp
      end if
   end subroutine factor_fit_statistics

   subroutine factor_regression_partition(r, predictors, outcomes, beta, residual_covariance, multiple_r2, info)
      real(dp), intent(in) :: r(:, :) !! Full correlation or covariance matrix.
      integer, intent(in) :: predictors(:) !! One-based predictor-variable indices.
      integer, intent(in) :: outcomes(:) !! One-based outcome-variable indices.
      real(dp), allocatable, intent(out) :: beta(:, :) !! Predictor-by-outcome standardized regression coefficients.
      real(dp), allocatable, intent(out) :: residual_covariance(:, :) !! Outcome residual covariance matrix.
      real(dp), allocatable, intent(out) :: multiple_r2(:) !! Squared multiple correlation for each outcome.
      integer, intent(out) :: info !! Zero when predictor block inversion succeeds.
      real(dp), allocatable :: rxx(:, :), rxy(:, :), ryy(:, :), inv(:, :)
      integer :: i, j, nx, ny

      nx = size(predictors)
      ny = size(outcomes)
      allocate(rxx(nx, nx), rxy(nx, ny), ryy(ny, ny))
      do j = 1, nx
         do i = 1, nx
            rxx(i, j) = r(predictors(i), predictors(j))
         end do
      end do
      do j = 1, ny
         do i = 1, nx
            rxy(i, j) = r(predictors(i), outcomes(j))
         end do
         do i = 1, ny
            ryy(i, j) = r(outcomes(i), outcomes(j))
         end do
      end do
      call pseudoinverse(rxx, inv, info=info)
      allocate(beta(nx, ny), residual_covariance(ny, ny), multiple_r2(ny))
      if (info /= 0) then
         beta = 0.0_dp
         residual_covariance = ryy
         multiple_r2 = 0.0_dp
         return
      end if
      beta = matmul(inv, rxy)
      residual_covariance = ryy - matmul(transpose(rxy), beta)
      do j = 1, ny
         if (abs(ryy(j, j)) > tiny(1.0_dp)) then
            multiple_r2(j) = 1.0_dp - residual_covariance(j, j) / ryy(j, j)
         else
            multiple_r2(j) = 0.0_dp
         end if
      end do
   end subroutine factor_regression_partition

   subroutine factor_extension_loadings(r_extension, r_factor, loadings, extended, info)
      real(dp), intent(in) :: r_extension(:, :) !! Correlations of extension variables with factor-analysis variables.
      real(dp), intent(in) :: r_factor(:, :) !! Correlation matrix used to fit the original factors.
      real(dp), intent(in) :: loadings(:, :) !! Original item-by-factor loading matrix.
      real(dp), allocatable, intent(out) :: extended(:, :) !! Extension-variable factor loadings.
      integer, intent(out) :: info !! Zero when pseudoinversion succeeds.
      real(dp), allocatable :: rinv(:, :), weights(:, :)

      call pseudoinverse(r_factor, rinv, info=info)
      allocate(extended(size(r_extension, 1), size(loadings, 2)))
      if (info /= 0) then
         extended = 0.0_dp
         return
      end if
      weights = matmul(rinv, loadings)
      extended = matmul(r_extension, weights)
   end subroutine factor_extension_loadings

   subroutine equamax_rotate(loadings, rotated, rotation, info, max_iter, tolerance)
      real(dp), intent(in) :: loadings(:, :) !! Loading matrix to rotate orthogonally.
      real(dp), allocatable, intent(out) :: rotated(:, :) !! Equamax-rotated loadings.
      real(dp), allocatable, intent(out) :: rotation(:, :) !! Orthogonal rotation matrix.
      integer, intent(out) :: info !! Zero after the deterministic pair-search optimizer completes.
      integer, intent(in), optional :: max_iter !! Maximum pair-search sweeps; default 60.
      real(dp), intent(in), optional :: tolerance !! Minimum accepted angle magnitude; default `1e-7`.
      real(dp) :: kappa

      kappa = real(size(loadings, 2), dp) / (2.0_dp * real(max(1, size(loadings, 1)), dp))
      call orthogonal_cf_rotate(loadings, kappa, .false., rotated, rotation, info, max_iter, tolerance)
   end subroutine equamax_rotate

   subroutine varimin_rotate(loadings, rotated, rotation, info, max_iter, tolerance)
      real(dp), intent(in) :: loadings(:, :) !! Loading matrix to rotate orthogonally.
      real(dp), allocatable, intent(out) :: rotated(:, :) !! Varimin-rotated loadings that minimize variance of squared loadings.
      real(dp), allocatable, intent(out) :: rotation(:, :) !! Orthogonal rotation matrix.
      integer, intent(out) :: info !! Zero after the deterministic pair-search optimizer completes.
      integer, intent(in), optional :: max_iter !! Maximum pair-search sweeps; default 60.
      real(dp), intent(in), optional :: tolerance !! Minimum accepted angle magnitude; default `1e-7`.

      call orthogonal_cf_rotate(loadings, 1.0_dp / real(max(1, size(loadings, 1)), dp), .true., &
         rotated, rotation, info, max_iter, tolerance)
   end subroutine varimin_rotate

   pure subroutine varimin_objective_gradient(loadings, objective, gradient)
      real(dp), intent(in) :: loadings(:, :) !! Loading matrix at which to evaluate the varimin criterion.
      real(dp), intent(out) :: objective !! Varimin objective from centered squared loadings.
      real(dp), intent(out) :: gradient(:, :) !! Gradient with respect to the loading matrix.
      real(dp), allocatable :: q(:, :), means(:)
      integer :: j

      allocate(q(size(loadings, 1), size(loadings, 2)), means(size(loadings, 2)))
      do j = 1, size(loadings, 2)
         means(j) = sum(loadings(:, j) * loadings(:, j)) / real(max(1, size(loadings, 1)), dp)
         q(:, j) = loadings(:, j) * loadings(:, j) - means(j)
      end do
      objective = 0.25_dp * sum(q * q)
      gradient = loadings * q
   end subroutine varimin_objective_gradient

   pure subroutine bimin_objective_gradient(loadings, objective, gradient)
      real(dp), intent(in) :: loadings(:, :) !! Loading matrix whose first column is the general factor.
      real(dp), intent(out) :: objective !! Jennrich-Bentler bi-factor quartimin objective.
      real(dp), intent(out) :: gradient(:, :) !! Gradient with zero first-factor column.
      real(dp), allocatable :: l2(:, :), cross(:, :)
      integer :: i, j, k

      gradient = 0.0_dp
      if (size(loadings, 2) <= 1) then
         objective = 0.0_dp
         return
      end if
      allocate(l2(size(loadings, 1), size(loadings, 2) - 1))
      allocate(cross(size(loadings, 1), size(loadings, 2) - 1))
      l2 = loadings(:, 2:) * loadings(:, 2:)
      cross = 0.0_dp
      do j = 1, size(l2, 2)
         do k = 1, size(l2, 2)
            if (k == j) cycle
            cross(:, j) = cross(:, j) + l2(:, k)
         end do
      end do
      objective = sum(l2 * cross)
      do j = 1, size(l2, 2)
         i = j + 1
         gradient(:, i) = 4.0_dp * loadings(:, i) * cross(:, j)
      end do
   end subroutine bimin_objective_gradient

   pure subroutine target_objective_gradient(loadings, target, objective, gradient)
      real(dp), intent(in) :: loadings(:, :) !! Current loading matrix.
      real(dp), intent(in) :: target(:, :) !! Target loadings; IEEE NaNs mark unconstrained cells.
      real(dp), intent(out) :: objective !! Sum of squared constrained-cell deviations.
      real(dp), intent(out) :: gradient(:, :) !! Gradient, with zero entries for unconstrained target cells.
      integer :: i, j

      objective = 0.0_dp
      gradient = 0.0_dp
      do j = 1, min(size(loadings, 2), size(target, 2))
         do i = 1, min(size(loadings, 1), size(target, 1))
            if (ieee_is_nan(target(i, j))) cycle
            objective = objective + (loadings(i, j) - target(i, j)) ** 2
            gradient(i, j) = 2.0_dp * (loadings(i, j) - target(i, j))
         end do
      end do
   end subroutine target_objective_gradient

   subroutine bifactor_rotate(loadings, rotated, rotation, info, max_iter, tolerance)
      real(dp), intent(in) :: loadings(:, :) !! Loading matrix whose first factor is retained as the general factor.
      real(dp), allocatable, intent(out) :: rotated(:, :) !! Orthogonally rotated bifactor loadings.
      real(dp), allocatable, intent(out) :: rotation(:, :) !! Orthogonal transformation matrix.
      integer, intent(out) :: info !! Zero after deterministic optimization.
      integer, intent(in), optional :: max_iter !! Maximum pair-search sweeps; default 60.
      real(dp), intent(in), optional :: tolerance !! Minimum accepted angle; default `1e-7`.

      call orthogonal_bimin_rotate(loadings, rotated, rotation, info, max_iter, tolerance)
   end subroutine bifactor_rotate

   subroutine biquartimin_rotate(loadings, rotated, transformation, phi, info, max_iter, tolerance)
      real(dp), intent(in) :: loadings(:, :) !! Loading matrix to rotate obliquely toward bifactor simple structure.
      real(dp), allocatable, intent(out) :: rotated(:, :) !! Obliquely rotated biquartimin loadings.
      real(dp), allocatable, intent(out) :: transformation(:, :) !! Nonsingular loading transformation matrix.
      real(dp), allocatable, intent(out) :: phi(:, :) !! Resulting factor-correlation matrix.
      integer, intent(out) :: info !! Zero when the final transformation can be pseudoinverted.
      integer, intent(in), optional :: max_iter !! Maximum coordinate-search sweeps; default 80.
      real(dp), intent(in), optional :: tolerance !! Minimum coordinate step; default `1e-5`.

      call oblique_bimin_coordinate_search(loadings, rotated, transformation, phi, info, max_iter, tolerance)
   end subroutine biquartimin_rotate

   subroutine target_orthogonal_rotate(loadings, target, rotated, rotation, info)
      real(dp), intent(in) :: loadings(:, :) !! Loading matrix to rotate.
      real(dp), intent(in) :: target(:, :) !! Complete orthogonal target matrix.
      real(dp), allocatable, intent(out) :: rotated(:, :) !! Orthogonally target-rotated loadings.
      real(dp), allocatable, intent(out) :: rotation(:, :) !! Orthogonal target rotation.
      integer, intent(out) :: info !! SVD status from the Procrustes solution.

      allocate(rotated(size(loadings, 1), size(loadings, 2)))
      allocate(rotation(size(loadings, 2), size(loadings, 2)))
      call procrustes_rotate(loadings, target, rotated, rotation, info)
   end subroutine target_orthogonal_rotate

   pure subroutine sort_factor_loadings(loadings, sorted, order, flip_signs)
      real(dp), intent(in) :: loadings(:, :) !! Item-by-factor loading matrix.
      real(dp), allocatable, intent(out) :: sorted(:, :) !! Factors ordered by descending sum of squared loadings.
      integer, allocatable, intent(out) :: order(:) !! Original one-based factor indices in sorted order.
      logical, intent(in), optional :: flip_signs !! Flip columns to make their sums nonnegative; default true.
      real(dp), allocatable :: ss(:)
      real(dp) :: ssv
      integer :: i, j, idx, k
      logical :: flip

      k = size(loadings, 2)
      allocate(sorted(size(loadings, 1), k), order(k), ss(k))
      do j = 1, k
         order(j) = j
         ss(j) = sum(loadings(:, j) * loadings(:, j))
      end do
      do i = 1, k - 1
         idx = i
         do j = i + 1, k
            if (ss(j) > ss(idx)) idx = j
         end do
         if (idx /= i) then
            ssv = ss(i)
            ss(i) = ss(idx)
            ss(idx) = ssv
            j = order(i)
            order(i) = order(idx)
            order(idx) = j
         end if
      end do
      do j = 1, k
         sorted(:, j) = loadings(:, order(j))
      end do
      flip = .true.
      if (present(flip_signs)) flip = flip_signs
      if (flip) then
         do j = 1, k
            if (sum(sorted(:, j)) < 0.0_dp) sorted(:, j) = -sorted(:, j)
         end do
      end if
   end subroutine sort_factor_loadings

   pure subroutine factor_to_irt_parameters(loadings, thresholds, discrimination, difficulty)
      real(dp), intent(in) :: loadings(:, :) !! Item-by-factor loadings on a latent-normal response model.
      real(dp), intent(in) :: thresholds(:, :) !! Item-by-threshold latent-normal cut points.
      real(dp), allocatable, intent(out) :: discrimination(:, :) !! Item-by-factor IRT discriminations `l/sqrt(1-l^2)`.
      real(dp), allocatable, intent(out) :: difficulty(:, :, :) !! Item-by-threshold-by-factor transformed difficulties.
      real(dp) :: denom
      integer :: i, j, f

      allocate(discrimination(size(loadings, 1), size(loadings, 2)))
      allocate(difficulty(size(thresholds, 1), size(thresholds, 2), size(loadings, 2)))
      do f = 1, size(loadings, 2)
         do i = 1, size(loadings, 1)
            denom = sqrt(max(tiny(1.0_dp), 1.0_dp - loadings(i, f) * loadings(i, f)))
            discrimination(i, f) = loadings(i, f) / denom
            do j = 1, size(thresholds, 2)
               difficulty(i, j, f) = thresholds(i, j) / denom
            end do
         end do
      end do
   end subroutine factor_to_irt_parameters

   pure integer function eigenvalue_factor_count(values, reference) result(n_factors)
      real(dp), intent(in) :: values(:) !! Observed descending eigenvalues.
      real(dp), intent(in), optional :: reference(:) !! Optional parallel-analysis reference eigenvalues.
      integer :: i, n

      n = size(values)
      if (present(reference)) n = min(n, size(reference))
      n_factors = 0
      do i = 1, n
         if (present(reference)) then
            if (values(i) > reference(i)) n_factors = n_factors + 1
         else
            if (values(i) > 1.0_dp) n_factors = n_factors + 1
         end if
      end do
   end function eigenvalue_factor_count

   subroutine optimize_uniqueness(r, n_factors, uniqueness, weighted, loadings, objective, info, max_iter, tolerance)
      real(dp), intent(in) :: r(:, :) !! Correlation matrix.
      integer, intent(in) :: n_factors !! Number of factors.
      real(dp), intent(inout) :: uniqueness(:) !! Starting and optimized uniqueness vector.
      logical, intent(in) :: weighted !! Use WLS residual weights when true.
      real(dp), allocatable, intent(out) :: loadings(:, :) !! Loading matrix at the optimized uniqueness vector.
      real(dp), intent(out) :: objective !! Final residual objective value.
      integer, intent(out) :: info !! Zero when all eigensolutions and weight inversions succeed.
      integer, intent(in), optional :: max_iter !! Maximum coordinate sweeps.
      real(dp), intent(in), optional :: tolerance !! Minimum coordinate step.
      real(dp), allocatable :: trial_loadings(:, :), best_loadings(:, :), weights(:)
      real(dp) :: best, candidate, old_value, step, tol
      integer :: i, iter, limit, eval_info
      logical :: improved

      limit = 80
      if (present(max_iter)) limit = max(1, max_iter)
      tol = 1.0e-5_dp
      if (present(tolerance)) tol = max(1.0e-10_dp, tolerance)
      allocate(weights(size(r, 1)))
      call residual_weights(r, weighted, weights, info)
      if (info /= 0) then
         allocate(loadings(size(r, 1), min(max(1, n_factors), size(r, 1))))
         loadings = 0.0_dp
         objective = huge(1.0_dp)
         return
      end if
      call factor_from_uniqueness(r, uniqueness, n_factors, best_loadings, best, weights, eval_info)
      info = eval_info
      if (info /= 0) then
         loadings = best_loadings
         objective = best
         return
      end if
      step = 0.10_dp
      do iter = 1, limit
         improved = .false.
         do i = 1, size(uniqueness)
            old_value = uniqueness(i)
            uniqueness(i) = min(1.0_dp, old_value + step)
            call factor_from_uniqueness(r, uniqueness, n_factors, trial_loadings, candidate, weights, eval_info)
            if (eval_info == 0 .and. candidate < best) then
               best = candidate
               best_loadings = trial_loadings
               improved = .true.
               cycle
            end if
            uniqueness(i) = max(0.005_dp, old_value - step)
            call factor_from_uniqueness(r, uniqueness, n_factors, trial_loadings, candidate, weights, eval_info)
            if (eval_info == 0 .and. candidate < best) then
               best = candidate
               best_loadings = trial_loadings
               improved = .true.
            else
               uniqueness(i) = old_value
            end if
         end do
         if (.not. improved) step = 0.5_dp * step
         if (step <= tol) exit
      end do
      loadings = best_loadings
      objective = best
      info = 0
   end subroutine optimize_uniqueness

   subroutine factor_from_uniqueness(r, uniqueness, n_factors, loadings, objective, weights, info)
      real(dp), intent(in) :: r(:, :) !! Correlation matrix.
      real(dp), intent(in) :: uniqueness(:) !! Current uniqueness vector.
      integer, intent(in) :: n_factors !! Number of factors.
      real(dp), allocatable, intent(out) :: loadings(:, :) !! Loading matrix from the leading eigensystem.
      real(dp), intent(out) :: objective !! Weighted off-diagonal residual sum of squares.
      real(dp), intent(in) :: weights(:) !! Per-variable residual weights.
      integer, intent(out) :: info !! Eigensolver status.
      real(dp), allocatable :: work(:, :), values(:), vectors(:, :), residual(:, :)
      integer :: i, j, k, p

      p = size(r, 1)
      k = min(max(1, n_factors), p)
      allocate(work(p, p))
      work = r
      do i = 1, p
         work(i, i) = 1.0_dp - uniqueness(i)
      end do
      call symmetric_eigen(work, values, vectors, info, descending=.true.)
      allocate(loadings(p, k))
      loadings = 0.0_dp
      if (info /= 0) then
         objective = huge(1.0_dp)
         return
      end if
      do j = 1, k
         loadings(:, j) = vectors(:, j) * sqrt(max(0.0_dp, values(j)))
      end do
      allocate(residual(p, p))
      residual = r - matmul(loadings, transpose(loadings))
      objective = 0.0_dp
      do j = 1, p
         do i = j + 1, p
            objective = objective + weights(i) * weights(j) * residual(i, j) * residual(i, j)
         end do
      end do
   end subroutine factor_from_uniqueness

   subroutine residual_weights(r, weighted, weights, info)
      real(dp), intent(in) :: r(:, :) !! Correlation matrix.
      logical, intent(in) :: weighted !! Request inverse-correlation diagonal weights when true.
      real(dp), intent(out) :: weights(:) !! Per-variable WLS weights.
      integer, intent(out) :: info !! Zero when requested inversion succeeds.
      real(dp), allocatable :: rinv(:, :)
      integer :: i

      weights = 1.0_dp
      info = 0
      if (.not. weighted) return
      call pseudoinverse(r, rinv, info=info)
      if (info /= 0) return
      do i = 1, size(weights)
         if (abs(rinv(i, i)) > tiny(1.0_dp)) weights(i) = 1.0_dp / abs(rinv(i, i))
      end do
   end subroutine residual_weights

   subroutine orthogonal_cf_rotate(loadings, kappa, minimize_varimax, rotated, rotation, info, max_iter, tolerance)
      real(dp), intent(in) :: loadings(:, :) !! Loading matrix.
      real(dp), intent(in) :: kappa !! Crawford-Ferguson row/column weighting parameter.
      logical, intent(in) :: minimize_varimax !! Minimize centered-square variance rather than the CF criterion when true.
      real(dp), allocatable, intent(out) :: rotated(:, :) !! Rotated loadings.
      real(dp), allocatable, intent(out) :: rotation(:, :) !! Orthogonal rotation matrix.
      integer, intent(out) :: info !! Zero after deterministic optimization.
      integer, intent(in), optional :: max_iter !! Maximum pair sweeps.
      real(dp), intent(in), optional :: tolerance !! Minimum useful pair angle.
      real(dp) :: angle, best_angle, best_value, candidate, span, tol
      real(dp) :: c, s
      integer :: i, j, iter, grid, k, limit
      logical :: changed

      k = size(loadings, 2)
      allocate(rotated(size(loadings, 1), k), rotation(k, k))
      rotated = loadings
      rotation = 0.0_dp
      do i = 1, k
         rotation(i, i) = 1.0_dp
      end do
      limit = 60
      if (present(max_iter)) limit = max(1, max_iter)
      tol = 1.0e-7_dp
      if (present(tolerance)) tol = max(1.0e-12_dp, tolerance)
      do iter = 1, limit
         changed = .false.
         do i = 1, k - 1
            do j = i + 1, k
               best_angle = 0.0_dp
               if (minimize_varimax) then
                  best_value = varimin_value(rotated)
               else
                  best_value = cf_value(rotated, kappa)
               end if
               span = 0.78539816339744830962_dp
               do grid = -16, 16
                  angle = span * real(grid, dp) / 16.0_dp
                  candidate = pair_objective(rotated, i, j, angle, kappa, minimize_varimax)
                  if (candidate < best_value) then
                     best_value = candidate
                     best_angle = angle
                  end if
               end do
               if (abs(best_angle) > tol) then
                  c = cos(best_angle)
                  s = sin(best_angle)
                  call apply_pair_rotation(rotated, i, j, c, s)
                  call apply_pair_rotation(rotation, i, j, c, s)
                  changed = .true.
               end if
            end do
         end do
         if (.not. changed) exit
      end do
      info = 0
   end subroutine orthogonal_cf_rotate

   subroutine orthogonal_bimin_rotate(loadings, rotated, rotation, info, max_iter, tolerance)
      real(dp), intent(in) :: loadings(:, :) !! Loading matrix.
      real(dp), allocatable, intent(out) :: rotated(:, :) !! Bifactor-rotated loadings.
      real(dp), allocatable, intent(out) :: rotation(:, :) !! Orthogonal transformation.
      integer, intent(out) :: info !! Zero after deterministic optimization.
      integer, intent(in), optional :: max_iter !! Maximum pair sweeps.
      real(dp), intent(in), optional :: tolerance !! Minimum useful angle.
      real(dp) :: angle, best_angle, best_value, candidate, c, s, tol
      integer :: i, j, iter, grid, k, limit
      logical :: changed

      k = size(loadings, 2)
      allocate(rotated(size(loadings, 1), k), rotation(k, k))
      rotated = loadings
      rotation = 0.0_dp
      do i = 1, k
         rotation(i, i) = 1.0_dp
      end do
      limit = 60
      if (present(max_iter)) limit = max(1, max_iter)
      tol = 1.0e-7_dp
      if (present(tolerance)) tol = max(1.0e-12_dp, tolerance)
      do iter = 1, limit
         changed = .false.
         do i = 1, k - 1
            do j = i + 1, k
               best_angle = 0.0_dp
               best_value = bimin_value(rotated)
               do grid = -16, 16
                  angle = 0.78539816339744830962_dp * real(grid, dp) / 16.0_dp
                  candidate = bimin_pair_value(rotated, i, j, angle)
                  if (candidate < best_value) then
                     best_value = candidate
                     best_angle = angle
                  end if
               end do
               if (abs(best_angle) > tol) then
                  c = cos(best_angle)
                  s = sin(best_angle)
                  call apply_pair_rotation(rotated, i, j, c, s)
                  call apply_pair_rotation(rotation, i, j, c, s)
                  changed = .true.
               end if
            end do
         end do
         if (.not. changed) exit
      end do
      info = 0
   end subroutine orthogonal_bimin_rotate

   subroutine oblique_bimin_coordinate_search(loadings, rotated, transformation, phi, info, max_iter, tolerance)
      real(dp), intent(in) :: loadings(:, :) !! Loading matrix.
      real(dp), allocatable, intent(out) :: rotated(:, :) !! Obliquely rotated loadings.
      real(dp), allocatable, intent(out) :: transformation(:, :) !! Oblique transformation.
      real(dp), allocatable, intent(out) :: phi(:, :) !! Factor-correlation matrix.
      integer, intent(out) :: info !! Zero when the transformation pseudoinverse succeeds.
      integer, intent(in), optional :: max_iter !! Maximum coordinate-search sweeps.
      real(dp), intent(in), optional :: tolerance !! Minimum coordinate step.
      real(dp), allocatable :: trial(:, :), inv(:, :)
      real(dp) :: best, candidate, old_value, step, tol
      integer :: i, j, iter, k, limit
      logical :: improved

      k = size(loadings, 2)
      allocate(transformation(k, k))
      transformation = 0.0_dp
      do i = 1, k
         transformation(i, i) = 1.0_dp
      end do
      best = bimin_value(matmul(loadings, transformation))
      step = 0.10_dp
      limit = 80
      if (present(max_iter)) limit = max(1, max_iter)
      tol = 1.0e-5_dp
      if (present(tolerance)) tol = max(1.0e-10_dp, tolerance)
      do iter = 1, limit
         improved = .false.
         do j = 1, k
            do i = 1, k
               if (i == j) cycle
               old_value = transformation(i, j)
               transformation(i, j) = old_value + step
               call normalize_transform_columns(transformation)
               trial = matmul(loadings, transformation)
               candidate = bimin_value(trial)
               if (candidate < best) then
                  best = candidate
                  improved = .true.
                  cycle
               end if
               transformation(i, j) = old_value - step
               call normalize_transform_columns(transformation)
               trial = matmul(loadings, transformation)
               candidate = bimin_value(trial)
               if (candidate < best) then
                  best = candidate
                  improved = .true.
               else
                  transformation(i, j) = old_value
                  call normalize_transform_columns(transformation)
               end if
            end do
         end do
         if (.not. improved) step = 0.5_dp * step
         if (step <= tol) exit
      end do
      rotated = matmul(loadings, transformation)
      call pseudoinverse(transformation, inv, info=info)
      allocate(phi(k, k))
      if (info == 0) then
         phi = matmul(inv, transpose(inv))
         call standardize_symmetric(phi)
      else
         phi = 0.0_dp
      end if
   end subroutine oblique_bimin_coordinate_search

   pure real(dp) function pair_objective(loadings, i, j, angle, kappa, minimize_varimax) result(value)
      real(dp), intent(in) :: loadings(:, :) !! Current loading matrix.
      integer, intent(in) :: i !! First factor column.
      integer, intent(in) :: j !! Second factor column.
      real(dp), intent(in) :: angle !! Trial Givens angle in radians.
      real(dp), intent(in) :: kappa !! Crawford-Ferguson criterion parameter.
      logical, intent(in) :: minimize_varimax !! Select the centered-square criterion when true.
      real(dp), allocatable :: trial(:, :)

      trial = loadings
      call apply_pair_rotation(trial, i, j, cos(angle), sin(angle))
      if (minimize_varimax) then
         value = varimin_value(trial)
      else
         value = cf_value(trial, kappa)
      end if
   end function pair_objective

   pure real(dp) function bimin_pair_value(loadings, i, j, angle) result(value)
      real(dp), intent(in) :: loadings(:, :) !! Current loading matrix.
      integer, intent(in) :: i !! First factor column.
      integer, intent(in) :: j !! Second factor column.
      real(dp), intent(in) :: angle !! Trial Givens angle in radians.
      real(dp), allocatable :: trial(:, :)

      trial = loadings
      call apply_pair_rotation(trial, i, j, cos(angle), sin(angle))
      value = bimin_value(trial)
   end function bimin_pair_value

   pure real(dp) function cf_value(loadings, kappa) result(value)
      real(dp), intent(in) :: loadings(:, :) !! Loading matrix.
      real(dp), intent(in) :: kappa !! Crawford-Ferguson criterion parameter.
      real(dp) :: row_term, col_term
      integer :: i, j, h, k

      row_term = 0.0_dp
      col_term = 0.0_dp
      do i = 1, size(loadings, 1)
         do j = 1, size(loadings, 2) - 1
            do k = j + 1, size(loadings, 2)
               row_term = row_term + loadings(i, j) ** 2 * loadings(i, k) ** 2
            end do
         end do
      end do
      do j = 1, size(loadings, 2)
         do i = 1, size(loadings, 1) - 1
            do h = i + 1, size(loadings, 1)
               col_term = col_term + loadings(i, j) ** 2 * loadings(h, j) ** 2
            end do
         end do
      end do
      value = (1.0_dp - kappa) * row_term + kappa * col_term
   end function cf_value

   pure real(dp) function varimin_value(loadings) result(value)
      real(dp), intent(in) :: loadings(:, :) !! Loading matrix.
      real(dp), allocatable :: gradient(:, :)
      allocate(gradient(size(loadings, 1), size(loadings, 2)))
      call varimin_objective_gradient(loadings, value, gradient)
   end function varimin_value

   pure real(dp) function bimin_value(loadings) result(value)
      real(dp), intent(in) :: loadings(:, :) !! Loading matrix.
      real(dp), allocatable :: gradient(:, :)
      allocate(gradient(size(loadings, 1), size(loadings, 2)))
      call bimin_objective_gradient(loadings, value, gradient)
   end function bimin_value

   pure subroutine apply_pair_rotation(a, i, j, c, s)
      real(dp), intent(inout) :: a(:, :) !! Matrix whose selected columns are rotated in place.
      integer, intent(in) :: i !! First selected column.
      integer, intent(in) :: j !! Second selected column.
      real(dp), intent(in) :: c !! Cosine of the Givens angle.
      real(dp), intent(in) :: s !! Sine of the Givens angle.
      real(dp), allocatable :: first(:), second(:)

      allocate(first(size(a, 1)), second(size(a, 1)))
      first = a(:, i)
      second = a(:, j)
      a(:, i) = c * first + s * second
      a(:, j) = -s * first + c * second
   end subroutine apply_pair_rotation

   pure subroutine normalize_transform_columns(t)
      real(dp), intent(inout) :: t(:, :) !! Oblique transformation normalized to unit-length columns.
      real(dp) :: norm
      integer :: j

      do j = 1, size(t, 2)
         norm = sqrt(sum(t(:, j) * t(:, j)))
         if (norm > tiny(1.0_dp)) t(:, j) = t(:, j) / norm
      end do
   end subroutine normalize_transform_columns

   pure subroutine standardize_symmetric(a)
      real(dp), intent(inout) :: a(:, :) !! Positive-diagonal symmetric matrix standardized to unit diagonal.
      real(dp), allocatable :: sd(:)
      integer :: i, j

      allocate(sd(size(a, 1)))
      do i = 1, size(a, 1)
         sd(i) = sqrt(max(tiny(1.0_dp), a(i, i)))
      end do
      do j = 1, size(a, 2)
         do i = 1, size(a, 1)
            a(i, j) = a(i, j) / (sd(i) * sd(j))
         end do
      end do
   end subroutine standardize_symmetric

end module psych_factor_extended
