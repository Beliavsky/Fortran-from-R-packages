module psych_validity_extended
   use r_kinds, only : dp
   use psych_basic, only : yule_inverse, phi_coefficient
   use psych_distributions, only : normal_quantile, normal_cdf
   use psych_polychoric, only : tetrachoric_table, phi_to_tetrachoric
   use psych_stats, only : correlation_matrix
   use psych_correlations, only : pseudoinverse
   use psych_factor, only : principal_axis_factor, principal_components_from_matrix
   use psych_types, only : pca_result
   implicit none
   private
   public :: item_validity_by_keys, predicted_validity_summary, effect_validity_by_keys
   public :: yule_bonett_interval, yule_correlation_matrix, yule_to_phi_values
   public :: yule_to_tetrachoric_value, yule_to_tetrachoric_values
   public :: phi_to_tetrachoric_values, regression_from_matrices, set_correlation_value
   public :: unidimensionality_statistics, vss_complexity_fit
   public :: cross_validation_from_correlation

contains

   pure subroutine cross_validation_from_correlation(weights, correlation, cross_validity, matrix_r, info)
      real(dp), intent(in) :: weights(:, :) !! Predictor-by-criterion weight matrix used to form criterion score estimates.
      real(dp), intent(in) :: correlation(:, :) !! Correlation matrix ordered as criteria first, followed by predictors.
      real(dp), allocatable, intent(out) :: cross_validity(:) !! Correlation of each weighted predictor score with its criterion.
      real(dp), allocatable, intent(out) :: matrix_r(:, :) !! Weighted-score correlations with every supplied criterion.
      integer, intent(out), optional :: info !! Zero for compatible dimensions; one when the supplied matrix is too small.
      real(dp), allocatable :: rxx(:, :), w(:)
      real(dp) :: variance_score, denom
      integer :: i, j, ncrit, npred

      npred = size(weights, 1)
      ncrit = size(weights, 2)
      allocate(cross_validity(ncrit), matrix_r(ncrit, ncrit))
      cross_validity = 0.0_dp
      matrix_r = 0.0_dp
      if (size(correlation, 1) < ncrit + npred .or. size(correlation, 2) < ncrit + npred) then
         if (present(info)) info = 1
         return
      end if
      if (present(info)) info = 0

      allocate(rxx(npred, npred), w(npred))
      rxx = correlation(ncrit + 1:ncrit + npred, ncrit + 1:ncrit + npred)
      do i = 1, ncrit
         w = weights(:, i)
         variance_score = dot_product(w, matmul(rxx, w))
         denom = sqrt(max(variance_score, tiny(1.0_dp)))
         do j = 1, ncrit
            matrix_r(i, j) = dot_product(w, correlation(j, ncrit + 1:ncrit + npred)) / denom
         end do
         cross_validity(i) = matrix_r(i, i)
      end do
   end subroutine cross_validation_from_correlation


   pure subroutine item_validity_by_keys(items, criteria, keys, validity)
      real(dp), intent(in) :: items(:, :) !! Subject-by-item data matrix used to form item-criterion correlations.
      real(dp), intent(in) :: criteria(:, :) !! Subject-by-criterion matrix aligned with `items`.
      real(dp), intent(in) :: keys(:, :) !! Item-by-scale signed key matrix; zero entries exclude an item.
      real(dp), allocatable, intent(out) :: validity(:, :) !! Scale-by-criterion mean signed item validities.
      real(dp), allocatable :: combined(:, :), correlations(:, :)
      real(dp) :: total
      integer :: i, j, k, n_used, n_item, n_criterion

      n_item = size(items, 2)
      n_criterion = size(criteria, 2)
      allocate(validity(size(keys, 2), n_criterion))
      validity = 0.0_dp
      if (size(items, 1) /= size(criteria, 1) .or. size(keys, 1) /= n_item) return
      allocate(combined(size(items, 1), n_item + n_criterion))
      allocate(correlations(n_item + n_criterion, n_item + n_criterion))
      combined(:, 1:n_item) = items
      combined(:, n_item + 1:) = criteria
      call correlation_matrix(combined, correlations)
      do k = 1, size(keys, 2)
         do j = 1, n_criterion
            total = 0.0_dp
            n_used = 0
            do i = 1, n_item
               if (abs(keys(i, k)) > tiny(1.0_dp)) then
                  total = total + sign(1.0_dp, keys(i, k)) * correlations(i, n_item + j)
                  n_used = n_used + 1
               end if
            end do
            if (n_used > 0) validity(k, j) = total / real(n_used, dp)
         end do
      end do
   end subroutine item_validity_by_keys

   pure subroutine predicted_validity_summary(item_validity, mean_r, n_items, predicted, asymptotic)
      real(dp), intent(in) :: item_validity(:, :) !! Scale-by-criterion mean item validities.
      real(dp), intent(in) :: mean_r(:) !! Mean inter-item correlation for each scale.
      integer, intent(in) :: n_items(:) !! Number of keyed items in each scale.
      real(dp), allocatable, intent(out) :: predicted(:, :) !! Spearman-Brown-style predicted scale validity.
      real(dp), allocatable, intent(out) :: asymptotic(:, :) !! Infinite-length limiting validity estimate.
      real(dp) :: denom
      integer :: i

      allocate(predicted(size(item_validity, 1), size(item_validity, 2)))
      allocate(asymptotic(size(item_validity, 1), size(item_validity, 2)))
      predicted = 0.0_dp
      asymptotic = 0.0_dp
      do i = 1, min(size(item_validity, 1), min(size(mean_r), size(n_items)))
         if (n_items(i) > 0) then
            denom = real(n_items(i), dp) + real(n_items(i) * max(0, n_items(i) - 1), dp) * mean_r(i)
            if (denom > tiny(1.0_dp)) then
               predicted(i, :) = item_validity(i, :) * real(n_items(i), dp) / sqrt(denom)
            end if
         end if
         if (mean_r(i) > tiny(1.0_dp)) asymptotic(i, :) = item_validity(i, :) / sqrt(mean_r(i))
      end do
   end subroutine predicted_validity_summary

   pure subroutine effect_validity_by_keys(effect, keys, validity)
      real(dp), intent(in) :: effect(:) !! Item effect sizes in the same order as rows of `keys`.
      real(dp), intent(in) :: keys(:, :) !! Item-by-scale signed key matrix.
      real(dp), allocatable, intent(out) :: validity(:) !! Mean signed effect size for each scale.
      real(dp) :: total
      integer :: i, k, n_used

      allocate(validity(size(keys, 2)))
      validity = 0.0_dp
      if (size(keys, 1) /= size(effect)) return
      do k = 1, size(keys, 2)
         total = 0.0_dp
         n_used = 0
         do i = 1, size(effect)
            if (abs(keys(i, k)) > tiny(1.0_dp)) then
               total = total + sign(1.0_dp, keys(i, k)) * effect(i)
               n_used = n_used + 1
            end if
         end do
         if (n_used > 0) validity(k) = total / real(n_used, dp)
      end do
   end subroutine effect_validity_by_keys

   pure subroutine yule_bonett_interval(table2x2, exponent, use_bonett_exponent, alpha, rho, se, lower, upper)
      real(dp), intent(in) :: table2x2(2, 2) !! Two-by-two nonnegative frequency table.
      real(dp), intent(in) :: exponent !! Generalized-Yule exponent used when Bonett adaptation is disabled.
      logical, intent(in) :: use_bonett_exponent !! Derive the exponent from the smallest marginal when true.
      real(dp), intent(in) :: alpha !! Two-sided tail probability, normally 0.05.
      real(dp), intent(out) :: rho !! Small-cell-adjusted generalized Yule association.
      real(dp), intent(out) :: se !! Bonett large-sample standard error.
      real(dp), intent(out) :: lower !! Lower confidence limit.
      real(dp), intent(out) :: upper !! Upper confidence limit.
      real(dp) :: p(2, 2), total, c_value, ad, bc, odds, variance_log, critical
      real(dp) :: row_sum(2), col_sum(2)

      total = sum(max(table2x2, 0.0_dp))
      if (total <= tiny(1.0_dp)) then
         rho = 0.0_dp
         se = 0.0_dp
         lower = 0.0_dp
         upper = 0.0_dp
         return
      end if
      p = max(table2x2, 0.0_dp) / total
      row_sum = sum(p, dim=2)
      col_sum = sum(p, dim=1)
      c_value = exponent
      if (use_bonett_exponent) c_value = 0.5_dp - (0.5_dp - min(minval(row_sum), minval(col_sum)))**2
      ad = (p(1, 1) * p(2, 2))**c_value
      bc = (p(1, 2) * p(2, 1))**c_value
      if (ad + bc <= tiny(1.0_dp)) then
         rho = 0.0_dp
      else
         odds = ((table2x2(1, 1) + 0.1_dp) * (table2x2(2, 2) + 0.1_dp)) / &
            ((table2x2(2, 1) + 0.1_dp) * (table2x2(1, 2) + 0.1_dp))
         rho = (odds**c_value - 1.0_dp) / (odds**c_value + 1.0_dp)
      end if
      variance_log = 1.0_dp / (table2x2(1, 1) + 0.1_dp) + 1.0_dp / (table2x2(2, 2) + 0.1_dp) + &
         1.0_dp / (table2x2(2, 1) + 0.1_dp) + 1.0_dp / (table2x2(1, 2) + 0.1_dp)
      se = sqrt(max(0.0_dp, (c_value * c_value / 4.0_dp) * (1.0_dp - rho * rho)**2 * variance_log))
      critical = normal_quantile(1.0_dp - 0.5_dp * alpha)
      lower = tanh(atanh(max(-1.0_dp + 1.0e-12_dp, min(1.0_dp - 1.0e-12_dp, rho))) - &
         critical * sqrt((c_value * c_value / 4.0_dp) * variance_log))
      upper = tanh(atanh(max(-1.0_dp + 1.0e-12_dp, min(1.0_dp - 1.0e-12_dp, rho))) + &
         critical * sqrt((c_value * c_value / 4.0_dp) * variance_log))
   end subroutine yule_bonett_interval

   pure subroutine yule_correlation_matrix(items, exponent, use_bonett_exponent, alpha, rho, ci, z_p)
      integer, intent(in) :: items(:, :) !! Subject-by-binary-item matrix with values 0 or 1.
      real(dp), intent(in) :: exponent !! Generalized-Yule exponent used for every pair.
      logical, intent(in) :: use_bonett_exponent !! Use pair-specific Bonett exponent when true.
      real(dp), intent(in) :: alpha !! Confidence-interval tail probability.
      real(dp), allocatable, intent(out) :: rho(:, :) !! Symmetric generalized-Yule correlation matrix.
      real(dp), allocatable, intent(out) :: ci(:, :) !! Lower triangle lower limits and upper triangle upper limits.
      real(dp), allocatable, intent(out) :: z_p(:, :) !! Lower triangle z statistics and upper triangle normal-tail p-values.
      real(dp) :: table2x2(2, 2), pair_rho, se, lower, upper, z
      integer :: i, j, s, a, b, p

      p = size(items, 2)
      allocate(rho(p, p), ci(p, p), z_p(p, p))
      rho = 0.0_dp
      ci = 0.0_dp
      z_p = 0.0_dp
      do i = 1, p
         do j = 1, i
            table2x2 = 0.0_dp
            do s = 1, size(items, 1)
               a = min(1, max(0, items(s, i))) + 1
               b = min(1, max(0, items(s, j))) + 1
               table2x2(a, b) = table2x2(a, b) + 1.0_dp
            end do
            call yule_bonett_interval(table2x2, exponent, use_bonett_exponent, alpha, pair_rho, se, lower, upper)
            rho(i, j) = pair_rho
            rho(j, i) = pair_rho
            ci(i, j) = lower
            ci(j, i) = upper
            if (se > tiny(1.0_dp)) then
               z = pair_rho / se
            else
               z = 0.0_dp
            end if
            z_p(i, j) = z
            z_p(j, i) = 1.0_dp - normal_cdf(abs(z))
         end do
      end do
   end subroutine yule_correlation_matrix

   pure subroutine yule_to_phi_values(q, margins, phi)
      real(dp), intent(in) :: q(:, :) !! Symmetric matrix of Yule-Q coefficients.
      real(dp), intent(in) :: margins(:) !! First-category marginal proportions for each variable.
      real(dp), allocatable, intent(out) :: phi(:, :) !! Phi correlations implied by Yule Q and margins.
      real(dp) :: table2x2(2, 2)
      integer :: i, j, p

      p = min(size(q, 1), min(size(q, 2), size(margins)))
      allocate(phi(p, p))
      phi = 0.0_dp
      do i = 1, p
         phi(i, i) = 1.0_dp
         do j = 1, i - 1
            call yule_inverse(q(i, j), [margins(i), margins(j)], table2x2)
            phi(i, j) = phi_coefficient(table2x2)
            phi(j, i) = phi(i, j)
         end do
      end do
   end subroutine yule_to_phi_values

   pure subroutine yule_to_tetrachoric_value(q, margins, rho)
      real(dp), intent(in) :: q !! Yule-Q coefficient in `[-1,1]`.
      real(dp), intent(in) :: margins(2) !! First-category marginal proportions of the two binary variables.
      real(dp), intent(out) :: rho !! Tetrachoric correlation implied by Q and the supplied margins.
      real(dp) :: table2x2(2, 2)

      call yule_inverse(q, margins, table2x2)
      call tetrachoric_table(table2x2, rho, 0.0_dp)
   end subroutine yule_to_tetrachoric_value

   pure subroutine yule_to_tetrachoric_values(q, margins, tetrachoric)
      real(dp), intent(in) :: q(:, :) !! Symmetric matrix of Yule-Q coefficients.
      real(dp), intent(in) :: margins(:) !! First-category marginal proportions for each variable.
      real(dp), allocatable, intent(out) :: tetrachoric(:, :) !! Tetrachoric correlations implied by Q and margins.
      real(dp) :: table2x2(2, 2), rho
      integer :: i, j, p

      p = min(size(q, 1), min(size(q, 2), size(margins)))
      allocate(tetrachoric(p, p))
      tetrachoric = 0.0_dp
      do i = 1, p
         tetrachoric(i, i) = 1.0_dp
         do j = 1, i - 1
            call yule_inverse(q(i, j), [margins(i), margins(j)], table2x2)
            call tetrachoric_table(table2x2, rho, 0.0_dp)
            tetrachoric(i, j) = rho
            tetrachoric(j, i) = rho
         end do
      end do
   end subroutine yule_to_tetrachoric_values

   pure subroutine phi_to_tetrachoric_values(phi, margins, tetrachoric)
      real(dp), intent(in) :: phi(:, :) !! Symmetric matrix of observed binary phi correlations.
      real(dp), intent(in) :: margins(:) !! First-category marginal proportions for each variable.
      real(dp), allocatable, intent(out) :: tetrachoric(:, :) !! Tetrachoric correlations implied by phi and margins.
      real(dp) :: rho
      integer :: i, j, p

      p = min(size(phi, 1), min(size(phi, 2), size(margins)))
      allocate(tetrachoric(p, p))
      tetrachoric = 0.0_dp
      do i = 1, p
         tetrachoric(i, i) = 1.0_dp
         do j = 1, i - 1
            call phi_to_tetrachoric(phi(i, j), margins(i), margins(j), rho)
            tetrachoric(i, j) = rho
            tetrachoric(j, i) = rho
         end do
      end do
   end subroutine phi_to_tetrachoric_values

   subroutine regression_from_matrices(x_matrix, y_matrix, xy_matrix, beta, r2, residual, set_r2)
      real(dp), intent(in) :: x_matrix(:, :) !! Predictor covariance or correlation matrix.
      real(dp), intent(in) :: y_matrix(:, :) !! Criterion covariance or correlation matrix.
      real(dp), intent(in) :: xy_matrix(:, :) !! Predictor-by-criterion covariance or correlation block.
      real(dp), allocatable, intent(out) :: beta(:, :) !! Regression weights solving X beta = XY.
      real(dp), allocatable, intent(out) :: r2(:) !! Multiple R-squared for each criterion.
      real(dp), allocatable, intent(out) :: residual(:, :) !! Residual criterion covariance matrix.
      real(dp), intent(out) :: set_r2 !! Cohen set correlation based on determinant reduction.
      real(dp), allocatable :: x_inverse(:, :), joint(:, :)
      real(dp) :: dx, dy, dj
      integer :: i, nx, ny, info

      nx = size(x_matrix, 1)
      ny = size(y_matrix, 1)
      allocate(beta(nx, ny), r2(ny), residual(ny, ny))
      call pseudoinverse(x_matrix, x_inverse, info=info)
      if (info /= 0) x_inverse = 0.0_dp
      beta = matmul(x_inverse, xy_matrix)
      residual = y_matrix - matmul(transpose(xy_matrix), beta)
      do i = 1, ny
         if (abs(y_matrix(i, i)) > tiny(1.0_dp)) then
            r2(i) = 1.0_dp - residual(i, i) / y_matrix(i, i)
         else
            r2(i) = 0.0_dp
         end if
      end do
      allocate(joint(nx + ny, nx + ny))
      joint = 0.0_dp
      joint(1:nx, 1:nx) = x_matrix
      joint(1:nx, nx + 1:) = xy_matrix
      joint(nx + 1:, 1:nx) = transpose(xy_matrix)
      joint(nx + 1:, nx + 1:) = y_matrix
      dx = determinant_local(x_matrix)
      dy = determinant_local(y_matrix)
      dj = determinant_local(joint)
      if (abs(dx * dy) > tiny(1.0_dp)) then
         set_r2 = 1.0_dp - dj / (dx * dy)
      else
         set_r2 = 0.0_dp
      end if
   end subroutine regression_from_matrices

   pure real(dp) function set_correlation_value(x_matrix, y_matrix, xy_matrix) result(value)
      real(dp), intent(in) :: x_matrix(:, :) !! Predictor covariance or correlation block.
      real(dp), intent(in) :: y_matrix(:, :) !! Criterion covariance or correlation block.
      real(dp), intent(in) :: xy_matrix(:, :) !! Predictor-by-criterion cross block.
      real(dp), allocatable :: joint(:, :)
      real(dp) :: dx, dy
      integer :: nx, ny

      nx = size(x_matrix, 1)
      ny = size(y_matrix, 1)
      allocate(joint(nx + ny, nx + ny))
      joint = 0.0_dp
      joint(1:nx, 1:nx) = x_matrix
      joint(1:nx, nx + 1:) = xy_matrix
      joint(nx + 1:, 1:nx) = transpose(xy_matrix)
      joint(nx + 1:, nx + 1:) = y_matrix
      dx = determinant_local(x_matrix)
      dy = determinant_local(y_matrix)
      if (abs(dx * dy) > tiny(1.0_dp)) then
         value = 1.0_dp - determinant_local(joint) / (dx * dy)
      else
         value = 0.0_dp
      end if
   end function set_correlation_value

   subroutine unidimensionality_statistics(r, n_factors, u, alpha_std, average_r, median_r, map_value, f1_f2)
      real(dp), intent(in) :: r(:, :) !! Correlation matrix for one candidate scale.
      integer, intent(in) :: n_factors !! Number of factors used for the first-to-second eigenvalue diagnostic.
      real(dp), intent(out) :: u !! Unidimensionality index combining one-factor fit and alpha-model residual fit.
      real(dp), intent(out) :: alpha_std !! Standardized Cronbach alpha after one-factor sign reflection.
      real(dp), intent(out) :: average_r !! Mean off-diagonal correlation after reflection.
      real(dp), intent(out) :: median_r !! Median lower-triangle correlation after reflection.
      real(dp), intent(out) :: map_value !! Velicer-style one-component mean squared partial residual.
      real(dp), intent(out) :: f1_f2 !! Ratio of first to second correlation-matrix eigenvalues.
      real(dp), allocatable :: loadings(:, :), uniqueness(:), model(:, :), reflected(:, :)
      real(dp), allocatable :: residual(:,:), diag_scale(:), offdiag(:)
      real(dp) :: fit, alpha_res, total_sq, residual_sq, total, trace_r
      type(pca_result) :: pca
      integer :: i, j, n, k, info

      n = size(r, 1)
      u = 0.0_dp
      alpha_std = 0.0_dp
      average_r = 0.0_dp
      median_r = 0.0_dp
      map_value = 0.0_dp
      f1_f2 = 0.0_dp
      if (n < 2 .or. size(r, 2) /= n) return
      allocate(loadings(n, 1), uniqueness(n))
      call principal_axis_factor(r, 1, loadings, uniqueness, info, max_iter=100, tolerance=1.0e-8_dp)
      allocate(model(n, n), reflected(n, n))
      model = matmul(loadings, transpose(loadings))
      reflected = r
      do i = 1, n
         if (loadings(i, 1) < 0.0_dp) then
            reflected(i, :) = -reflected(i, :)
            reflected(:, i) = -reflected(:, i)
         end if
      end do
      trace_r = sum([(reflected(i, i), i = 1, n)])
      total = sum(reflected)
      if (n > 1 .and. abs(total) > tiny(1.0_dp)) alpha_std = (1.0_dp - trace_r / total) * real(n, dp) / real(n - 1, dp)
      average_r = (total - trace_r) / real(n * (n - 1), dp)
      allocate(offdiag(n * (n - 1) / 2))
      k = 0
      do i = 2, n
         do j = 1, i - 1
            k = k + 1
            offdiag(k) = reflected(i, j)
         end do
      end do
      call sort_real_in_place(offdiag)
      if (size(offdiag) > 0) then
         if (mod(size(offdiag), 2) == 1) then
            median_r = offdiag((size(offdiag) + 1) / 2)
         else
            median_r = 0.5_dp * (offdiag(size(offdiag) / 2) + offdiag(size(offdiag) / 2 + 1))
         end if
      end if
      residual_sq = 0.0_dp
      total_sq = 0.0_dp
      do i = 2, n
         do j = 1, i - 1
            residual_sq = residual_sq + (reflected(i, j) - average_r)**2
            total_sq = total_sq + reflected(i, j)**2
         end do
      end do
      if (total_sq > tiny(1.0_dp)) then
         alpha_res = residual_sq / total_sq
      else
         alpha_res = 1.0_dp
      end if
      fit = 1.0_dp
      total_sq = 0.0_dp
      residual_sq = 0.0_dp
      do i = 2, n
         do j = 1, i - 1
            total_sq = total_sq + r(i, j)**2
            residual_sq = residual_sq + (r(i, j) - model(i, j))**2
         end do
      end do
      if (total_sq > tiny(1.0_dp)) fit = 1.0_dp - residual_sq / total_sq
      u = fit * (1.0_dp - alpha_res)
      call principal_components_from_matrix(r, min(max(2, n_factors), n), pca)
      if (size(pca%values) >= 2) then
         if (abs(pca%values(2)) > tiny(1.0_dp)) f1_f2 = pca%values(1) / pca%values(2)
      end if
      allocate(residual(n, n), diag_scale(n))
      residual = r - matmul(pca%loadings(:, 1:1), transpose(pca%loadings(:, 1:1)))
      do i = 1, n
         diag_scale(i) = 1.0_dp / sqrt(max(residual(i, i), 1.0e-12_dp))
      end do
      map_value = 0.0_dp
      do i = 1, n
         do j = 1, n
            if (i /= j) map_value = map_value + (diag_scale(i) * residual(i, j) * diag_scale(j))**2
         end do
      end do
      map_value = map_value / real(n * (n - 1), dp)
   end subroutine unidimensionality_statistics

   pure subroutine vss_complexity_fit(r, loadings, phi, max_complexity, fit, residual_sum)
      real(dp), intent(in) :: r(:, :) !! Target correlation matrix.
      real(dp), intent(in) :: loadings(:, :) !! Unrotated or rotated factor loading matrix.
      real(dp), intent(in) :: phi(:, :) !! Factor-correlation matrix; identity for orthogonal factors.
      integer, intent(in) :: max_complexity !! Maximum retained nonzero loadings per item.
      real(dp), allocatable, intent(out) :: fit(:) !! VSS fit for complexities one through `max_complexity`.
      real(dp), allocatable, intent(out) :: residual_sum(:) !! Off-diagonal residual sum of squares by complexity.
      real(dp), allocatable :: simple(:, :), model(:, :), work(:)
      real(dp) :: total_original
      integer :: c, i, j, k, idx, n_keep

      allocate(fit(max(0, max_complexity)), residual_sum(max(0, max_complexity)))
      fit = 0.0_dp
      residual_sum = 0.0_dp
      if (max_complexity <= 0) return
      total_original = 0.0_dp
      do i = 2, size(r, 1)
         do j = 1, min(i - 1, size(r, 2))
            total_original = total_original + 2.0_dp * r(i, j)**2
         end do
      end do
      allocate(simple(size(loadings, 1), size(loadings, 2)), work(size(loadings, 2)))
      allocate(model(size(r, 1), size(r, 2)))
      do c = 1, max_complexity
         simple = 0.0_dp
         do i = 1, size(loadings, 1)
            work = abs(loadings(i, :))
            n_keep = min(c, size(loadings, 2))
            do k = 1, n_keep
               idx = maxloc(work, dim=1)
               simple(i, idx) = loadings(i, idx)
               work(idx) = -1.0_dp
            end do
         end do
         model = matmul(matmul(simple, phi), transpose(simple))
         residual_sum(c) = 0.0_dp
         do i = 1, size(r, 1)
            do j = 1, size(r, 2)
               if (i /= j) residual_sum(c) = residual_sum(c) + (r(i, j) - model(i, j))**2
            end do
         end do
         if (total_original > tiny(1.0_dp)) fit(c) = 1.0_dp - residual_sum(c) / total_original
      end do
   end subroutine vss_complexity_fit

   pure real(dp) function determinant_local(a) result(value)
      real(dp), intent(in) :: a(:, :) !! Square matrix whose determinant is required.
      real(dp), allocatable :: work(:, :), row_tmp(:)
      real(dp) :: pivot, factor
      integer :: i, j, k, n, pivot_row

      n = size(a, 1)
      if (size(a, 2) /= n) then
         value = 0.0_dp
         return
      end if
      allocate(work(n, n), row_tmp(n))
      work = a
      value = 1.0_dp
      do k = 1, n
         pivot_row = k - 1 + maxloc(abs(work(k:n, k)), dim=1)
         if (abs(work(pivot_row, k)) <= tiny(1.0_dp)) then
            value = 0.0_dp
            return
         end if
         if (pivot_row /= k) then
            row_tmp = work(k, :)
            work(k, :) = work(pivot_row, :)
            work(pivot_row, :) = row_tmp
            value = -value
         end if
         pivot = work(k, k)
         value = value * pivot
         do i = k + 1, n
            factor = work(i, k) / pivot
            do j = k + 1, n
               work(i, j) = work(i, j) - factor * work(k, j)
            end do
         end do
      end do
   end function determinant_local

   pure subroutine sort_real_in_place(x)
      real(dp), intent(inout) :: x(:) !! Real vector sorted into ascending order in place.
      real(dp) :: value
      integer :: i, j

      do i = 2, size(x)
         value = x(i)
         j = i - 1
         do while (j >= 1)
            if (x(j) <= value) exit
            x(j + 1) = x(j)
            j = j - 1
         end do
         x(j + 1) = value
      end do
   end subroutine sort_real_in_place

end module psych_validity_extended
