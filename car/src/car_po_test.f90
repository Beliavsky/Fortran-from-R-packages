module car_po_test
   use car_kinds, only: dp
   use car_utils, only: chi_square_survival, matrix_inverse
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_is_nan, ieee_quiet_nan, ieee_value
   implicit none
   private

   public :: po_test

contains

   subroutine po_test(x, y_level, overall_chisq, overall_df, overall_p, predictor_chisq, predictor_p, &
                      predictor_df, status, weights, beta_threshold, max_iter, tol)
      real(dp), intent(in) :: x(:, :) !! Model matrix with an intercept in column one and predictors in remaining columns.
      integer, intent(in) :: y_level(:) !! Ordered response codes 1..K, with every level represented at least once.
      real(dp), intent(out) :: overall_chisq !! Wald chi-square statistic for equality of slopes across cumulative logits.
      integer, intent(out) :: overall_df !! Degrees of freedom, (K-2) times the number of non-intercept predictors.
      real(dp), intent(out) :: overall_p !! Upper-tail chi-square probability for overall_chisq.
      real(dp), allocatable, intent(out) :: predictor_chisq(:) !! Predictor-specific chi-square statistics in x-column order after the intercept.
      real(dp), allocatable, intent(out) :: predictor_p(:) !! Upper-tail chi-square probabilities for predictor_chisq.
      integer, intent(out) :: predictor_df !! Degrees of freedom K-2 for each predictor-specific test.
      integer, intent(out) :: status !! Zero on success; nonzero for invalid inputs, failed logistic fits, or singular covariance matrices.
      real(dp), intent(in), optional :: weights(:) !! Optional nonnegative case weights; omitted values are treated as one.
      real(dp), allocatable, intent(out), optional :: beta_threshold(:, :) !! Fitted slope matrix with shape (P, K-1), excluding intercepts.
      integer, intent(in), optional :: max_iter !! Maximum logistic IRLS iterations per threshold; default 100.
      real(dp), intent(in), optional :: tol !! Positive relative coefficient-convergence tolerance; default 1e-10.

      real(dp), allocatable :: beta_all(:, :)
      real(dp), allocatable :: beta_flat(:)
      real(dp), allocatable :: covariance(:, :)
      real(dp), allocatable :: dmat(:, :)
      real(dp), allocatable :: dvcov(:, :)
      real(dp), allocatable :: dvcov_inv(:, :)
      real(dp), allocatable :: fitted(:, :)
      real(dp), allocatable :: q(:)
      real(dp), allocatable :: v_jj(:, :)
      real(dp), allocatable :: v_el(:, :)
      real(dp), allocatable :: v_mid(:, :)
      real(dp), allocatable :: w_case(:)
      real(dp), allocatable :: w_cross(:)
      real(dp), allocatable :: w_diag(:)
      real(dp), allocatable :: binary(:)
      real(dp), allocatable :: dd(:, :)
      real(dp), allocatable :: q_pred(:)
      real(dp), allocatable :: pred_cov(:, :)
      real(dp), allocatable :: pred_cov_inv(:, :)
      integer :: el
      integer :: fit_status
      integer :: i
      integer :: inv_status
      integer :: j
      integer :: k
      integer :: limit
      integer :: n
      integer :: n_threshold
      integer :: p
      integer :: row0
      integer :: col0
      real(dp) :: tolerance

      n = size(x, 1)
      p = size(x, 2) - 1
      k = 0
      if (size(y_level) > 0) k = maxval(y_level)
      n_threshold = k - 1
      predictor_df = max(0, k - 2)
      overall_df = max(0, p * predictor_df)
      overall_chisq = ieee_value(0.0_dp, ieee_quiet_nan)
      overall_p = ieee_value(0.0_dp, ieee_quiet_nan)
      allocate(predictor_chisq(max(0, p)), predictor_p(max(0, p)))
      if (p > 0) then
         predictor_chisq = ieee_value(0.0_dp, ieee_quiet_nan)
         predictor_p = ieee_value(0.0_dp, ieee_quiet_nan)
      end if
      status = 0

      limit = 100
      if (present(max_iter)) limit = max_iter
      tolerance = 1.0e-10_dp
      if (present(tol)) tolerance = tol

      if (n < 1 .or. p < 1 .or. size(y_level) /= n .or. k < 3 .or. limit < 1 .or. tolerance <= 0.0_dp) then
         status = 1
         return
      end if
      if (any(y_level < 1) .or. any(y_level > k)) then
         status = 1
         return
      end if
      do i = 1, k
         if (count(y_level == i) == 0) then
            status = 1
            return
         end if
      end do
      if (any(ieee_is_nan(x)) .or. any(.not. ieee_is_finite(x))) then
         status = 2
         return
      end if

      allocate(w_case(n))
      w_case = 1.0_dp
      if (present(weights)) then
         if (size(weights) /= n .or. any(weights < 0.0_dp) .or. any(ieee_is_nan(weights)) .or. &
             any(.not. ieee_is_finite(weights)) .or. sum(weights) <= 0.0_dp) then
            status = 2
            return
         end if
         w_case = weights
      end if

      allocate(beta_all(p, n_threshold), fitted(n, n_threshold), binary(n))
      do j = 1, n_threshold
         where (y_level > j)
            binary = 1.0_dp
         elsewhere
            binary = 0.0_dp
         end where
         call logistic_fit(x, binary, w_case, beta_all(:, j), fitted(:, j), fit_status, limit, tolerance)
         if (fit_status /= 0) then
            status = 3
            return
         end if
      end do
      if (present(beta_threshold)) then
         allocate(beta_threshold(p, n_threshold))
         beta_threshold = beta_all
      end if

      allocate(covariance(p * n_threshold, p * n_threshold))
      allocate(w_cross(n), w_diag(n), v_jj(p + 1, p + 1), v_el(p + 1, p + 1), v_mid(p + 1, p + 1))
      covariance = 0.0_dp
      do el = 1, n_threshold
         do j = 1, el
            w_cross = (fitted(:, el) - fitted(:, j) * fitted(:, el)) * w_case
            w_diag = (fitted(:, j) - fitted(:, j) * fitted(:, j)) * w_case
            call weighted_information_inverse(x, w_diag, v_jj, inv_status)
            if (inv_status /= 0) then
               status = 4
               return
            end if
            w_diag = (fitted(:, el) - fitted(:, el) * fitted(:, el)) * w_case
            call weighted_information_inverse(x, w_diag, v_el, inv_status)
            if (inv_status /= 0) then
               status = 4
               return
            end if
            v_mid = matmul(transpose(x), spread(w_cross, 2, size(x, 2)) * x)
            v_mid = matmul(v_jj, matmul(v_mid, v_el))
            row0 = (j - 1) * p
            col0 = (el - 1) * p
            covariance(row0 + 1:row0 + p, col0 + 1:col0 + p) = v_mid(2:p + 1, 2:p + 1)
            covariance(col0 + 1:col0 + p, row0 + 1:row0 + p) = v_mid(2:p + 1, 2:p + 1)
         end do
      end do

      allocate(beta_flat(p * n_threshold))
      beta_flat = reshape(beta_all, [p * n_threshold])
      allocate(dmat(overall_df, p * n_threshold))
      dmat = 0.0_dp
      do j = 1, predictor_df
         row0 = (j - 1) * p
         do i = 1, p
            dmat(row0 + i, i) = 1.0_dp
            dmat(row0 + i, j * p + i) = -1.0_dp
         end do
      end do

      allocate(q(overall_df), dvcov(overall_df, overall_df), dvcov_inv(overall_df, overall_df))
      q = matmul(dmat, beta_flat)
      dvcov = matmul(dmat, matmul(covariance, transpose(dmat)))
      call matrix_inverse(dvcov, dvcov_inv, inv_status)
      if (inv_status /= 0) then
         status = 4
         return
      end if
      overall_chisq = dot_product(q, matmul(dvcov_inv, q))
      overall_p = chi_square_survival(overall_chisq, overall_df)

      allocate(dd(predictor_df, p * n_threshold), q_pred(predictor_df))
      allocate(pred_cov(predictor_df, predictor_df), pred_cov_inv(predictor_df, predictor_df))
      do i = 1, p
         dd = 0.0_dp
         do j = 1, predictor_df
            dd(j, i) = 1.0_dp
            dd(j, j * p + i) = -1.0_dp
         end do
         q_pred = matmul(dd, beta_flat)
         pred_cov = matmul(dd, matmul(covariance, transpose(dd)))
         call matrix_inverse(pred_cov, pred_cov_inv, inv_status)
         if (inv_status /= 0) then
            status = 4
            return
         end if
         predictor_chisq(i) = dot_product(q_pred, matmul(pred_cov_inv, q_pred))
         predictor_p(i) = chi_square_survival(predictor_chisq(i), predictor_df)
      end do

   contains

      subroutine logistic_fit(design, response, case_weight, slopes, probability, fit_status_local, iteration_limit, fit_tol)
         real(dp), intent(in) :: design(:, :) !! Logistic-regression design matrix, including an intercept column.
         real(dp), intent(in) :: response(:) !! Binary response values containing only zero or one.
         real(dp), intent(in) :: case_weight(:) !! Nonnegative case weights with one value per observation.
         real(dp), intent(out) :: slopes(:) !! Fitted non-intercept coefficients in design-column order.
         real(dp), intent(out) :: probability(:) !! Fitted Bernoulli probabilities for each observation.
         integer, intent(out) :: fit_status_local !! Zero on convergence; nonzero for singularity or iteration failure.
         integer, intent(in) :: iteration_limit !! Maximum number of IRLS coefficient updates.
         real(dp), intent(in) :: fit_tol !! Positive relative coefficient convergence threshold.
         real(dp), allocatable :: beta(:)
         real(dp), allocatable :: beta_new(:)
         real(dp), allocatable :: eta(:)
         real(dp), allocatable :: info(:, :)
         real(dp), allocatable :: info_inv(:, :)
         real(dp), allocatable :: rhs_local(:)
         real(dp), allocatable :: work_response(:)
         real(dp), allocatable :: work_weight(:)
         real(dp) :: diff
         real(dp) :: mean_response
         real(dp) :: scale
         integer :: iter
         integer :: m
         integer :: qcols

         m = size(design, 1)
         qcols = size(design, 2)
         allocate(beta(qcols), beta_new(qcols), eta(m), work_response(m), work_weight(m))
         allocate(info(qcols, qcols), info_inv(qcols, qcols), rhs_local(qcols))
         beta = 0.0_dp
         mean_response = sum(case_weight * response) / sum(case_weight)
         mean_response = min(1.0_dp - 1.0e-8_dp, max(1.0e-8_dp, mean_response))
         beta(1) = log(mean_response / (1.0_dp - mean_response))
         fit_status_local = 5

         do iter = 1, iteration_limit
            eta = matmul(design, beta)
            call logistic_probabilities(eta, probability)
            probability = min(1.0_dp - 1.0e-12_dp, max(1.0e-12_dp, probability))
            work_weight = case_weight * probability * (1.0_dp - probability)
            if (sum(work_weight) <= 0.0_dp) then
               fit_status_local = 1
               return
            end if
            work_response = eta + (response - probability) / (probability * (1.0_dp - probability))
            info = matmul(transpose(design), spread(work_weight, 2, qcols) * design)
            rhs_local = matmul(transpose(design), work_weight * work_response)
            call matrix_inverse(info, info_inv, inv_status)
            if (inv_status /= 0) then
               fit_status_local = 2
               return
            end if
            beta_new = matmul(info_inv, rhs_local)
            diff = maxval(abs(beta_new - beta))
            scale = 1.0_dp + maxval(abs(beta))
            beta = beta_new
            if (diff <= fit_tol * scale) then
               eta = matmul(design, beta)
               call logistic_probabilities(eta, probability)
               slopes = beta(2:qcols)
               fit_status_local = 0
               return
            end if
         end do
      end subroutine logistic_fit

      pure subroutine logistic_probabilities(eta, probability)
         real(dp), intent(in) :: eta(:) !! Linear predictors for a logistic model.
         real(dp), intent(out) :: probability(:) !! Numerically stable inverse-logit values corresponding to eta.
         integer :: ii

         do ii = 1, size(eta)
            if (eta(ii) >= 0.0_dp) then
               probability(ii) = 1.0_dp / (1.0_dp + exp(-eta(ii)))
            else
               probability(ii) = exp(eta(ii)) / (1.0_dp + exp(eta(ii)))
            end if
         end do
      end subroutine logistic_probabilities

      subroutine weighted_information_inverse(design, weight, inverse_info, inverse_status)
         real(dp), intent(in) :: design(:, :) !! Logistic-regression design matrix including an intercept.
         real(dp), intent(in) :: weight(:) !! Diagonal working weights associated with design rows.
         real(dp), intent(out) :: inverse_info(:, :) !! Inverse of transpose(design) times diag(weight) times design.
         integer, intent(out) :: inverse_status !! Zero on success; nonzero when the weighted information matrix is singular.
         real(dp), allocatable :: information(:, :)

         allocate(information(size(design, 2), size(design, 2)))
         information = matmul(transpose(design), spread(weight, 2, size(design, 2)) * design)
         call matrix_inverse(information, inverse_info, inverse_status)
      end subroutine weighted_information_inverse

   end subroutine po_test

end module car_po_test
