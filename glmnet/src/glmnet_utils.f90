! SPDX-License-Identifier: GPL-2.0-only
module glmnet_utils
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_value, ieee_quiet_nan
   use glmnet_kinds, only : dp, glmnet_eps, glmnet_huge
   use glmnet_status, only : glmnet_success, glmnet_invalid_argument, &
      glmnet_nonfinite_input, glmnet_all_predictors_constant
   use glmnet_types, only : glmnet_control_type, glmnet_sparse_csc
   implicit none
   private
   public :: soft_threshold, group_soft_threshold, logistic, log1pexp
   public :: normalize_weights, prepare_design, make_lambda_sequence
   public :: count_nonzero_rows, clamp_value, all_finite_vector, all_finite_matrix
   public :: sparse_to_dense, dense_to_sparse, deterministic_fold_ids
   public :: weighted_mean, weighted_variance, linear_interpolate_coefficients
   public :: sort_indices_real, unique_sorted_real, safe_log, safe_exp
   public :: projected_simplex
contains
   pure elemental function soft_threshold(x, threshold) result(value)
      real(dp), intent(in) :: x, threshold
      real(dp) :: value
      value = sign(max(abs(x) - max(threshold, 0.0_dp), 0.0_dp), x)
   end function soft_threshold

   pure subroutine group_soft_threshold(x, threshold, value)
      real(dp), intent(in) :: x(:), threshold
      real(dp), intent(out) :: value(size(x))
      real(dp) :: norm_x, scale
      norm_x = sqrt(max(dot_product(x, x), 0.0_dp))
      if (norm_x <= max(threshold, 0.0_dp)) then
         value = 0.0_dp
      else
         scale = 1.0_dp - max(threshold, 0.0_dp) / norm_x
         value = scale * x
      end if
   end subroutine group_soft_threshold

   pure elemental function logistic(x) result(value)
      real(dp), intent(in) :: x
      real(dp) :: value
      if (x >= 0.0_dp) then
         value = 1.0_dp / (1.0_dp + exp(-min(x, 700.0_dp)))
      else
         value = exp(max(x, -700.0_dp)) / (1.0_dp + exp(max(x, -700.0_dp)))
      end if
   end function logistic

   pure elemental function log1pexp(x) result(value)
      real(dp), intent(in) :: x
      real(dp) :: value
      if (x > 35.0_dp) then
         value = x
      else if (x < -35.0_dp) then
         value = exp(x)
      else
         value = log(1.0_dp + exp(x))
      end if
   end function log1pexp

   pure elemental function safe_log(x) result(value)
      real(dp), intent(in) :: x
      real(dp) :: value
      value = log(max(x, tiny(1.0_dp)))
   end function safe_log

   pure elemental function safe_exp(x) result(value)
      real(dp), intent(in) :: x
      real(dp) :: value
      value = exp(max(min(x, log(glmnet_huge) - 2.0_dp), log(tiny(1.0_dp)) + 2.0_dp))
   end function safe_exp

   pure elemental function clamp_value(x, lower, upper) result(value)
      real(dp), intent(in) :: x, lower, upper
      real(dp) :: value
      value = min(max(x, lower), upper)
   end function clamp_value

   pure function all_finite_vector(x) result(ok)
      real(dp), intent(in) :: x(:)
      logical :: ok
      integer :: i
      ok = .true.
      do i = 1, size(x)
         if (.not. ieee_is_finite(x(i))) then
            ok = .false.
            return
         end if
      end do
   end function all_finite_vector

   pure function all_finite_matrix(x) result(ok)
      real(dp), intent(in) :: x(:,:)
      logical :: ok
      integer :: i, j
      ok = .true.
      do j = 1, size(x, 2)
         do i = 1, size(x, 1)
            if (.not. ieee_is_finite(x(i, j))) then
               ok = .false.
               return
            end if
         end do
      end do
   end function all_finite_matrix

   subroutine normalize_weights(n, weights_in, weights, status)
      integer, intent(in) :: n
      real(dp), intent(in), optional :: weights_in(:)
      real(dp), allocatable, intent(out) :: weights(:)
      integer, intent(out) :: status
      real(dp) :: total
      status = glmnet_success
      allocate(weights(n))
      if (present(weights_in)) then
         if (size(weights_in) /= n) then
            status = glmnet_invalid_argument
            weights = 0.0_dp
            return
         end if
         weights = weights_in
      else
         weights = 1.0_dp
      end if
      if (.not. all_finite_vector(weights) .or. any(weights < 0.0_dp)) then
         status = glmnet_nonfinite_input
         weights = 0.0_dp
         return
      end if
      total = sum(weights)
      if (total <= glmnet_eps) then
         status = glmnet_invalid_argument
         weights = 0.0_dp
         return
      end if
      weights = weights / total
   end subroutine normalize_weights

   pure function weighted_mean(x, weights) result(value)
      real(dp), intent(in) :: x(:), weights(:)
      real(dp) :: value
      value = sum(weights * x) / max(sum(weights), glmnet_eps)
   end function weighted_mean

   pure function weighted_variance(x, weights, mean_value) result(value)
      real(dp), intent(in) :: x(:), weights(:)
      real(dp), intent(in), optional :: mean_value
      real(dp) :: value, mu
      if (present(mean_value)) then
         mu = mean_value
      else
         mu = weighted_mean(x, weights)
      end if
      value = sum(weights * (x - mu) ** 2) / max(sum(weights), glmnet_eps)
   end function weighted_variance

   subroutine prepare_design(x, weights, control, x_work, x_mean, x_scale, usable, status)
      real(dp), intent(in) :: x(:,:), weights(:)
      type(glmnet_control_type), intent(in) :: control
      real(dp), allocatable, intent(out) :: x_work(:,:), x_mean(:), x_scale(:)
      logical, allocatable, intent(out) :: usable(:)
      integer, intent(out) :: status
      integer :: n, p, j
      real(dp) :: variance
      n = size(x, 1)
      p = size(x, 2)
      status = glmnet_success
      allocate(x_work(n, p), x_mean(p), x_scale(p), usable(p))
      if (.not. all_finite_matrix(x)) then
         status = glmnet_nonfinite_input
         x_work = 0.0_dp
         x_mean = 0.0_dp
         x_scale = 1.0_dp
         usable = .false.
         return
      end if
      do j = 1, p
         if (control%intercept) then
            x_mean(j) = weighted_mean(x(:, j), weights)
         else
            x_mean(j) = 0.0_dp
         end if
         variance = sum(weights * (x(:, j) - x_mean(j)) ** 2)
         usable(j) = variance > 100.0_dp * glmnet_eps
         if (control%standardize .and. usable(j)) then
            x_scale(j) = sqrt(variance)
         else
            x_scale(j) = 1.0_dp
         end if
         if (usable(j)) then
            x_work(:, j) = (x(:, j) - x_mean(j)) / x_scale(j)
         else
            x_work(:, j) = 0.0_dp
         end if
      end do
      if (.not. any(usable)) status = glmnet_all_predictors_constant
   end subroutine prepare_design

   subroutine make_lambda_sequence(lambda_max, n, p, control, lambda_in, lambda, status)
      real(dp), intent(in) :: lambda_max
      integer, intent(in) :: n, p
      type(glmnet_control_type), intent(in) :: control
      real(dp), intent(in), optional :: lambda_in(:)
      real(dp), allocatable, intent(out) :: lambda(:)
      integer, intent(out) :: status
      integer :: l, nl
      real(dp) :: ratio, lmax, log_min
      status = glmnet_success
      if (present(lambda_in)) then
         if (size(lambda_in) < 1) then
            status = glmnet_invalid_argument
            allocate(lambda(0))
            return
         end if
         if (.not. all_finite_vector(lambda_in) .or. any(lambda_in < 0.0_dp)) then
            status = glmnet_invalid_argument
            allocate(lambda(0))
            return
         end if
         allocate(lambda(size(lambda_in)))
         lambda = lambda_in
         call sort_descending(lambda)
         return
      end if
      nl = max(control%nlambda, 1)
      allocate(lambda(nl))
      lmax = max(lambda_max, 100.0_dp * glmnet_eps)
      if (control%lambda_min_ratio > 0.0_dp) then
         ratio = control%lambda_min_ratio
      else if (n > p) then
         ratio = 1.0e-4_dp
      else
         ratio = 1.0e-2_dp
      end if
      ratio = min(max(ratio, 1.0e-8_dp), 1.0_dp)
      if (nl == 1) then
         lambda(1) = lmax
      else
         log_min = log(ratio)
         do l = 1, nl
            lambda(l) = lmax * exp(real(l - 1, dp) * log_min / real(nl - 1, dp))
         end do
      end if
   contains
      subroutine sort_descending(values)
         real(dp), intent(inout) :: values(:)
         integer :: i, j
         real(dp) :: key
         do i = 2, size(values)
            key = values(i)
            j = i - 1
            do while (j >= 1)
               if (values(j) >= key) exit
               values(j + 1) = values(j)
               j = j - 1
            end do
            values(j + 1) = key
         end do
      end subroutine sort_descending
   end subroutine make_lambda_sequence

   pure function count_nonzero_rows(beta, tolerance) result(count)
      real(dp), intent(in) :: beta(:,:)
      real(dp), intent(in), optional :: tolerance
      integer :: count, j
      real(dp) :: tol
      tol = 1.0e-12_dp
      if (present(tolerance)) tol = max(tolerance, 0.0_dp)
      count = 0
      do j = 1, size(beta, 1)
         if (maxval(abs(beta(j, :))) > tol) count = count + 1
      end do
   end function count_nonzero_rows

   subroutine sparse_to_dense(sparse, dense, status)
      type(glmnet_sparse_csc), intent(in) :: sparse
      real(dp), allocatable, intent(out) :: dense(:,:)
      integer, intent(out) :: status
      integer :: j, k, i, first, last
      status = glmnet_success
      if (sparse%nrow < 1 .or. sparse%ncol < 1) then
         status = glmnet_invalid_argument
         allocate(dense(0, 0))
         return
      end if
      if (.not. allocated(sparse%values) .or. .not. allocated(sparse%row_index) .or. &
          .not. allocated(sparse%col_pointer)) then
         status = glmnet_invalid_argument
         allocate(dense(0, 0))
         return
      end if
      if (size(sparse%values) /= size(sparse%row_index) .or. &
          size(sparse%col_pointer) /= sparse%ncol + 1) then
         status = glmnet_invalid_argument
         allocate(dense(0, 0))
         return
      end if
      allocate(dense(sparse%nrow, sparse%ncol))
      dense = 0.0_dp
      do j = 1, sparse%ncol
         first = sparse%col_pointer(j)
         last = sparse%col_pointer(j + 1) - 1
         do k = first, last
            if (k < 1 .or. k > size(sparse%values)) then
               status = glmnet_invalid_argument
               return
            end if
            i = sparse%row_index(k)
            if (i < 1 .or. i > sparse%nrow) then
               status = glmnet_invalid_argument
               return
            end if
            dense(i, j) = dense(i, j) + sparse%values(k)
         end do
      end do
   end subroutine sparse_to_dense

   subroutine dense_to_sparse(dense, sparse, tolerance)
      real(dp), intent(in) :: dense(:,:)
      type(glmnet_sparse_csc), intent(out) :: sparse
      real(dp), intent(in), optional :: tolerance
      real(dp) :: tol
      integer :: i, j, k, nnz
      tol = 0.0_dp
      if (present(tolerance)) tol = max(tolerance, 0.0_dp)
      sparse%nrow = size(dense, 1)
      sparse%ncol = size(dense, 2)
      nnz = count(abs(dense) > tol)
      allocate(sparse%values(nnz), sparse%row_index(nnz), sparse%col_pointer(sparse%ncol + 1))
      k = 1
      sparse%col_pointer(1) = 1
      do j = 1, sparse%ncol
         do i = 1, sparse%nrow
            if (abs(dense(i, j)) > tol) then
               sparse%values(k) = dense(i, j)
               sparse%row_index(k) = i
               k = k + 1
            end if
         end do
         sparse%col_pointer(j + 1) = k
      end do
   end subroutine dense_to_sparse

   subroutine deterministic_fold_ids(n, nfolds, fold_id, seed)
      integer, intent(in) :: n, nfolds
      integer, allocatable, intent(out) :: fold_id(:)
      integer, intent(in), optional :: seed
      integer(kind=8) :: state
      integer :: i, j, tmp
      integer, allocatable :: order(:)
      allocate(fold_id(n), order(n))
      do i = 1, n
         order(i) = i
      end do
      state = 88172645463325252_8
      if (present(seed)) state = int(abs(seed) + 1, kind=8)
      do i = n, 2, -1
         state = ieor(state, shiftl(state, 13))
         state = ieor(state, shiftr(state, 7))
         state = ieor(state, shiftl(state, 17))
         j = 1 + int(modulo(iand(state, int(z'7FFFFFFFFFFFFFFF', kind=8)), int(i, kind=8)))
         tmp = order(i)
         order(i) = order(j)
         order(j) = tmp
      end do
      do i = 1, n
         fold_id(order(i)) = 1 + modulo(i - 1, max(nfolds, 1))
      end do
   end subroutine deterministic_fold_ids

   subroutine linear_interpolate_coefficients(lambda, intercept, beta, s, a, b)
      real(dp), intent(in) :: lambda(:), intercept(:,:), beta(:,:,:), s
      real(dp), intent(out) :: a(size(intercept, 1))
      real(dp), intent(out) :: b(size(beta, 1), size(beta, 2))
      integer :: left, right, l
      real(dp) :: fraction, log_s, denom
      if (s >= lambda(1)) then
         a = intercept(:, 1)
         b = beta(:, :, 1)
         return
      end if
      if (s <= lambda(size(lambda))) then
         a = intercept(:, size(lambda))
         b = beta(:, :, size(lambda))
         return
      end if
      left = 1
      right = 2
      do l = 1, size(lambda) - 1
         if (lambda(l) >= s .and. s >= lambda(l + 1)) then
            left = l
            right = l + 1
            exit
         end if
      end do
      log_s = log(max(s, tiny(1.0_dp)))
      denom = log(max(lambda(right), tiny(1.0_dp))) - log(max(lambda(left), tiny(1.0_dp)))
      if (abs(denom) <= glmnet_eps) then
         fraction = 0.0_dp
      else
         fraction = (log_s - log(max(lambda(left), tiny(1.0_dp)))) / denom
      end if
      fraction = min(max(fraction, 0.0_dp), 1.0_dp)
      a = (1.0_dp - fraction) * intercept(:, left) + fraction * intercept(:, right)
      b = (1.0_dp - fraction) * beta(:, :, left) + fraction * beta(:, :, right)
   end subroutine linear_interpolate_coefficients

   subroutine sort_indices_real(values, indices, descending)
      real(dp), intent(in) :: values(:)
      integer, allocatable, intent(out) :: indices(:)
      logical, intent(in), optional :: descending
      logical :: desc
      integer :: i, j, key
      desc = .false.
      if (present(descending)) desc = descending
      allocate(indices(size(values)))
      do i = 1, size(values)
         indices(i) = i
      end do
      do i = 2, size(values)
         key = indices(i)
         j = i - 1
         do while (j >= 1)
            if (desc) then
               if (values(indices(j)) >= values(key)) exit
            else
               if (values(indices(j)) <= values(key)) exit
            end if
            indices(j + 1) = indices(j)
            j = j - 1
         end do
         indices(j + 1) = key
      end do
   end subroutine sort_indices_real

   subroutine unique_sorted_real(values, unique_values, tolerance)
      real(dp), intent(in) :: values(:)
      real(dp), allocatable, intent(out) :: unique_values(:)
      real(dp), intent(in), optional :: tolerance
      integer, allocatable :: order(:)
      real(dp), allocatable :: work(:)
      real(dp) :: tol
      integer :: i, count_unique
      tol = 0.0_dp
      if (present(tolerance)) tol = max(tolerance, 0.0_dp)
      if (size(values) == 0) then
         allocate(unique_values(0))
         return
      end if
      call sort_indices_real(values, order)
      allocate(work(size(values)))
      count_unique = 1
      work(1) = values(order(1))
      do i = 2, size(values)
         if (abs(values(order(i)) - work(count_unique)) > tol) then
            count_unique = count_unique + 1
            work(count_unique) = values(order(i))
         end if
      end do
      allocate(unique_values(count_unique))
      unique_values = work(:count_unique)
   end subroutine unique_sorted_real

   subroutine projected_simplex(x, target_sum, result)
      real(dp), intent(in) :: x(:), target_sum
      real(dp), intent(out) :: result(size(x))
      integer, allocatable :: order(:)
      real(dp) :: cumulative, theta
      integer :: j, rho
      if (target_sum <= 0.0_dp) then
         result = 0.0_dp
         return
      end if
      call sort_indices_real(x, order, descending=.true.)
      cumulative = 0.0_dp
      rho = 1
      theta = 0.0_dp
      do j = 1, size(x)
         cumulative = cumulative + x(order(j))
         theta = (cumulative - target_sum) / real(j, dp)
         if (x(order(j)) - theta > 0.0_dp) rho = j
      end do
      cumulative = sum(x(order(:rho)))
      theta = (cumulative - target_sum) / real(rho, dp)
      result = max(x - theta, 0.0_dp)
   end subroutine projected_simplex
end module glmnet_utils
