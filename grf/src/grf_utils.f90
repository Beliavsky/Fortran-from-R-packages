module grf_utils
   use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
   use grf_kinds, only : dp
   implicit none
   private

   public :: weighted_mean
   public :: weighted_variance
   public :: weighted_covariance
   public :: weighted_quantile
   public :: solve_linear_system
   public :: invert_matrix
   public :: weighted_linear_regression
   public :: sort_indices_by_values
   public :: unique_sorted
   public :: standard_error_mean
   public :: objective_bayes_debias
   public :: cluster_robust_linear_regression
   public :: cluster_standard_error_mean

contains

   pure real(dp) function weighted_mean(x, w) result(value)
      real(dp), intent(in) :: x(:) !! Values whose weighted arithmetic mean is requested.
      real(dp), intent(in) :: w(:) !! Nonnegative observation weights with the same length as x.
      real(dp) :: sw

      sw = sum(w)
      if (size(x) == 0 .or. sw <= 0.0_dp) then
         value = 0.0_dp
      else
         value = sum(w * x) / sw
      end if
   end function weighted_mean

   pure real(dp) function weighted_variance(x, w) result(value)
      real(dp), intent(in) :: x(:) !! Values whose weighted population variance is requested.
      real(dp), intent(in) :: w(:) !! Nonnegative observation weights with the same length as x.
      real(dp) :: mu
      real(dp) :: sw

      sw = sum(w)
      if (size(x) == 0 .or. sw <= 0.0_dp) then
         value = 0.0_dp
         return
      end if
      mu = sum(w * x) / sw
      value = sum(w * (x - mu)**2) / sw
   end function weighted_variance

   pure real(dp) function weighted_covariance(x, y, w) result(value)
      real(dp), intent(in) :: x(:) !! First variable in the weighted population covariance.
      real(dp), intent(in) :: y(:) !! Second variable in the weighted population covariance.
      real(dp), intent(in) :: w(:) !! Nonnegative observation weights shared by x and y.
      real(dp) :: mx
      real(dp) :: my
      real(dp) :: sw

      sw = sum(w)
      if (size(x) == 0 .or. sw <= 0.0_dp) then
         value = 0.0_dp
         return
      end if
      mx = sum(w * x) / sw
      my = sum(w * y) / sw
      value = sum(w * (x - mx) * (y - my)) / sw
   end function weighted_covariance

   pure subroutine sort_indices_by_values(values, indices)
      real(dp), intent(in) :: values(:) !! Numeric values used as sort keys; NaNs are placed before finite values.
      integer, intent(out) :: indices(:) !! Permutation of 1:size(values) in nondecreasing key order.
      integer :: i
      integer :: j
      integer :: key
      real(dp) :: key_value
      real(dp) :: current

      if (size(indices) /= size(values)) return
      do i = 1, size(values)
         indices(i) = i
      end do
      do i = 2, size(values)
         key = indices(i)
         key_value = values(key)
         j = i - 1
         do while (j >= 1)
            current = values(indices(j))
            if (ieee_is_nan(key_value)) then
               if (ieee_is_nan(current)) exit
            else
               if (ieee_is_nan(current)) exit
               if (current <= key_value) exit
            end if
            indices(j + 1) = indices(j)
            j = j - 1
         end do
         indices(j + 1) = key
      end do
   end subroutine sort_indices_by_values

   pure subroutine weighted_quantile(values, weights, probability, value)
      real(dp), intent(in) :: values(:) !! Sample values entering the weighted empirical distribution.
      real(dp), intent(in) :: weights(:) !! Nonnegative sample weights with the same length as values.
      real(dp), intent(in) :: probability !! Requested probability in the closed interval [0,1].
      real(dp), intent(out) :: value !! Left-continuous weighted empirical quantile at probability.
      integer, allocatable :: order(:)
      real(dp) :: total_weight
      real(dp) :: target
      real(dp) :: cumulative
      integer :: i

      if (size(values) == 0) then
         value = 0.0_dp
         return
      end if
      allocate(order(size(values)))
      call sort_indices_by_values(values, order)
      total_weight = sum(max(weights, 0.0_dp))
      if (total_weight <= 0.0_dp) then
         value = values(order(max(1, min(size(values), nint(1.0_dp + probability * real(size(values) - 1, dp))))))
         return
      end if
      target = min(max(probability, 0.0_dp), 1.0_dp) * total_weight
      cumulative = 0.0_dp
      value = values(order(size(values)))
      do i = 1, size(values)
         cumulative = cumulative + max(weights(order(i)), 0.0_dp)
         if (cumulative >= target) then
            value = values(order(i))
            exit
         end if
      end do
   end subroutine weighted_quantile

   pure subroutine unique_sorted(values, unique_values)
      real(dp), intent(in) :: values(:) !! Input numeric vector; NaNs are omitted from the returned finite set.
      real(dp), allocatable, intent(out) :: unique_values(:) !! Sorted distinct finite values from the input.
      integer, allocatable :: order(:)
      real(dp), allocatable :: work(:)
      integer :: i
      integer :: count

      if (size(values) == 0) then
         allocate(unique_values(0))
         return
      end if
      allocate(order(size(values)))
      call sort_indices_by_values(values, order)
      allocate(work(size(values)))
      count = 0
      do i = 1, size(values)
         if (ieee_is_nan(values(order(i)))) cycle
         if (count == 0) then
            count = 1
            work(count) = values(order(i))
         else if (abs(values(order(i)) - work(count)) > 0.0_dp) then
            count = count + 1
            work(count) = values(order(i))
         end if
      end do
      allocate(unique_values(count))
      if (count > 0) unique_values = work(1:count)
   end subroutine unique_sorted

   pure subroutine solve_linear_system(a, b, x, info)
      real(dp), intent(in) :: a(:,:) !! Square coefficient matrix of the linear system.
      real(dp), intent(in) :: b(:) !! Right-hand side vector whose length equals size(a,1).
      real(dp), intent(out) :: x(:) !! Solution vector; set to zero if the matrix is singular.
      integer, intent(out) :: info !! Zero on success; positive pivot index if numerical singularity is detected.
      real(dp), allocatable :: aa(:,:)
      real(dp), allocatable :: bb(:)
      real(dp) :: pivot_value
      real(dp) :: factor
      real(dp) :: scale
      real(dp) :: tmp
      integer :: i
      integer :: j
      integer :: k
      integer :: pivot
      integer :: n

      n = size(b)
      x = 0.0_dp
      info = 0
      if (size(a,1) /= n .or. size(a,2) /= n .or. size(x) /= n) then
         info = -1
         return
      end if
      if (n == 0) return
      allocate(aa(n,n), bb(n))
      aa = a
      bb = b
      scale = max(1.0_dp, maxval(abs(aa)))
      do k = 1, n
         pivot = k
         pivot_value = abs(aa(k,k))
         do i = k + 1, n
            if (abs(aa(i,k)) > pivot_value) then
               pivot = i
               pivot_value = abs(aa(i,k))
            end if
         end do
         if (pivot_value <= 100.0_dp * epsilon(1.0_dp) * scale) then
            info = k
            return
         end if
         if (pivot /= k) then
            do j = k, n
               tmp = aa(k,j)
               aa(k,j) = aa(pivot,j)
               aa(pivot,j) = tmp
            end do
            tmp = bb(k)
            bb(k) = bb(pivot)
            bb(pivot) = tmp
         end if
         do i = k + 1, n
            factor = aa(i,k) / aa(k,k)
            aa(i,k) = 0.0_dp
            aa(i,k+1:n) = aa(i,k+1:n) - factor * aa(k,k+1:n)
            bb(i) = bb(i) - factor * bb(k)
         end do
      end do
      do i = n, 1, -1
         if (i < n) then
            x(i) = (bb(i) - dot_product(aa(i,i+1:n), x(i+1:n))) / aa(i,i)
         else
            x(i) = bb(i) / aa(i,i)
         end if
      end do
   end subroutine solve_linear_system

   pure subroutine invert_matrix(a, inverse, info)
      real(dp), intent(in) :: a(:,:) !! Square matrix to invert with repeated pivoted linear solves.
      real(dp), intent(out) :: inverse(:,:) !! Matrix inverse, or zeros if a solve fails.
      integer, intent(out) :: info !! Zero on success; otherwise the first failing pivot index.
      real(dp), allocatable :: rhs(:)
      real(dp), allocatable :: column(:)
      integer :: j
      integer :: n
      integer :: local_info

      n = size(a,1)
      inverse = 0.0_dp
      info = 0
      if (size(a,2) /= n .or. size(inverse,1) /= n .or. size(inverse,2) /= n) then
         info = -1
         return
      end if
      allocate(rhs(n), column(n))
      do j = 1, n
         rhs = 0.0_dp
         rhs(j) = 1.0_dp
         call solve_linear_system(a, rhs, column, local_info)
         if (local_info /= 0) then
            info = local_info
            inverse = 0.0_dp
            return
         end if
         inverse(:,j) = column
      end do
   end subroutine invert_matrix

   pure subroutine weighted_linear_regression(x, y, w, beta, std_err, info, ridge)
      real(dp), intent(in) :: x(:,:) !! Design matrix with observations in rows and regressors in columns.
      real(dp), intent(in) :: y(:) !! Response vector with one value per row of x.
      real(dp), intent(in) :: w(:) !! Nonnegative observation weights with one value per row of x.
      real(dp), intent(out) :: beta(:) !! Estimated weighted least-squares coefficients.
      real(dp), intent(out) :: std_err(:) !! HC0 sandwich standard errors for the fitted coefficients.
      integer, intent(out) :: info !! Zero on success; nonzero when dimensions are invalid or the normal matrix is singular.
      real(dp), intent(in), optional :: ridge !! Optional nonnegative diagonal ridge penalty added to the normal matrix.
      real(dp), allocatable :: xtwx(:,:)
      real(dp), allocatable :: meat(:,:)
      real(dp), allocatable :: inverse(:,:)
      real(dp), allocatable :: rhs(:)
      real(dp), allocatable :: residual(:)
      real(dp) :: penalty
      integer :: i
      integer :: j
      integer :: k
      integer :: n
      integer :: p
      integer :: inv_info

      n = size(x,1)
      p = size(x,2)
      beta = 0.0_dp
      std_err = 0.0_dp
      info = 0
      if (size(y) /= n .or. size(w) /= n .or. size(beta) /= p .or. size(std_err) /= p) then
         info = -1
         return
      end if
      if (n == 0 .or. p == 0) then
         info = -2
         return
      end if
      penalty = 0.0_dp
      if (present(ridge)) penalty = max(ridge, 0.0_dp)
      allocate(xtwx(p,p), rhs(p), inverse(p,p), residual(n), meat(p,p))
      xtwx = 0.0_dp
      rhs = 0.0_dp
      do i = 1, n
         do j = 1, p
            rhs(j) = rhs(j) + w(i) * x(i,j) * y(i)
            do k = 1, p
               xtwx(j,k) = xtwx(j,k) + w(i) * x(i,j) * x(i,k)
            end do
         end do
      end do
      do j = 1, p
         xtwx(j,j) = xtwx(j,j) + penalty
      end do
      call solve_linear_system(xtwx, rhs, beta, info)
      if (info /= 0) return
      residual = y - matmul(x, beta)
      meat = 0.0_dp
      do i = 1, n
         do j = 1, p
            do k = 1, p
               meat(j,k) = meat(j,k) + (w(i) * residual(i))**2 * x(i,j) * x(i,k)
            end do
         end do
      end do
      call invert_matrix(xtwx, inverse, inv_info)
      if (inv_info /= 0) then
         info = inv_info
         return
      end if
      meat = matmul(inverse, matmul(meat, transpose(inverse)))
      do j = 1, p
         std_err(j) = sqrt(max(meat(j,j), 0.0_dp))
      end do
   end subroutine weighted_linear_regression

   pure real(dp) function standard_error_mean(x, w) result(value)
      real(dp), intent(in) :: x(:) !! Values whose weighted mean standard error is estimated.
      real(dp), intent(in) :: w(:) !! Nonnegative observation weights corresponding to x.
      real(dp) :: mean_value
      real(dp) :: sw
      real(dp) :: effective_n

      sw = sum(w)
      if (size(x) <= 1 .or. sw <= 0.0_dp) then
         value = 0.0_dp
         return
      end if
      mean_value = sum(w * x) / sw
      effective_n = sw * sw / max(sum(w * w), tiny(1.0_dp))
      if (effective_n <= 1.0_dp) then
         value = 0.0_dp
      else
         value = sqrt(sum(w * (x - mean_value)**2) / sw / effective_n)
      end if
   end function standard_error_mean

   pure real(dp) function cluster_standard_error_mean(x, w, clusters) result(value)
      real(dp), intent(in) :: x(:) !! Values whose weighted mean is the estimand.
      real(dp), intent(in) :: w(:) !! Nonnegative observation weights matching x.
      integer, intent(in) :: clusters(:) !! Integer cluster labels; equal labels contribute one aggregated influence score.
      real(dp), allocatable :: score(:)
      integer, allocatable :: labels(:)
      real(dp) :: mu
      real(dp) :: sw
      integer :: g
      integer :: i
      integer :: ng
      integer :: found

      value = 0.0_dp
      if (size(x) <= 1 .or. size(w) /= size(x) .or. size(clusters) /= size(x)) return
      sw = sum(w)
      if (sw <= 0.0_dp) return
      allocate(labels(size(x)), score(size(x)))
      ng = 0
      score = 0.0_dp
      mu = weighted_mean(x, w)
      do i = 1, size(x)
         found = 0
         do g = 1, ng
            if (labels(g) == clusters(i)) then
               found = g
               exit
            end if
         end do
         if (found == 0) then
            ng = ng + 1
            labels(ng) = clusters(i)
            found = ng
         end if
         score(found) = score(found) + w(i) * (x(i) - mu)
      end do
      if (ng <= 1) return
      value = sqrt(max(0.0_dp, real(ng,dp) / real(ng - 1,dp) * sum(score(1:ng)**2) / (sw * sw)))
   end function cluster_standard_error_mean

   pure subroutine cluster_robust_linear_regression(x, y, w, clusters, beta, std_err, info, hc3)
      real(dp), intent(in) :: x(:,:) !! Design matrix with observations in rows and regressors in columns.
      real(dp), intent(in) :: y(:) !! Response vector matching rows of x.
      real(dp), intent(in) :: w(:) !! Nonnegative observation weights matching rows of x.
      integer, intent(in) :: clusters(:) !! Integer cluster labels used to aggregate sandwich scores.
      real(dp), intent(out) :: beta(:) !! Weighted least-squares coefficient estimates.
      real(dp), intent(out) :: std_err(:) !! Cluster-robust sandwich standard errors.
      integer, intent(out) :: info !! Zero on success; negative values indicate dimensions or a singular normal matrix.
      logical, intent(in), optional :: hc3 !! If true, apply HC3-style leverage adjustment before cluster aggregation.
      real(dp), allocatable :: bread_matrix(:,:)
      real(dp), allocatable :: bread(:,:)
      real(dp), allocatable :: rhs(:)
      real(dp), allocatable :: residual(:)
      real(dp), allocatable :: scores(:,:)
      real(dp), allocatable :: meat(:,:)
      integer, allocatable :: labels(:)
      real(dp) :: leverage
      real(dp) :: adjustment
      real(dp) :: correction
      integer :: n
      integer :: p
      integer :: ng
      integer :: found
      integer :: i
      integer :: j
      integer :: k
      integer :: g
      logical :: use_hc3

      n = size(x,1)
      p = size(x,2)
      beta = 0.0_dp
      std_err = 0.0_dp
      info = 0
      if (size(y) /= n .or. size(w) /= n .or. size(clusters) /= n .or. size(beta) /= p .or. size(std_err) /= p) then
         info = -1
         return
      end if
      if (n <= p .or. p < 1) then
         info = -2
         return
      end if
      allocate(bread_matrix(p,p), bread(p,p), rhs(p), residual(n))
      bread_matrix = 0.0_dp
      rhs = 0.0_dp
      do i = 1, n
         do j = 1, p
            rhs(j) = rhs(j) + w(i) * x(i,j) * y(i)
            do k = 1, p
               bread_matrix(j,k) = bread_matrix(j,k) + w(i) * x(i,j) * x(i,k)
            end do
         end do
      end do
      call solve_linear_system(bread_matrix, rhs, beta, info)
      if (info /= 0) return
      call invert_matrix(bread_matrix, bread, info)
      if (info /= 0) return
      residual = y - matmul(x, beta)
      allocate(labels(n), scores(n,p))
      scores = 0.0_dp
      ng = 0
      use_hc3 = .false.
      if (present(hc3)) use_hc3 = hc3
      do i = 1, n
         found = 0
         do g = 1, ng
            if (labels(g) == clusters(i)) then
               found = g
               exit
            end if
         end do
         if (found == 0) then
            ng = ng + 1
            labels(ng) = clusters(i)
            found = ng
         end if
         adjustment = 1.0_dp
         if (use_hc3) then
            leverage = w(i) * dot_product(x(i,:), matmul(bread, x(i,:)))
            adjustment = 1.0_dp / max(1.0_dp - leverage, sqrt(epsilon(1.0_dp)))
         end if
         scores(found,:) = scores(found,:) + w(i) * x(i,:) * residual(i) * adjustment
      end do
      if (ng <= 1) then
         info = -3
         return
      end if
      allocate(meat(p,p))
      meat = 0.0_dp
      do g = 1, ng
         do j = 1, p
            do k = 1, p
               meat(j,k) = meat(j,k) + scores(g,j) * scores(g,k)
            end do
         end do
      end do
      correction = real(ng,dp) / real(ng - 1,dp)
      if (n > p) correction = correction * real(n - 1,dp) / real(n - p,dp)
      meat = correction * matmul(bread, matmul(meat, transpose(bread)))
      do j = 1, p
         std_err(j) = sqrt(max(0.0_dp, meat(j,j)))
      end do
   end subroutine cluster_robust_linear_regression

   pure real(dp) function objective_bayes_debias(var_between, group_noise, num_good_groups) result(value)
      real(dp), intent(in) :: var_between !! Raw between-CI-group second moment of the influence pseudo-outcome.
      real(dp), intent(in) :: group_noise !! Estimated finite-group inflation of var_between from within-group Monte Carlo noise.
      real(dp), intent(in) :: num_good_groups !! Number of complete nonempty CI groups contributing to the variance calculation.
      real(dp), parameter :: one_over_sqrt_two_pi = 0.39894228040143267794_dp
      real(dp), parameter :: one_over_sqrt_two = 0.70710678118654752440_dp
      real(dp) :: initial_estimate
      real(dp) :: initial_se
      real(dp) :: ratio
      real(dp) :: numerator
      real(dp) :: denominator

      if (num_good_groups <= 0.0_dp) then
         value = 0.0_dp
         return
      end if
      initial_estimate = var_between - group_noise
      initial_se = max(abs(var_between), abs(group_noise)) * sqrt(2.0_dp / num_good_groups)
      if (initial_se <= 1.0e-10_dp) then
         value = 0.0_dp
         return
      end if
      ratio = initial_estimate / initial_se
      numerator = exp(-0.5_dp * ratio * ratio) * one_over_sqrt_two_pi
      denominator = 0.5_dp * erfc(-ratio * one_over_sqrt_two)
      if (denominator <= tiny(1.0_dp)) then
         value = max(initial_estimate, 0.0_dp)
      else
         value = initial_estimate + initial_se * numerator / denominator
         value = max(value, 0.0_dp)
      end if
   end function objective_bayes_debias

end module grf_utils
