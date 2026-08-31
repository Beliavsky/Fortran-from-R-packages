! Substantive-model-compatible numerical kernels translated from jomo's SMC engines.
! Covers the reusable computations shared by jomo1smc*, jomo1ransmc*, and jomo2smc*.
! Upstream jomo 2.7-6 by Matteo Quartagno and James Carpenter; License: GPL-2.
! Modern Fortran translation, 2026. Distributed under GPL-2.0-only.
module jomo_substantive
   use jomo_kinds, only : dp
   use jomo_rng, only : rng_state, rng_uniform, rng_normal, rng_chisq
   use jomo_linalg, only : inverse_spd
   use jomo_distributions, only : mvnormal_sample
   implicit none
   private

   public :: normal_cdf
   public :: linear_loglik
   public :: binary_probit_loglik
   public :: ordinal_probit_loglik
   public :: cox_partial_loglik_ordered
   public :: sample_binary_probit_latent
   public :: sample_ordinal_probit_latent
   public :: update_ordinal_thresholds
   public :: sample_gaussian_coefficients
   public :: sample_linear_variance
   public :: cox_coordinate_newton

contains

   pure elemental function normal_cdf(x) result(value)
      real(dp), intent(in) :: x !! Standard-normal variate at which the cumulative probability is evaluated.
      real(dp) :: value

      value = 0.5_dp * erfc(-x / sqrt(2.0_dp))
   end function normal_cdf

   pure function linear_loglik(y, x, beta, variance) result(value)
      real(dp), intent(in) :: y(:) !! Continuous substantive-model response vector.
      real(dp), intent(in) :: x(:, :) !! Substantive fixed-effect design matrix, shape n by p.
      real(dp), intent(in) :: beta(:) !! Regression coefficient vector of length p.
      real(dp), intent(in) :: variance !! Positive Gaussian residual variance.
      real(dp) :: value
      real(dp), allocatable :: residual(:)

      if (size(x, 1) /= size(y) .or. size(x, 2) /= size(beta)) error stop "linear_loglik: shape mismatch"
      if (variance <= 0.0_dp) then
         value = -huge(1.0_dp)
         return
      end if
      residual = y - matmul(x, beta)
      value = -0.5_dp * sum(residual * residual) / variance - &
         0.5_dp * real(size(y), dp) * log(variance)
   end function linear_loglik

   pure function binary_probit_loglik(y, x, beta) result(value)
      integer, intent(in) :: y(:) !! Binary substantive outcome encoded exactly as 1 or 2.
      real(dp), intent(in) :: x(:, :) !! Substantive fixed-effect design matrix, shape n by p.
      real(dp), intent(in) :: beta(:) !! Probit regression coefficient vector of length p.
      real(dp) :: value
      integer :: i
      real(dp) :: eta
      real(dp) :: probability

      if (size(x, 1) /= size(y) .or. size(x, 2) /= size(beta)) error stop "binary_probit_loglik: shape mismatch"
      value = 0.0_dp
      do i = 1, size(y)
         eta = dot_product(x(i, :), beta)
         select case (y(i))
         case (1)
            probability = normal_cdf(-eta)
         case (2)
            probability = normal_cdf(eta)
         case default
            value = -huge(1.0_dp)
            return
         end select
         value = value + log(max(probability, tiny(1.0_dp)))
      end do
   end function binary_probit_loglik

   pure function ordinal_probit_loglik(y, x, beta, thresholds) result(value)
      integer, intent(in) :: y(:) !! Ordered categories encoded 1..K.
      real(dp), intent(in) :: x(:, :) !! Substantive fixed-effect design matrix, shape n by p.
      real(dp), intent(in) :: beta(:) !! Probit location coefficient vector of length p.
      real(dp), intent(in) :: thresholds(:) !! Strictly increasing K-1 latent-normal cut points.
      real(dp) :: value
      integer :: i
      integer :: k
      real(dp) :: eta
      real(dp) :: probability

      if (size(x, 1) /= size(y) .or. size(x, 2) /= size(beta)) error stop "ordinal_probit_loglik: shape mismatch"
      k = size(thresholds) + 1
      if (k < 2) error stop "ordinal_probit_loglik: at least two categories are required"
      if (size(thresholds) > 1) then
         if (any(thresholds(2:) <= thresholds(:size(thresholds) - 1))) &
            error stop "ordinal_probit_loglik: thresholds must be strictly increasing"
      end if
      value = 0.0_dp
      do i = 1, size(y)
         if (y(i) < 1 .or. y(i) > k) then
            value = -huge(1.0_dp)
            return
         end if
         eta = dot_product(x(i, :), beta)
         if (y(i) == 1) then
            probability = normal_cdf(thresholds(1) - eta)
         else if (y(i) == k) then
            probability = 1.0_dp - normal_cdf(thresholds(k - 1) - eta)
         else
            probability = normal_cdf(thresholds(y(i)) - eta) - normal_cdf(thresholds(y(i) - 1) - eta)
         end if
         value = value + log(max(probability, tiny(1.0_dp)))
      end do
   end function ordinal_probit_loglik

   pure function cox_partial_loglik_ordered(event, x, beta) result(value)
      logical, intent(in) :: event(:) !! Event indicators for rows sorted by increasing event/censoring time.
      real(dp), intent(in) :: x(:, :) !! Cox design matrix in the same increasing-time row order, shape n by p.
      real(dp), intent(in) :: beta(:) !! Cox log-hazard coefficient vector of length p.
      real(dp) :: value
      integer :: i
      real(dp) :: risk_sum
      real(dp), allocatable :: eta(:)

      if (size(x, 1) /= size(event) .or. size(x, 2) /= size(beta)) &
         error stop "cox_partial_loglik_ordered: shape mismatch"
      eta = matmul(x, beta)
      value = 0.0_dp
      risk_sum = 0.0_dp
      do i = size(event), 1, -1
         risk_sum = risk_sum + exp(eta(i))
         if (event(i)) value = value + eta(i) - log(risk_sum)
      end do
   end function cox_partial_loglik_ordered

   subroutine sample_binary_probit_latent(rng, y, x, beta, latent, max_tries)
      type(rng_state), intent(inout) :: rng !! Mutable generator state used for latent-normal rejection draws.
      integer, intent(in) :: y(:) !! Binary substantive outcome encoded exactly as 1 or 2.
      real(dp), intent(in) :: x(:, :) !! Probit fixed-effect design matrix, shape n by p.
      real(dp), intent(in) :: beta(:) !! Probit regression coefficients of length p.
      real(dp), intent(out) :: latent(:) !! Sampled unit-variance latent-normal outcomes, length n.
      integer, intent(in), optional :: max_tries !! Maximum rejection attempts per row; defaults to upstream value 10000.
      integer :: i
      integer :: tries
      integer :: limit
      real(dp) :: eta
      real(dp) :: draw
      logical :: accepted

      if (size(x, 1) /= size(y) .or. size(x, 2) /= size(beta) .or. size(latent) /= size(y)) &
         error stop "sample_binary_probit_latent: shape mismatch"
      limit = 10000
      if (present(max_tries)) limit = max_tries
      if (limit <= 0) error stop "sample_binary_probit_latent: max_tries must be positive"
      do i = 1, size(y)
         if (y(i) /= 1 .and. y(i) /= 2) error stop "sample_binary_probit_latent: y must contain only 1 or 2"
         eta = dot_product(x(i, :), beta)
         accepted = .false.
         do tries = 1, limit
            draw = rng_normal(rng, eta, 1.0_dp)
            if ((y(i) == 1 .and. draw < 0.0_dp) .or. (y(i) == 2 .and. draw > 0.0_dp)) then
               accepted = .true.
               exit
            end if
         end do
         if (.not. accepted) error stop "sample_binary_probit_latent: rejection sampler failed"
         latent(i) = draw
      end do
   end subroutine sample_binary_probit_latent

   subroutine sample_ordinal_probit_latent(rng, y, x, beta, thresholds, latent, max_tries)
      type(rng_state), intent(inout) :: rng !! Mutable generator state used for latent-normal rejection draws.
      integer, intent(in) :: y(:) !! Ordered substantive outcome encoded 1..K.
      real(dp), intent(in) :: x(:, :) !! Probit location design matrix, shape n by p.
      real(dp), intent(in) :: beta(:) !! Probit location coefficients of length p.
      real(dp), intent(in) :: thresholds(:) !! Strictly increasing K-1 latent-normal cut points.
      real(dp), intent(out) :: latent(:) !! Sampled unit-variance latent-normal outcomes, length n.
      integer, intent(in), optional :: max_tries !! Maximum rejection attempts per row; defaults to upstream value 10000.
      integer :: i
      integer :: k
      integer :: tries
      integer :: limit
      real(dp) :: eta
      real(dp) :: draw
      real(dp) :: lower
      real(dp) :: upper
      logical :: accepted

      if (size(x, 1) /= size(y) .or. size(x, 2) /= size(beta) .or. size(latent) /= size(y)) &
         error stop "sample_ordinal_probit_latent: shape mismatch"
      k = size(thresholds) + 1
      if (k < 2) error stop "sample_ordinal_probit_latent: at least two categories are required"
      if (size(thresholds) > 1) then
         if (any(thresholds(2:) <= thresholds(:size(thresholds) - 1))) &
            error stop "sample_ordinal_probit_latent: thresholds must be strictly increasing"
      end if
      limit = 10000
      if (present(max_tries)) limit = max_tries
      if (limit <= 0) error stop "sample_ordinal_probit_latent: max_tries must be positive"
      do i = 1, size(y)
         if (y(i) < 1 .or. y(i) > k) error stop "sample_ordinal_probit_latent: category out of range"
         eta = dot_product(x(i, :), beta)
         lower = -huge(1.0_dp)
         upper = huge(1.0_dp)
         if (y(i) > 1) lower = thresholds(y(i) - 1)
         if (y(i) < k) upper = thresholds(y(i))
         accepted = .false.
         do tries = 1, limit
            draw = rng_normal(rng, eta, 1.0_dp)
            if (draw > lower .and. draw < upper) then
               accepted = .true.
               exit
            end if
         end do
         if (.not. accepted) error stop "sample_ordinal_probit_latent: rejection sampler failed"
         latent(i) = draw
      end do
   end subroutine sample_ordinal_probit_latent

   subroutine update_ordinal_thresholds(rng, y, latent, thresholds, lower_bound, upper_bound)
      type(rng_state), intent(inout) :: rng !! Mutable generator state used for uniform threshold updates.
      integer, intent(in) :: y(:) !! Ordered category labels 1..K associated with latent values.
      real(dp), intent(in) :: latent(:) !! Current latent-normal outcome values.
      real(dp), intent(inout) :: thresholds(:) !! K-1 ordered cut points updated in place.
      real(dp), intent(in), optional :: lower_bound !! Finite lower search bound; defaults to upstream value -10.
      real(dp), intent(in), optional :: upper_bound !! Finite upper search bound; defaults to upstream value 10.
      integer :: i
      integer :: j
      integer :: k
      real(dp) :: lo
      real(dp) :: hi
      real(dp) :: global_lo
      real(dp) :: global_hi

      if (size(y) /= size(latent)) error stop "update_ordinal_thresholds: length mismatch"
      k = size(thresholds) + 1
      global_lo = -10.0_dp
      global_hi = 10.0_dp
      if (present(lower_bound)) global_lo = lower_bound
      if (present(upper_bound)) global_hi = upper_bound
      if (global_lo >= global_hi) error stop "update_ordinal_thresholds: invalid bounds"
      do j = 1, k - 1
         lo = global_lo
         hi = global_hi
         do i = 1, size(y)
            if (y(i) <= j) lo = max(lo, latent(i))
            if (y(i) > j) hi = min(hi, latent(i))
         end do
         if (lo >= hi) error stop "update_ordinal_thresholds: latent state is incompatible with category ordering"
         thresholds(j) = lo + (hi - lo) * rng_uniform(rng)
      end do
      if (size(thresholds) > 1) then
         if (any(thresholds(2:) <= thresholds(:size(thresholds) - 1))) &
            error stop "update_ordinal_thresholds: threshold update lost ordering"
      end if
   end subroutine update_ordinal_thresholds

   subroutine sample_gaussian_coefficients(rng, latent, x, variance, beta, info)
      type(rng_state), intent(inout) :: rng !! Mutable generator state used for the Gaussian regression coefficient draw.
      real(dp), intent(in) :: latent(:) !! Gaussian or latent-Gaussian substantive response vector.
      real(dp), intent(in) :: x(:, :) !! Full-rank fixed-effect design matrix, shape n by p.
      real(dp), intent(in) :: variance !! Positive residual variance; use 1 for probit models.
      real(dp), intent(out) :: beta(:) !! Draw from the flat-prior Gaussian coefficient posterior, length p.
      integer, intent(out) :: info !! Zero on success; nonzero if the posterior precision is singular or the draw fails.
      real(dp), allocatable :: precision(:, :)
      real(dp), allocatable :: covariance(:, :)
      real(dp), allocatable :: rhs(:)
      real(dp), allocatable :: mean(:)

      info = 0
      if (size(x, 1) /= size(latent) .or. size(x, 2) /= size(beta) .or. variance <= 0.0_dp) then
         info = -1
         return
      end if
      precision = matmul(transpose(x), x) / variance
      allocate(covariance(size(beta), size(beta)))
      rhs = matmul(transpose(x), latent) / variance
      call inverse_spd(precision, covariance, info)
      if (info /= 0) return
      mean = matmul(covariance, rhs)
      call mvnormal_sample(rng, mean, covariance, beta, info)
   end subroutine sample_gaussian_coefficients

   subroutine sample_linear_variance(rng, residual, prior_scale, variance)
      type(rng_state), intent(inout) :: rng !! Mutable generator state used for the one-dimensional Wishart-equivalent draw.
      real(dp), intent(in) :: residual(:) !! Current Gaussian substantive-model residual vector.
      real(dp), intent(in) :: prior_scale !! Nonnegative residual-scale contribution corresponding to upstream varY.prior.
      real(dp), intent(out) :: variance !! Sampled positive residual variance.
      real(dp) :: scale
      real(dp) :: chisq

      if (prior_scale < 0.0_dp) error stop "sample_linear_variance: prior_scale must be nonnegative"
      scale = sum(residual * residual) + prior_scale
      if (scale <= 0.0_dp) error stop "sample_linear_variance: posterior scale must be positive"
      chisq = rng_chisq(rng, real(size(residual) + 1, dp))
      variance = scale / chisq
   end subroutine sample_linear_variance

   subroutine cox_coordinate_newton(event, x, beta, max_iter, tolerance, converged)
      logical, intent(in) :: event(:) !! Event indicators for observations sorted by increasing survival time.
      real(dp), intent(in) :: x(:, :) !! Cox design matrix in increasing-time row order, shape n by p.
      real(dp), intent(inout) :: beta(:) !! Cox coefficient vector updated in place by upstream-style coordinate Newton steps.
      integer, intent(in), optional :: max_iter !! Maximum Newton iterations per coefficient; defaults to upstream value 15.
      real(dp), intent(in), optional :: tolerance !! Absolute coefficient convergence tolerance; defaults to 1e-5.
      logical, intent(out), optional :: converged !! True if every coefficient met the requested tolerance before its limit.
      integer :: coordinate
      integer :: i
      integer :: iter
      integer :: limit
      real(dp) :: tol
      real(dp) :: eta
      real(dp) :: exp_eta
      real(dp) :: risk0
      real(dp) :: risk1
      real(dp) :: risk2
      real(dp) :: score
      real(dp) :: hessian
      real(dp) :: old_beta
      real(dp) :: new_beta
      logical :: all_converged
      logical :: coordinate_converged

      if (size(x, 1) /= size(event) .or. size(x, 2) /= size(beta)) error stop "cox_coordinate_newton: shape mismatch"
      limit = 15
      if (present(max_iter)) limit = max_iter
      tol = 1.0e-5_dp
      if (present(tolerance)) tol = tolerance
      if (limit <= 0 .or. tol <= 0.0_dp) error stop "cox_coordinate_newton: invalid control parameter"
      all_converged = .true.
      do coordinate = 1, size(beta)
         coordinate_converged = .false.
         do iter = 1, limit
            old_beta = beta(coordinate)
            score = 0.0_dp
            hessian = 0.0_dp
            risk0 = 0.0_dp
            risk1 = 0.0_dp
            risk2 = 0.0_dp
            do i = size(event), 1, -1
               eta = dot_product(x(i, :), beta)
               exp_eta = exp(eta)
               risk0 = risk0 + exp_eta
               risk1 = risk1 + x(i, coordinate) * exp_eta
               risk2 = risk2 + x(i, coordinate) * x(i, coordinate) * exp_eta
               if (event(i)) then
                  score = score + x(i, coordinate) - risk1 / risk0
                  hessian = hessian - (risk2 / risk0 - (risk1 / risk0) * (risk1 / risk0))
               end if
            end do
            if (abs(hessian) <= tiny(1.0_dp)) exit
            new_beta = old_beta - score / hessian
            beta(coordinate) = new_beta
            if (abs(new_beta - old_beta) < tol) then
               coordinate_converged = .true.
               exit
            end if
         end do
         if (.not. coordinate_converged) all_converged = .false.
      end do
      if (present(converged)) converged = all_converged
   end subroutine cox_coordinate_newton

end module jomo_substantive
