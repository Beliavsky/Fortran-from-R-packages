! SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only
!
! Multivariate phylogenetic correlation model translated from ape R/corphylo.R.
! Formula/data-frame plumbing and stochastic SANN optimization are omitted.
module ape_corphylo
   use r_kinds, only : dp
   use r_linalg, only : spd_inverse_logdet
   use ape_types, only : phylo_tree
   use ape_topology, only : phylogenetic_vcv
   use ape_optimize, only : bounded_problem, bounded_bfgs
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   implicit none
   private

   real(dp), parameter :: pi = 3.1415926535897932384626433832795029_dp
   real(dp), parameter :: two_pi = 2.0_dp * pi

   type, public :: corphylo_result
      real(dp), allocatable :: correlation(:, :)
      real(dp), allocatable :: d(:)
      real(dp), allocatable :: coefficients(:)
      real(dp), allocatable :: coefficient_se(:)
      real(dp), allocatable :: coefficient_covariance(:, :)
      real(dp), allocatable :: process_covariance(:, :)
      real(dp), allocatable :: observation_covariance(:, :)
      real(dp), allocatable :: phylogenetic_covariance(:, :)
      real(dp) :: log_likelihood = -huge(1.0_dp)
      real(dp) :: aic = huge(1.0_dp)
      real(dp) :: bic = huge(1.0_dp)
      integer :: iterations = 0
      logical :: converged = .false.
      logical :: reml = .true.
   end type corphylo_result

   type, extends(bounded_problem) :: corphylo_problem
      real(dp), allocatable :: xx(:)
      real(dp), allocatable :: uu(:, :)
      real(dp), allocatable :: mm(:)
      real(dp), allocatable :: tau(:, :)
      real(dp), allocatable :: vphy(:, :)
      integer :: n = 0
      integer :: p = 0
      logical :: reml = .true.
      logical :: constrain_d = .false.
   contains
      procedure :: value => corphylo_problem_value
      procedure :: gradient => corphylo_problem_gradient
   end type corphylo_problem

   public :: corphylo_fit
   public :: corphylo_objective

contains

   subroutine corphylo_fit(x, tree, result, info, measurement_error, covariates, reml, constrain_d, max_iter, tolerance)
      !! Fits ape's multivariate OU `corphylo` likelihood with deterministic bounded BFGS.
      real(dp), intent(in) :: x(:, :) !! Species-by-trait matrix in numeric tip order; each trait must have positive sample SD.
      type(phylo_tree), intent(in) :: tree !! Rooted tree with branch lengths and one tip for each row of `x`.
      type(corphylo_result), intent(out) :: result !! Trait correlations, OU d values, coefficients, covariance, and IC statistics.
      integer, intent(out) :: info !! Zero on success or nonzero for invalid input, singular covariance, or optimizer failure.
      real(dp), intent(in), optional :: measurement_error(:, :) !! Species-by-trait standard errors; default all zeros.
      real(dp), intent(in), optional :: covariates(:, :, :) !! Species-by-predictor-by-trait covariates; constant columns
         !! are ignored.
      logical, intent(in), optional :: reml !! Use REML objective if true; default true as in upstream `corphylo`.
      logical, intent(in), optional :: constrain_d !! Constrain each OU d to `(0,1)` by a logit parameterization; default false.
      integer, intent(in), optional :: max_iter !! Maximum bounded-BFGS iterations; default 2000.
      real(dp), intent(in), optional :: tolerance !! Projected-gradient tolerance; default `1e-7`.
      type(corphylo_problem) :: problem
      real(dp), allocatable :: coefficient_scale(:)
      real(dp), allocatable :: covariance_b(:, :)
      real(dp), allocatable :: lower(:)
      real(dp), allocatable :: lstart(:, :)
      real(dp), allocatable :: parameters(:)
      real(dp), allocatable :: upper(:)
      real(dp), allocatable :: x_mean(:)
      real(dp), allocatable :: x_sd(:)
      real(dp) :: objective
      real(dp) :: tol
      integer :: iterations
      integer :: n_l
      integer :: n_parameters
      integer :: opt_info
      integer :: status

      result = corphylo_result()
      info = 0
      call setup_corphylo(x, tree, measurement_error, covariates, problem, x_mean, x_sd, coefficient_scale, status)
      if (status /= 0) then
         info = status
         return
      end if
      problem%reml = .true.
      if (present(reml)) problem%reml = reml
      problem%constrain_d = .false.
      if (present(constrain_d)) problem%constrain_d = constrain_d
      tol = 1.0e-7_dp
      if (present(tolerance)) tol = tolerance
      if (tol <= 0.0_dp .or. .not. ieee_is_finite(tol)) then
         info = 4
         return
      end if

      n_l = problem%p * (problem%p + 1) / 2
      n_parameters = n_l + problem%p
      allocate(parameters(n_parameters), lower(n_parameters), upper(n_parameters))
      call initial_cholesky(problem%xx, problem%n, problem%p, lstart, status)
      if (status /= 0) then
         info = 10 + status
         return
      end if
      call pack_lower(lstart, parameters(1:n_l))
      lower(1:n_l) = -20.0_dp
      upper(1:n_l) = 20.0_dp
      if (problem%constrain_d) then
         parameters(n_l + 1:) = 0.0_dp
         lower(n_l + 1:) = -10.0_dp
         upper(n_l + 1:) = 10.0_dp
      else
         parameters(n_l + 1:) = 0.5_dp
         lower(n_l + 1:) = 1.0e-6_dp
         upper(n_l + 1:) = 10.0_dp
      end if
      call bounded_bfgs(problem, parameters, lower, upper, objective, opt_info, iterations, max_iter, tol)
      if (opt_info /= 0 .and. opt_info /= 7) then
         info = 20 + opt_info
         return
      end if
      call corphylo_finish(problem, parameters, x_mean, coefficient_scale, result, covariance_b, status)
      if (status /= 0) then
         info = 30 + status
         return
      end if
      result%iterations = iterations
      result%converged = opt_info == 0
      result%reml = problem%reml
   end subroutine corphylo_fit

   function corphylo_objective(x, tree, l_elements, d, info, measurement_error, covariates, reml) result(value)
      !! Evaluates upstream `corphylo.LL` for fixed covariance-factor and OU parameters after standardizing raw inputs.
      real(dp), intent(in) :: x(:, :) !! Species-by-trait matrix in numeric tip order.
      type(phylo_tree), intent(in) :: tree !! Rooted tree defining the normalized phylogenetic covariance.
      real(dp), intent(in) :: l_elements(:) !! Lower-triangle elements of L in R column order, with `R = transpose(L)*L`.
      real(dp), intent(in) :: d(:) !! Positive OU d values, one for each trait.
      integer, intent(out) :: info !! Zero on success or nonzero for invalid dimensions/covariance factorization.
      real(dp), intent(in), optional :: measurement_error(:, :) !! Species-by-trait measurement standard errors.
      real(dp), intent(in), optional :: covariates(:, :, :) !! Species-by-predictor-by-trait covariates.
      logical, intent(in), optional :: reml !! Evaluate REML if true; default true.
      real(dp) :: value
      type(corphylo_problem) :: problem
      real(dp), allocatable :: coefficient_scale(:)
      real(dp), allocatable :: parameters(:)
      real(dp), allocatable :: x_mean(:)
      real(dp), allocatable :: x_sd(:)
      integer :: n_l
      integer :: status

      value = huge(1.0_dp) / 100.0_dp
      info = 0
      call setup_corphylo(x, tree, measurement_error, covariates, problem, x_mean, x_sd, coefficient_scale, status)
      if (status /= 0) then
         info = status
         return
      end if
      problem%reml = .true.
      if (present(reml)) problem%reml = reml
      problem%constrain_d = .false.
      n_l = problem%p * (problem%p + 1) / 2
      if (size(l_elements) /= n_l .or. size(d) /= problem%p .or. any(d <= 0.0_dp) .or. any(d > 10.0_dp)) then
         info = 4
         return
      end if
      allocate(parameters(n_l + problem%p))
      parameters(1:n_l) = l_elements
      parameters(n_l + 1:) = d
      value = problem%value(parameters)
      if (.not. ieee_is_finite(value) .or. value >= huge(1.0_dp) / 1000.0_dp) info = 5
   end function corphylo_objective

   function corphylo_problem_value(self, x) result(value)
      !! Evaluates the profiled `corphylo.LL` objective for optimizer parameters.
      class(corphylo_problem), intent(inout) :: self !! Standardized corphylo data and normalized tree covariance.
      real(dp), intent(in) :: x(:) !! Packed lower-Cholesky elements followed by d or logit-d parameters.
      real(dp) :: value
      real(dp), allocatable :: coefficients(:)
      real(dp), allocatable :: covariance_b(:, :)
      real(dp), allocatable :: c(:, :)
      real(dp), allocatable :: d(:)
      real(dp), allocatable :: r(:, :)
      real(dp), allocatable :: v(:, :)
      integer :: status

      call evaluate_corphylo(self, x, value, coefficients, covariance_b, r, d, c, v, status)
      if (status /= 0 .or. .not. ieee_is_finite(value)) value = huge(1.0_dp) / 100.0_dp
   end function corphylo_problem_value

   subroutine corphylo_problem_gradient(self, x, gradient)
      !! Computes a robust central finite-difference gradient of the corphylo objective.
      class(corphylo_problem), intent(inout) :: self !! Corphylo objective differentiated numerically.
      real(dp), intent(in) :: x(:) !! Packed model parameter vector.
      real(dp), intent(out) :: gradient(:) !! Numerical derivative for each packed parameter.
      real(dp), allocatable :: xm(:)
      real(dp), allocatable :: xp(:)
      real(dp) :: f0
      real(dp) :: fm
      real(dp) :: fp
      real(dp) :: h
      integer :: i

      f0 = self%value(x)
      xm = x
      xp = x
      do i = 1, size(x)
         h = 1.0e-6_dp * max(1.0_dp, abs(x(i)))
         xm = x
         xp = x
         xm(i) = x(i) - h
         xp(i) = x(i) + h
         fm = self%value(xm)
         fp = self%value(xp)
         if (fm < huge(1.0_dp) / 1000.0_dp .and. fp < huge(1.0_dp) / 1000.0_dp) then
            gradient(i) = (fp - fm) / (2.0_dp * h)
         else if (fp < huge(1.0_dp) / 1000.0_dp) then
            gradient(i) = (fp - f0) / h
         else if (fm < huge(1.0_dp) / 1000.0_dp) then
            gradient(i) = (f0 - fm) / h
         else
            gradient(i) = 0.0_dp
         end if
      end do
   end subroutine corphylo_problem_gradient

   subroutine setup_corphylo(x, tree, measurement_error, covariates, problem, x_mean, x_sd, coefficient_scale, info)
      !! Standardizes inputs and constructs XX, UU, MM, tau, and normalized Vphy as in upstream corphylo.
      real(dp), intent(in) :: x(:, :) !! Species-by-trait raw observations.
      type(phylo_tree), intent(in) :: tree !! Rooted tree with one tip per species row.
      real(dp), intent(in), optional :: measurement_error(:, :) !! Raw measurement standard errors matching `x`.
      real(dp), intent(in), optional :: covariates(:, :, :) !! Raw covariates with shape `(n_species,n_predictor,n_trait)`.
      type(corphylo_problem), intent(out) :: problem !! Standardized data and matrices used by the likelihood.
      real(dp), allocatable, intent(out) :: x_mean(:) !! Raw trait means used for coefficient back-transformation.
      real(dp), allocatable, intent(out) :: x_sd(:) !! Raw trait sample standard deviations.
      real(dp), allocatable, intent(out) :: coefficient_scale(:) !! Multipliers for coefficient covariance back-transformation.
      integer, intent(out) :: info !! Zero on success or nonzero for invalid dimensions or singular Vphy.
      real(dp), allocatable :: inverse_vphy(:, :)
      real(dp), allocatable :: sem(:, :)
      real(dp), allocatable :: standardized(:, :)
      real(dp), allocatable :: u_mean(:, :)
      real(dp), allocatable :: u_sd(:, :)
      real(dp), allocatable :: vphy(:, :)
      real(dp) :: logdet
      integer :: active
      integer :: col
      integer :: i
      integer :: j
      integer :: n
      integer :: p
      integer :: q
      integer :: status
      integer :: trait

      info = 0
      n = size(x, 1)
      p = size(x, 2)
      if (n < 2 .or. p < 1 .or. tree%n_tip /= n .or. .not. tree%valid() .or. .not. tree%has_lengths()) then
         info = 1
         return
      end if
      if (any(tree%edge_length < 0.0_dp) .or. any(.not. ieee_is_finite(x))) then
         info = 2
         return
      end if
      allocate(x_mean(p), x_sd(p), standardized(n, p))
      do trait = 1, p
         x_mean(trait) = sum(x(:, trait)) / real(n, dp)
         x_sd(trait) = sample_sd(x(:, trait))
         if (x_sd(trait) <= 0.0_dp .or. .not. ieee_is_finite(x_sd(trait))) then
            info = 3
            return
         end if
         standardized(:, trait) = (x(:, trait) - x_mean(trait)) / x_sd(trait)
      end do
      allocate(sem(n, p))
      sem = 0.0_dp
      if (present(measurement_error)) then
         if (size(measurement_error, 1) /= n .or. size(measurement_error, 2) /= p) then
            info = 4
            return
         end if
         if (any(measurement_error < 0.0_dp) .or. any(.not. ieee_is_finite(measurement_error))) then
            info = 5
            return
         end if
         do trait = 1, p
            sem(:, trait) = measurement_error(:, trait) / x_sd(trait)
         end do
      end if
      problem%n = n
      problem%p = p
      problem%xx = reshape(standardized, [n * p])
      problem%mm = reshape(sem**2, [n * p])

      q = 0
      active = 0
      if (present(covariates)) then
         if (size(covariates, 1) /= n .or. size(covariates, 3) /= p) then
            info = 6
            return
         end if
         q = size(covariates, 2)
         allocate(u_mean(q, p), u_sd(q, p))
         do trait = 1, p
            do j = 1, q
               if (any(.not. ieee_is_finite(covariates(:, j, trait)))) then
                  info = 7
                  return
               end if
               u_mean(j, trait) = sum(covariates(:, j, trait)) / real(n, dp)
               u_sd(j, trait) = sample_sd(covariates(:, j, trait))
               if (u_sd(j, trait) > 0.0_dp) active = active + 1
            end do
         end do
      else
         allocate(u_mean(0, p), u_sd(0, p))
      end if
      allocate(problem%uu(n * p, p + active), coefficient_scale(p + active))
      problem%uu = 0.0_dp
      do trait = 1, p
         problem%uu((trait - 1) * n + 1:trait * n, trait) = 1.0_dp
         coefficient_scale(trait) = x_sd(trait)
      end do
      col = p
      if (present(covariates)) then
         do trait = 1, p
            do j = 1, q
               if (u_sd(j, trait) <= 0.0_dp) cycle
               col = col + 1
               problem%uu((trait - 1) * n + 1:trait * n, col) = &
                  (covariates(:, j, trait) - u_mean(j, trait)) / u_sd(j, trait)
               coefficient_scale(col) = x_sd(trait) / u_sd(j, trait)
            end do
         end do
      end if

      call phylogenetic_vcv(tree, vphy, status, correlation=.false.)
      if (status /= 0 .or. maxval(vphy) <= 0.0_dp) then
         info = 10 + status
         return
      end if
      vphy = vphy / maxval(vphy)
      call spd_inverse_logdet(vphy, inverse_vphy, logdet, status)
      if (status /= 0) then
         info = 20 + status
         return
      end if
      vphy = vphy / exp(logdet / real(n, dp))
      problem%vphy = vphy
      allocate(problem%tau(n, n))
      do i = 1, n
         do j = 1, n
            problem%tau(i, j) = vphy(j, j) - vphy(i, j)
         end do
      end do
   end subroutine setup_corphylo

   subroutine evaluate_corphylo(problem, parameters, objective, coefficients, covariance_b, r, d, c, v, info)
      !! Evaluates the profiled multivariate OU covariance likelihood and regression coefficients.
      type(corphylo_problem), intent(in) :: problem !! Standardized data and normalized phylogenetic covariance.
      real(dp), intent(in) :: parameters(:) !! Packed lower-Cholesky elements followed by d/logit-d values.
      real(dp), intent(out) :: objective !! Upstream `corphylo.LL` value excluding Gaussian constants.
      real(dp), allocatable, intent(out) :: coefficients(:) !! Profiled standardized regression coefficients.
      real(dp), allocatable, intent(out) :: covariance_b(:, :) !! Conditional covariance of standardized coefficients.
      real(dp), allocatable, intent(out) :: r(:, :) !! Trait process covariance `transpose(L)*L`.
      real(dp), allocatable, intent(out) :: d(:) !! Positive OU d values after optional logistic transformation.
      real(dp), allocatable, intent(out) :: c(:, :) !! Stacked process covariance before measurement error.
      real(dp), allocatable, intent(out) :: v(:, :) !! Full covariance including measurement-error variances.
      integer, intent(out) :: info !! Zero on success or nonzero for invalid parameters or singular matrices.
      real(dp), allocatable :: inverse_normal(:, :)
      real(dp), allocatable :: inverse_v(:, :)
      real(dp), allocatable :: l(:, :)
      real(dp), allocatable :: normal(:, :)
      real(dp), allocatable :: residual(:)
      real(dp), allocatable :: rhs(:)
      real(dp) :: denom
      real(dp) :: logdet_normal
      real(dp) :: logdet_v
      real(dp) :: product_d
      real(dp) :: quadratic
      integer :: a
      integer :: b
      integer :: i
      integer :: j
      integer :: n_l
      integer :: status

      objective = huge(1.0_dp) / 100.0_dp
      info = 0
      n_l = problem%p * (problem%p + 1) / 2
      if (size(parameters) /= n_l + problem%p) then
         info = 1
         return
      end if
      allocate(l(problem%p, problem%p), d(problem%p))
      call unpack_lower(parameters(1:n_l), l)
      r = matmul(transpose(l), l)
      if (problem%constrain_d) then
         if (maxval(abs(parameters(n_l + 1:))) > 10.0_dp) then
            info = 2
            return
         end if
         d = 1.0_dp / (1.0_dp + exp(-parameters(n_l + 1:)))
      else
         d = parameters(n_l + 1:)
         if (any(d <= 0.0_dp) .or. maxval(d) > 10.0_dp) then
            info = 2
            return
         end if
      end if
      allocate(c(problem%n * problem%p, problem%n * problem%p))
      c = 0.0_dp
      do a = 1, problem%p
         do b = 1, problem%p
            product_d = d(a) * d(b)
            denom = 1.0_dp - product_d
            if (product_d <= 0.0_dp .or. abs(denom) <= 100.0_dp * epsilon(1.0_dp)) then
               info = 3
               return
            end if
            do i = 1, problem%n
               do j = 1, problem%n
                  c((a - 1) * problem%n + i, (b - 1) * problem%n + j) = r(a, b) * &
                     exp(log(d(a)) * problem%tau(i, j) + log(d(b)) * problem%tau(j, i)) * &
                     (1.0_dp - exp(log(product_d) * problem%vphy(i, j))) / denom
               end do
            end do
         end do
      end do
      v = c
      do i = 1, size(problem%mm)
         v(i, i) = v(i, i) + problem%mm(i)
      end do
      call spd_inverse_logdet(v, inverse_v, logdet_v, status)
      if (status /= 0) then
         info = 4
         return
      end if
      normal = matmul(transpose(problem%uu), matmul(inverse_v, problem%uu))
      call spd_inverse_logdet(normal, inverse_normal, logdet_normal, status)
      if (status /= 0) then
         info = 5
         return
      end if
      rhs = matmul(transpose(problem%uu), matmul(inverse_v, problem%xx))
      coefficients = matmul(inverse_normal, rhs)
      covariance_b = inverse_normal
      residual = problem%xx - matmul(problem%uu, coefficients)
      quadratic = dot_product(residual, matmul(inverse_v, residual))
      objective = 0.5_dp * (logdet_v + quadratic)
      if (problem%reml) objective = objective + 0.5_dp * logdet_normal
   end subroutine evaluate_corphylo

   subroutine corphylo_finish(problem, parameters, x_mean, coefficient_scale, result, covariance_b, info)
      !! Converts a fitted standardized corphylo solution to the principal upstream result quantities.
      type(corphylo_problem), intent(in) :: problem !! Standardized fitted problem definition.
      real(dp), intent(in) :: parameters(:) !! Optimized packed covariance and OU parameters.
      real(dp), intent(in) :: x_mean(:) !! Raw trait means for intercept back-transformation.
      real(dp), intent(in) :: coefficient_scale(:) !! Raw-scale multiplier for each design coefficient.
      type(corphylo_result), intent(inout) :: result !! Result object populated by this routine.
      real(dp), allocatable, intent(out) :: covariance_b(:, :) !! Standardized coefficient covariance from likelihood profile.
      integer, intent(out) :: info !! Zero on success or nonzero for singular fitted covariance.
      real(dp), allocatable :: coefficients(:)
      real(dp), allocatable :: c(:, :)
      real(dp), allocatable :: d(:)
      real(dp), allocatable :: r(:, :)
      real(dp), allocatable :: v(:, :)
      real(dp) :: objective
      real(dp) :: txx
      real(dp) :: variance
      integer :: i
      integer :: k
      integer :: status

      call evaluate_corphylo(problem, parameters, objective, coefficients, covariance_b, r, d, c, v, status)
      if (status /= 0) then
         info = status
         return
      end if
      info = 0
      allocate(result%correlation(problem%p, problem%p))
      do i = 1, problem%p
         if (r(i, i) <= 0.0_dp) then
            info = 6
            return
         end if
      end do
      do i = 1, problem%p
         result%correlation(i, :) = r(i, :) / sqrt(r(i, i) * diagonal_vector(r))
      end do
      result%d = d
      result%coefficients = coefficients
      do i = 1, problem%p
         result%coefficients(i) = result%coefficients(i) + x_mean(i)
      end do
      do i = problem%p + 1, size(result%coefficients)
         result%coefficients(i) = result%coefficients(i) * coefficient_scale(i)
      end do
      allocate(result%coefficient_covariance(size(coefficients), size(coefficients)))
      do i = 1, size(coefficients)
         result%coefficient_covariance(i, :) = coefficient_scale(i) * covariance_b(i, :) * coefficient_scale
      end do
      allocate(result%coefficient_se(size(coefficients)))
      do i = 1, size(coefficients)
         variance = result%coefficient_covariance(i, i)
         result%coefficient_se(i) = sqrt(max(0.0_dp, variance))
      end do
      result%process_covariance = r
      result%observation_covariance = v
      result%phylogenetic_covariance = c
      if (problem%reml) then
         txx = dot_product(problem%xx, problem%xx)
         if (txx <= 0.0_dp) then
            info = 7
            return
         end if
         result%log_likelihood = -0.5_dp * real(problem%n * problem%p - size(problem%uu, 2), dp) * log(two_pi) + &
            0.5_dp * log(txx) - objective
      else
         result%log_likelihood = -0.5_dp * real(problem%n * problem%p, dp) * log(two_pi) - objective
      end if
      k = size(parameters) + size(problem%uu, 2)
      result%aic = -2.0_dp * result%log_likelihood + 2.0_dp * real(k, dp)
      result%bic = -2.0_dp * result%log_likelihood + real(k, dp) * (log(real(problem%n, dp)) - log(pi))
   end subroutine corphylo_finish

   subroutine initial_cholesky(xx, n, p, lower, info)
      !! Forms a lower Cholesky start from the standardized trait sample covariance.
      real(dp), intent(in) :: xx(:) !! Column-major standardized trait data with length `n*p`.
      integer, intent(in) :: n !! Number of species rows.
      integer, intent(in) :: p !! Number of trait columns.
      real(dp), allocatable, intent(out) :: lower(:, :) !! Lower Cholesky factor of the standardized trait covariance.
      integer, intent(out) :: info !! Zero on success or nonzero if the sample covariance is not positive definite.
      real(dp), allocatable :: covariance(:, :)
      real(dp), allocatable :: data(:, :)
      real(dp), allocatable :: centered(:, :)
      real(dp) :: s
      integer :: i
      integer :: j
      integer :: k

      data = reshape(xx, [n, p])
      centered = data
      allocate(covariance(p, p), lower(p, p))
      covariance = matmul(transpose(centered), centered) / real(n - 1, dp)
      lower = 0.0_dp
      info = 0
      do i = 1, p
         do j = 1, i
            s = covariance(i, j)
            do k = 1, j - 1
               s = s - lower(i, k) * lower(j, k)
            end do
            if (i == j) then
               if (s <= 0.0_dp) then
                  info = 1
                  return
               end if
               lower(i, j) = sqrt(s)
            else
               lower(i, j) = s / lower(j, j)
            end if
         end do
      end do
   end subroutine initial_cholesky

   pure subroutine pack_lower(lower, packed)
      !! Packs a lower-triangular matrix in R's column-major `lower.tri(...,diag=TRUE)` order.
      real(dp), intent(in) :: lower(:, :) !! Square lower-triangular matrix.
      real(dp), intent(out) :: packed(:) !! Packed lower triangle with length `p*(p+1)/2`.
      integer :: col
      integer :: index
      integer :: row

      index = 0
      do col = 1, size(lower, 2)
         do row = col, size(lower, 1)
            index = index + 1
            packed(index) = lower(row, col)
         end do
      end do
   end subroutine pack_lower

   pure subroutine unpack_lower(packed, lower)
      !! Unpacks R column-major lower-triangle elements to a square lower matrix.
      real(dp), intent(in) :: packed(:) !! Packed lower-triangle elements.
      real(dp), intent(out) :: lower(:, :) !! Square lower-triangular matrix receiving the elements.
      integer :: col
      integer :: index
      integer :: row

      lower = 0.0_dp
      index = 0
      do col = 1, size(lower, 2)
         do row = col, size(lower, 1)
            index = index + 1
            lower(row, col) = packed(index)
         end do
      end do
   end subroutine unpack_lower

   pure real(dp) function sample_sd(values) result(sd)
      !! Computes R-compatible sample standard deviation with denominator `n-1`.
      real(dp), intent(in) :: values(:) !! Finite sample values with at least two entries.
      real(dp) :: mean

      mean = sum(values) / real(size(values), dp)
      sd = sqrt(sum((values - mean)**2) / real(size(values) - 1, dp))
   end function sample_sd

   pure function diagonal_vector(matrix) result(diagonal)
      !! Extracts a square matrix's main diagonal.
      real(dp), intent(in) :: matrix(:, :) !! Square matrix whose diagonal is requested.
      real(dp) :: diagonal(min(size(matrix, 1), size(matrix, 2)))
      integer :: i

      do i = 1, size(diagonal)
         diagonal(i) = matrix(i, i)
      end do
   end function diagonal_vector

end module ape_corphylo
