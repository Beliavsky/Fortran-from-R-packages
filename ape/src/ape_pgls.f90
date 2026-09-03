! SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only
!
! Phylogenetic correlation matrices translated from ape R/PGLS.R
! (Copyright 2004-2021 Julien Dutheil and 2006-2017 Emmanuel Paradis),
! with deterministic GLS/profile-likelihood infrastructure replacing nlme plumbing.
module ape_pgls
   use r_kinds, only : dp
   use r_linalg, only : solve_spd, spd_inverse_logdet
   use ape_types, only : phylo_tree
   use ape_topology, only : phylogenetic_vcv
   use ape_tree_algorithms, only : dist_nodes, node_depth_count
   use ape_optimize, only : bounded_problem, bounded_bfgs
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   implicit none
   private

   real(dp), parameter :: two_pi = 2.0_dp * acos(-1.0_dp)

   type, public :: pgls_result
      real(dp), allocatable :: coefficients(:)
      real(dp), allocatable :: standard_error(:)
      real(dp) :: sigma2 = 0.0_dp
      real(dp) :: log_likelihood = -huge(1.0_dp)
      real(dp) :: residual_sum = huge(1.0_dp)
      real(dp) :: correlation_parameter = 0.0_dp
      integer :: iterations = 0
      logical :: converged = .false.
   end type pgls_result

   type, extends(bounded_problem) :: pgls_profile_problem
      type(phylo_tree) :: tree
      real(dp), allocatable :: y(:)
      real(dp), allocatable :: design(:, :)
      character(len=12) :: model = 'brownian'
      logical :: reml = .false.
   contains
      procedure :: value => pgls_profile_value
      procedure :: gradient => pgls_profile_gradient
   end type pgls_profile_problem

   public :: cor_brownian
   public :: cor_martins
   public :: grafen_tree
   public :: cor_grafen
   public :: cor_pagel
   public :: cor_blomberg
   public :: pgls_fit
   public :: pgls_fit_model

contains

   subroutine cor_brownian(tree, matrix, info, correlation)
      !! Returns ape's Brownian phylogenetic covariance or correlation matrix.
      type(phylo_tree), intent(in) :: tree !! Rooted tree with branch lengths defining shared evolutionary history.
      real(dp), allocatable, intent(out) :: matrix(:, :) !! Tip covariance/correlation matrix in numeric tip order.
      integer, intent(out) :: info !! Zero on success or the underlying phylogenetic-VCV status code.
      logical, intent(in), optional :: correlation !! Normalize to unit diagonal when true; default true as `corMatrix` does.
      logical :: as_correlation

      as_correlation = .true.
      if (present(correlation)) as_correlation = correlation
      call phylogenetic_vcv(tree, matrix, info, as_correlation)
   end subroutine cor_brownian

   subroutine cor_martins(tree, alpha, matrix, info)
      !! Returns the Martins exponential phylogenetic correlation matrix `exp(-alpha * cophenetic)`.
      type(phylo_tree), intent(in) :: tree !! Rooted tree with branch lengths used for tip-to-tip distances.
      real(dp), intent(in) :: alpha !! Nonnegative exponential decay parameter.
      real(dp), allocatable, intent(out) :: matrix(:, :) !! Unit-diagonal tip correlation matrix.
      integer, intent(out) :: info !! Zero on success; nonzero for invalid alpha or node-distance calculation.
      real(dp), allocatable :: distance(:, :)
      integer :: status

      info = 0
      if (.not. ieee_is_finite(alpha) .or. alpha < 0.0_dp) then
         allocate(matrix(0, 0))
         info = 1
         return
      end if
      call dist_nodes(tree, distance, status)
      if (status /= 0) then
         allocate(matrix(0, 0))
         info = 10 + status
         return
      end if
      allocate(matrix(tree%n_tip, tree%n_tip))
      matrix = exp(-alpha * distance(1:tree%n_tip, 1:tree%n_tip))
   end subroutine cor_martins

   subroutine grafen_tree(tree, rho, transformed, info)
      !! Applies ape's Grafen branch-length transformation at power `rho`.
      type(phylo_tree), intent(in) :: tree !! Rooted tree topology to transform; input branch lengths are ignored.
      real(dp), intent(in) :: rho !! Nonnegative Grafen power parameter.
      type(phylo_tree), intent(out) :: transformed !! Tree copy with Grafen edge lengths.
      integer, intent(out) :: info !! Zero on success; nonzero for invalid tree, tip count, or power.
      integer, allocatable :: depth(:)
      real(dp), allocatable :: height(:)
      real(dp) :: scale
      integer :: e

      transformed = tree
      info = 0
      if (.not. tree%valid() .or. tree%n_tip < 2) then
         info = 1
         return
      end if
      if (.not. ieee_is_finite(rho) .or. rho < 0.0_dp) then
         info = 2
         return
      end if
      call node_depth_count(tree, 1, depth, info)
      if (info /= 0) return
      allocate(height(tree%total_nodes()))
      scale = real(tree%n_tip - 1, dp)
      height = (real(depth - 1, dp) / scale)**rho
      if (.not. allocated(transformed%edge_length)) allocate(transformed%edge_length(tree%nedge()))
      do e = 1, tree%nedge()
         transformed%edge_length(e) = height(tree%edge(e, 1)) - height(tree%edge(e, 2))
      end do
   end subroutine grafen_tree

   subroutine cor_grafen(tree, rho, matrix, info, correlation)
      !! Returns ape's Grafen-transformed Brownian covariance or correlation matrix.
      type(phylo_tree), intent(in) :: tree !! Rooted tree topology used by the Grafen transformation.
      real(dp), intent(in) :: rho !! Nonnegative Grafen branch-length power.
      real(dp), allocatable, intent(out) :: matrix(:, :) !! Tip covariance/correlation matrix after transformation.
      integer, intent(out) :: info !! Zero on success or a transformation/VCV status code.
      logical, intent(in), optional :: correlation !! Normalize to unit diagonal when true; default true.
      type(phylo_tree) :: transformed
      logical :: as_correlation
      integer :: status

      as_correlation = .true.
      if (present(correlation)) as_correlation = correlation
      call grafen_tree(tree, rho, transformed, status)
      if (status /= 0) then
         allocate(matrix(0, 0))
         info = status
         return
      end if
      call phylogenetic_vcv(transformed, matrix, info, as_correlation)
   end subroutine cor_grafen

   subroutine cor_pagel(tree, lambda, matrix, info, correlation)
      !! Returns ape's Pagel-lambda Brownian covariance/correlation matrix.
      type(phylo_tree), intent(in) :: tree !! Rooted tree with branch lengths defining the Brownian base matrix.
      real(dp), intent(in) :: lambda !! Pagel lambda in the inclusive interval zero through one.
      real(dp), allocatable, intent(out) :: matrix(:, :) !! Matrix with off-diagonals multiplied by lambda.
      integer, intent(out) :: info !! Zero on success; nonzero for invalid lambda or VCV construction.
      logical, intent(in), optional :: correlation !! Normalize the Brownian base matrix first; default true.
      real(dp), allocatable :: diagonal(:)
      logical :: as_correlation
      integer :: i

      info = 0
      if (.not. ieee_is_finite(lambda) .or. lambda < 0.0_dp .or. lambda > 1.0_dp) then
         allocate(matrix(0, 0))
         info = 1
         return
      end if
      as_correlation = .true.
      if (present(correlation)) as_correlation = correlation
      call phylogenetic_vcv(tree, matrix, info, as_correlation)
      if (info /= 0) return
      allocate(diagonal(tree%n_tip))
      do i = 1, tree%n_tip
         diagonal(i) = matrix(i, i)
      end do
      matrix = lambda * matrix
      do i = 1, tree%n_tip
         matrix(i, i) = diagonal(i)
      end do
   end subroutine cor_pagel

   subroutine cor_blomberg(tree, g, matrix, info, correlation)
      !! Returns ape's Blomberg-g transformed Brownian covariance or correlation matrix.
      type(phylo_tree), intent(in) :: tree !! Rooted tree whose root distances are transformed.
      real(dp), intent(in) :: g !! Strictly positive Blomberg g parameter.
      real(dp), allocatable, intent(out) :: matrix(:, :) !! Tip covariance/correlation matrix after branch transformation.
      integer, intent(out) :: info !! Zero on success; nonzero for invalid g, distances, or VCV construction.
      logical, intent(in), optional :: correlation !! Normalize to unit diagonal when true; default true.
      type(phylo_tree) :: transformed
      real(dp), allocatable :: distance(:, :)
      real(dp), allocatable :: root_distance(:)
      logical :: as_correlation
      integer :: e
      integer :: root
      integer :: status

      info = 0
      if (.not. ieee_is_finite(g) .or. g <= 0.0_dp) then
         allocate(matrix(0, 0))
         info = 1
         return
      end if
      call dist_nodes(tree, distance, status)
      if (status /= 0) then
         allocate(matrix(0, 0))
         info = 10 + status
         return
      end if
      root = tree%root()
      allocate(root_distance(tree%total_nodes()))
      root_distance = distance(root, :)**(1.0_dp / g)
      transformed = tree
      if (.not. allocated(transformed%edge_length)) allocate(transformed%edge_length(tree%nedge()))
      do e = 1, tree%nedge()
         transformed%edge_length(e) = root_distance(tree%edge(e, 2)) - root_distance(tree%edge(e, 1))
      end do
      as_correlation = .true.
      if (present(correlation)) as_correlation = correlation
      call phylogenetic_vcv(transformed, matrix, info, as_correlation)
   end subroutine cor_blomberg

   subroutine pgls_fit(y, design, covariance, result, info, reml)
      !! Fits generalized least squares with a supplied positive-definite phylogenetic covariance matrix.
      real(dp), intent(in) :: y(:) !! Numeric response vector with one observation per covariance row.
      real(dp), intent(in) :: design(:, :) !! Full-rank design matrix with observations in rows and predictors in columns.
      real(dp), intent(in) :: covariance(:, :) !! Positive-definite covariance/correlation matrix for the observations.
      type(pgls_result), intent(out) :: result !! GLS coefficients, SEs, profiled variance, residual sum, and log likelihood.
      integer, intent(out) :: info !! Zero on success or nonzero for dimensions, rank, or linear-algebra failure.
      logical, intent(in), optional :: reml !! Use the residual maximum likelihood profile when true; default false.
      real(dp), allocatable :: covariance_inverse(:, :)
      real(dp), allocatable :: coefficient_covariance(:, :)
      real(dp), allocatable :: normal_matrix(:, :)
      real(dp), allocatable :: normal_inverse(:, :)
      real(dp), allocatable :: rhs(:)
      real(dp), allocatable :: residual(:)
      real(dp), allocatable :: weighted_design(:, :)
      real(dp), allocatable :: weighted_y(:)
      real(dp) :: covariance_logdet
      real(dp) :: normal_logdet
      real(dp) :: denominator
      integer :: i
      integer :: n
      integer :: p
      integer :: status
      logical :: use_reml

      result = pgls_result()
      info = 0
      n = size(y)
      p = size(design, 2)
      if (n < 1 .or. p < 1 .or. size(design, 1) /= n) then
         info = 1
         return
      end if
      if (size(covariance, 1) /= n .or. size(covariance, 2) /= n .or. p >= n) then
         info = 2
         return
      end if
      if (.not. all(ieee_is_finite(y)) .or. .not. all(ieee_is_finite(design)) &
         .or. .not. all(ieee_is_finite(covariance))) then
         info = 3
         return
      end if
      call spd_inverse_logdet(covariance, covariance_inverse, covariance_logdet, status)
      if (status /= 0) then
         info = 10 + status
         return
      end if
      allocate(weighted_design(n, p), weighted_y(n), normal_matrix(p, p), rhs(p))
      weighted_design = matmul(covariance_inverse, design)
      weighted_y = matmul(covariance_inverse, y)
      normal_matrix = matmul(transpose(design), weighted_design)
      rhs = matmul(transpose(design), weighted_y)
      allocate(result%coefficients(p), result%standard_error(p))
      call solve_spd(normal_matrix, rhs, result%coefficients, status)
      if (status /= 0) then
         info = 20 + status
         return
      end if
      allocate(residual(n))
      residual = y - matmul(design, result%coefficients)
      result%residual_sum = dot_product(residual, matmul(covariance_inverse, residual))
      use_reml = .false.
      if (present(reml)) use_reml = reml
      if (use_reml) then
         denominator = real(n - p, dp)
      else
         denominator = real(n, dp)
      end if
      if (result%residual_sum <= 0.0_dp) then
         info = 30
         return
      end if
      result%sigma2 = result%residual_sum / denominator
      call spd_inverse_logdet(normal_matrix, normal_inverse, normal_logdet, status)
      if (status /= 0) then
         info = 40 + status
         return
      end if
      allocate(coefficient_covariance(p, p))
      coefficient_covariance = result%sigma2 * normal_inverse
      do i = 1, p
         result%standard_error(i) = sqrt(max(0.0_dp, coefficient_covariance(i, i)))
      end do
      if (use_reml) then
         result%log_likelihood = -0.5_dp * (denominator * (log(two_pi * result%sigma2) + 1.0_dp) &
            + covariance_logdet + normal_logdet)
      else
         result%log_likelihood = -0.5_dp * (real(n, dp) * (log(two_pi * result%sigma2) + 1.0_dp) &
            + covariance_logdet)
      end if
      result%converged = .true.
   end subroutine pgls_fit

   subroutine pgls_fit_model(tree, y, design, model, result, info, parameter, reml, lower_bound, upper_bound, &
      max_iter, tolerance)
      !! Fits GLS under an ape phylogenetic correlation model, profiling one correlation parameter when absent.
      type(phylo_tree), intent(in) :: tree !! Rooted tree defining the phylogenetic correlation model.
      real(dp), intent(in) :: y(:) !! Response vector in numeric tip order.
      real(dp), intent(in) :: design(:, :) !! Design matrix whose rows correspond to numeric tip order.
      character(len=*), intent(in) :: model !! `brownian`, `martins`, `grafen`, `pagel`, or `blomberg`.
      type(pgls_result), intent(out) :: result !! Fitted GLS result and estimated/fixed correlation parameter.
      integer, intent(out) :: info !! Zero on success or nonzero for model, optimization, or GLS failure.
      real(dp), intent(in), optional :: parameter !! Fixed model parameter; omitted to profile non-Brownian models.
      logical, intent(in), optional :: reml !! Use REML rather than ML when true; default false.
      real(dp), intent(in), optional :: lower_bound !! Optional lower bound for a profiled correlation parameter.
      real(dp), intent(in), optional :: upper_bound !! Optional upper bound for a profiled correlation parameter.
      integer, intent(in), optional :: max_iter !! Maximum projected-BFGS iterations for parameter profiling.
      real(dp), intent(in), optional :: tolerance !! Optimizer tolerance; default `1e-8`.
      type(pgls_profile_problem) :: problem
      real(dp), allocatable :: covariance(:, :)
      real(dp) :: lower(1)
      real(dp) :: objective
      real(dp) :: theta(1)
      real(dp) :: upper(1)
      integer :: iterations
      integer :: opt_info
      integer :: status
      logical :: use_reml
      character(len=12) :: key

      result = pgls_result()
      info = 0
      key = model_key(model)
      use_reml = .false.
      if (present(reml)) use_reml = reml
      if (key == 'invalid') then
         info = 1
         return
      end if
      if (key == 'brownian') then
         call cor_brownian(tree, covariance, status, correlation=.true.)
         if (status /= 0) then
            info = 10 + status
            return
         end if
         call pgls_fit(y, design, covariance, result, info, use_reml)
         return
      end if
      if (present(parameter)) then
         call correlation_matrix(tree, key, parameter, covariance, status)
         if (status /= 0) then
            info = 20 + status
            return
         end if
         call pgls_fit(y, design, covariance, result, info, use_reml)
         if (info == 0) result%correlation_parameter = parameter
         return
      end if

      call default_parameter_bounds(key, lower(1), upper(1), theta(1))
      if (present(lower_bound)) lower(1) = lower_bound
      if (present(upper_bound)) upper(1) = upper_bound
      if (.not. ieee_is_finite(lower(1)) .or. .not. ieee_is_finite(upper(1)) .or. lower(1) >= upper(1)) then
         info = 2
         return
      end if
      theta(1) = max(lower(1), min(upper(1), theta(1)))
      problem%tree = tree
      problem%y = y
      problem%design = design
      problem%model = key
      problem%reml = use_reml
      call bounded_bfgs(problem, theta, lower, upper, objective, opt_info, iterations, max_iter, tolerance)
      if (opt_info /= 0 .and. opt_info /= 7) then
         info = 30 + opt_info
         return
      end if
      call correlation_matrix(tree, key, theta(1), covariance, status)
      if (status /= 0) then
         info = 40 + status
         return
      end if
      call pgls_fit(y, design, covariance, result, status, use_reml)
      if (status /= 0) then
         info = 50 + status
         return
      end if
      result%correlation_parameter = theta(1)
      result%iterations = iterations
      result%converged = opt_info == 0
      info = 0
   end subroutine pgls_fit_model

   function pgls_profile_value(self, x) result(value)
      class(pgls_profile_problem), intent(inout) :: self !! Profile-likelihood problem containing tree, data, and model.
      real(dp), intent(in) :: x(:) !! Candidate one-element phylogenetic correlation-parameter vector.
      real(dp) :: value
      type(pgls_result) :: fit
      real(dp), allocatable :: covariance(:, :)
      integer :: status

      value = huge(1.0_dp) / 100.0_dp
      if (size(x) /= 1) return
      call correlation_matrix(self%tree, self%model, x(1), covariance, status)
      if (status /= 0) return
      call pgls_fit(self%y, self%design, covariance, fit, status, self%reml)
      if (status /= 0 .or. .not. ieee_is_finite(fit%log_likelihood)) return
      value = -2.0_dp * fit%log_likelihood
   end function pgls_profile_value

   subroutine pgls_profile_gradient(self, x, gradient)
      class(pgls_profile_problem), intent(inout) :: self !! Profile-likelihood problem used for finite differences.
      real(dp), intent(in) :: x(:) !! Candidate one-element correlation-parameter vector.
      real(dp), intent(out) :: gradient(:) !! Finite-difference gradient of minus twice the profiled log likelihood.
      real(dp) :: f0
      real(dp) :: fm
      real(dp) :: fp
      real(dp) :: h
      real(dp) :: xm(1)
      real(dp) :: xp(1)

      h = epsilon(1.0_dp)**(1.0_dp / 3.0_dp) * max(1.0_dp, abs(x(1)))
      f0 = self%value(x)
      xp = x
      xp(1) = x(1) + h
      fp = self%value(xp)
      if (x(1) > h) then
         xm = x
         xm(1) = x(1) - h
         fm = self%value(xm)
         gradient(1) = (fp - fm) / (2.0_dp * h)
      else
         gradient(1) = (fp - f0) / h
      end if
   end subroutine pgls_profile_gradient

   subroutine correlation_matrix(tree, model, parameter, matrix, info)
      !! Dispatches one named ape phylogenetic correlation structure at a fixed parameter.
      type(phylo_tree), intent(in) :: tree !! Rooted tree defining the requested correlation structure.
      character(len=*), intent(in) :: model !! Canonical lowercase model key.
      real(dp), intent(in) :: parameter !! Fixed correlation-model parameter.
      real(dp), allocatable, intent(out) :: matrix(:, :) !! Unit-diagonal tip correlation matrix.
      integer, intent(out) :: info !! Zero on success or a model-specific status code.

      select case (trim(model))
      case ('martins')
         call cor_martins(tree, parameter, matrix, info)
      case ('grafen')
         call cor_grafen(tree, parameter, matrix, info, correlation=.true.)
      case ('pagel')
         call cor_pagel(tree, parameter, matrix, info, correlation=.true.)
      case ('blomberg')
         call cor_blomberg(tree, parameter, matrix, info, correlation=.true.)
      case default
         allocate(matrix(0, 0))
         info = 1
      end select
   end subroutine correlation_matrix

   pure subroutine default_parameter_bounds(model, lower, upper, start)
      !! Supplies conservative finite bounds for deterministic profiling of one ape correlation parameter.
      character(len=*), intent(in) :: model !! Canonical model key.
      real(dp), intent(out) :: lower !! Default inclusive lower bound.
      real(dp), intent(out) :: upper !! Default inclusive upper bound.
      real(dp), intent(out) :: start !! Default interior starting value.

      select case (trim(model))
      case ('pagel')
         lower = 0.0_dp
         upper = 1.0_dp
         start = 0.5_dp
      case ('martins')
         lower = 0.0_dp
         upper = 100.0_dp
         start = 0.1_dp
      case ('grafen')
         lower = 1.0e-6_dp
         upper = 100.0_dp
         start = 1.0_dp
      case ('blomberg')
         lower = 1.0e-6_dp
         upper = 100.0_dp
         start = 1.0_dp
      case default
         lower = 0.0_dp
         upper = 1.0_dp
         start = 0.5_dp
      end select
   end subroutine default_parameter_bounds

   pure function model_key(model) result(key)
      !! Normalizes supported PGLS model names to lowercase canonical keys.
      character(len=*), intent(in) :: model !! User model name, matched case-insensitively.
      character(len=12) :: key
      character(len=len(model)) :: lower
      integer :: code
      integer :: i

      lower = model
      do i = 1, len(model)
         code = iachar(lower(i:i))
         if (code >= iachar('A') .and. code <= iachar('Z')) lower(i:i) = achar(code + 32)
      end do
      select case (trim(adjustl(lower)))
      case ('brownian', 'bm')
         key = 'brownian'
      case ('martins')
         key = 'martins'
      case ('grafen')
         key = 'grafen'
      case ('pagel')
         key = 'pagel'
      case ('blomberg')
         key = 'blomberg'
      case default
         key = 'invalid'
      end select
   end function model_key

end module ape_pgls
