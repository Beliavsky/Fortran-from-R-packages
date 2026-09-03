! SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only
!
! Continuous Brownian-motion likelihood branches translated from ape R/ace.R.
! Upstream copyright and provenance are documented in NOTICE.md.
module ape_ace
   use r_kinds, only : dp
   use r_linalg, only : solve_spd, spd_inverse_logdet
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   use ape_types, only : phylo_tree, child_counts
   use ape_topology, only : phylogenetic_vcv
   use ape_pgls, only : pgls_result, pgls_fit, cor_brownian, cor_martins, cor_grafen, grafen_tree
   use ape_tree_algorithms, only : ace_pic, dist_nodes, node_depth_edgelength, mrca
   implicit none
   private

   type, public :: ace_continuous_result
      real(dp), allocatable :: estimates(:)
      real(dp), allocatable :: standard_error(:)
      real(dp) :: sigma2 = 0.0_dp
      real(dp) :: sigma2_standard_error = 0.0_dp
      real(dp) :: log_likelihood = 0.0_dp
      real(dp) :: residual_log_likelihood = 0.0_dp
   end type ace_continuous_result

   public :: ace_continuous_ml
   public :: ace_continuous_reml
   public :: ace_continuous_gls

contains

   subroutine ace_continuous_ml(tree, phenotype, result, info, kappa)
      !! Fits ape's continuous Brownian `ace(..., method="ML")` objective in closed form.
      type(phylo_tree), intent(in) :: tree !! Rooted fully dichotomous tree with positive branch lengths.
      real(dp), intent(in) :: phenotype(:) !! One continuous trait value per tip in numeric tip order.
      type(ace_continuous_result), intent(out) :: result !! ML ancestral states, standard errors, variance, and log-likelihood.
      integer, intent(out) :: info !! Zero on success; nonzero for invalid tree/data or failed linear algebra.
      real(dp), intent(in), optional :: kappa !! Positive or zero branch-length power applied before fitting; default one.
      type(phylo_tree) :: work
      real(dp), allocatable :: inverse_laplacian(:, :)
      real(dp), allocatable :: laplacian(:, :)
      real(dp), allocatable :: rhs(:)
      real(dp) :: logdet
      real(dp) :: residual_sum
      integer :: i
      integer :: m
      integer :: status

      result = ace_continuous_result()
      call prepare_continuous_problem(tree, phenotype, work, laplacian, rhs, info, kappa)
      if (info /= 0) return
      m = work%n_node
      allocate(result%estimates(m), result%standard_error(m))
      call solve_spd(laplacian, rhs, result%estimates, status)
      if (status /= 0) then
         info = 10 + status
         return
      end if
      residual_sum = edge_residual_sum(work, phenotype, result%estimates)
      result%sigma2 = residual_sum / (2.0_dp * real(m, dp))
      if (result%sigma2 <= 0.0_dp) then
         info = 20
         return
      end if
      result%log_likelihood = -residual_sum / (2.0_dp * result%sigma2) &
         - real(m, dp) * log(result%sigma2)
      result%sigma2_standard_error = result%sigma2 / sqrt(2.0_dp * real(m, dp))
      call spd_inverse_logdet(laplacian, inverse_laplacian, logdet, status)
      if (status /= 0) then
         info = 30 + status
         return
      end if
      do i = 1, m
         result%standard_error(i) = sqrt(0.5_dp * result%sigma2 * inverse_laplacian(i, i))
      end do
   end subroutine ace_continuous_ml

   subroutine ace_continuous_reml(tree, phenotype, result, info, kappa)
      !! Fits ape's continuous Brownian `ace(..., method="REML")` computational objective in closed form.
      type(phylo_tree), intent(in) :: tree !! Rooted fully dichotomous tree with positive branch lengths.
      real(dp), intent(in) :: phenotype(:) !! One continuous trait value per tip in numeric tip order.
      type(ace_continuous_result), intent(out) :: result !! REML variance and ancestral-state estimates/standard errors.
      integer, intent(out) :: info !! Zero on success; nonzero for invalid tree/data or failed linear algebra.
      real(dp), intent(in), optional :: kappa !! Positive or zero branch-length power applied before fitting; default one.
      type(phylo_tree) :: work
      real(dp), allocatable :: covariance(:, :)
      real(dp), allocatable :: inverse_laplacian(:, :)
      real(dp), allocatable :: laplacian(:, :)
      real(dp), allocatable :: pic_estimates(:)
      real(dp), allocatable :: pic_variance(:)
      real(dp), allocatable :: rhs(:)
      real(dp), allocatable :: solved(:)
      real(dp), allocatable :: centered(:)
      real(dp) :: logdet
      real(dp) :: residual_sum
      real(dp) :: root_mean
      integer :: i
      integer :: m
      integer :: n
      integer :: status

      result = ace_continuous_result()
      call prepare_continuous_problem(tree, phenotype, work, laplacian, rhs, info, kappa)
      if (info /= 0) return
      n = work%n_tip
      m = work%n_node
      call ace_pic(work, phenotype, pic_estimates, pic_variance, scaled=.true., info=status)
      if (status /= 0) then
         info = 10 + status
         return
      end if
      root_mean = pic_estimates(work%root() - work%n_tip)
      call phylogenetic_vcv(work, covariance, status)
      if (status /= 0) then
         info = 20 + status
         return
      end if
      allocate(centered(n), solved(n))
      centered = phenotype - root_mean
      call solve_spd(covariance, centered, solved, status)
      if (status /= 0) then
         info = 30 + status
         return
      end if
      result%sigma2 = dot_product(centered, solved) / real(n, dp)
      if (result%sigma2 <= 0.0_dp) then
         info = 40
         return
      end if
      result%sigma2_standard_error = result%sigma2 * sqrt(2.0_dp / real(n, dp))
      allocate(result%estimates(m), result%standard_error(m))
      call solve_spd(laplacian, rhs, result%estimates, status)
      if (status /= 0) then
         info = 50 + status
         return
      end if
      residual_sum = edge_residual_sum(work, phenotype, result%estimates)
      result%residual_log_likelihood = -residual_sum / (2.0_dp * result%sigma2) &
         - real(m, dp) * log(result%sigma2)
      call spd_inverse_logdet(laplacian, inverse_laplacian, logdet, status)
      if (status /= 0) then
         info = 60 + status
         return
      end if
      do i = 1, m
         result%standard_error(i) = sqrt(result%sigma2 * inverse_laplacian(i, i))
      end do
   end subroutine ace_continuous_reml

   subroutine ace_continuous_gls(tree, phenotype, model, parameter, result, info)
      !! Reproduces ape's continuous `ace(..., method="GLS")` ancestral-state calculation.
      type(phylo_tree), intent(in) :: tree !! Rooted tree with branch lengths used by the chosen correlation structure.
      real(dp), intent(in) :: phenotype(:) !! Continuous trait values in numeric tip order.
      character(len=*), intent(in) :: model !! Supported ape GLS structure: `brownian`, `martins`, or `grafen`.
      real(dp), intent(in), optional :: parameter !! Martins alpha or Grafen rho; ignored for Brownian.
      type(ace_continuous_result), intent(out) :: result !! Ancestral GLS estimates, conditional SEs, and profiled variance.
      integer, intent(out) :: info !! Zero on success or nonzero for invalid model/data or linear algebra failure.
      type(phylo_tree) :: transformed
      type(pgls_result) :: gls
      real(dp), allocatable :: correlation(:, :)
      real(dp), allocatable :: covariance(:, :)
      real(dp), allocatable :: covariance_inverse(:, :)
      real(dp), allocatable :: full_matrix(:, :)
      real(dp), allocatable :: ones(:, :)
      real(dp), allocatable :: conditional(:, :)
      real(dp), allocatable :: cross(:, :)
      real(dp), allocatable :: centered(:)
      real(dp), allocatable :: solved(:)
      real(dp) :: alpha
      real(dp) :: logdet
      integer :: i
      integer :: m
      integer :: n
      integer :: status
      character(len=12) :: key

      result = ace_continuous_result()
      info = 0
      n = tree%n_tip
      m = tree%n_node
      if (.not. tree%valid() .or. .not. tree%has_lengths() .or. size(phenotype) /= n) then
         info = 1
         return
      end if
      if (.not. all(ieee_is_finite(phenotype)) .or. .not. all(ieee_is_finite(tree%edge_length))) then
         info = 2
         return
      end if
      key = lowercase_model(model)
      select case (trim(key))
      case ('brownian')
         call cor_brownian(tree, covariance, status, correlation=.false.)
         if (status == 0) call cor_brownian(tree, correlation, status, correlation=.true.)
         if (status == 0) call full_brownian_covariance(tree, full_matrix, status)
      case ('martins')
         if (.not. present(parameter)) then
            info = 3
            return
         end if
         alpha = parameter
         call cor_martins(tree, alpha, covariance, status)
         if (status == 0) correlation = covariance
         if (status == 0) then
            call dist_nodes(tree, full_matrix, status)
            if (status == 0) full_matrix = alpha * full_matrix
         end if
      case ('grafen')
         if (.not. present(parameter)) then
            info = 3
            return
         end if
         call grafen_tree(tree, parameter, transformed, status)
         if (status == 0) call cor_grafen(tree, parameter, covariance, status, correlation=.false.)
         if (status == 0) call cor_grafen(tree, parameter, correlation, status, correlation=.true.)
         if (status == 0) call full_brownian_covariance(tree, full_matrix, status)
      case default
         info = 4
         return
      end select
      if (status /= 0) then
         info = 10 + status
         return
      end if

      allocate(ones(n, 1))
      ones = 1.0_dp
      call pgls_fit(phenotype, ones, correlation, gls, status)
      if (status /= 0) then
         info = 20 + status
         return
      end if
      call spd_inverse_logdet(covariance, covariance_inverse, logdet, status)
      if (status /= 0) then
         info = 30 + status
         return
      end if
      allocate(cross(m, n), conditional(m, m), centered(n), solved(n))
      cross = full_matrix(n + 1:n + m, 1:n)
      conditional = full_matrix(n + 1:n + m, n + 1:n + m) &
         - matmul(matmul(cross, covariance_inverse), transpose(cross))
      centered = phenotype - gls%coefficients(1)
      solved = matmul(covariance_inverse, centered)
      allocate(result%estimates(m), result%standard_error(m))
      result%estimates = matmul(cross, solved) + gls%coefficients(1)
      do i = 1, m
         if (conditional(i, i) >= 0.0_dp) then
            result%standard_error(i) = sqrt(conditional(i, i))
         else
            result%standard_error(i) = quiet_nan()
         end if
      end do
      result%sigma2 = gls%sigma2
      result%log_likelihood = gls%log_likelihood
   end subroutine ace_continuous_gls

   subroutine full_brownian_covariance(tree, covariance, info)
      !! Builds Brownian root-to-MRCA covariance for all tips and internal nodes.
      type(phylo_tree), intent(in) :: tree !! Rooted tree with edge lengths.
      real(dp), allocatable, intent(out) :: covariance(:, :) !! All-node Brownian covariance in ape node-number order.
      integer, intent(out) :: info !! Zero on success or nonzero if node depths/MRCAs cannot be computed.
      real(dp), allocatable :: depth(:)
      integer :: ancestor
      integer :: i
      integer :: j
      integer :: n_total

      call node_depth_edgelength(tree, depth, info)
      n_total = tree%total_nodes()
      allocate(covariance(n_total, n_total))
      covariance = 0.0_dp
      if (info /= 0) return
      do i = 1, n_total
         do j = i, n_total
            ancestor = mrca(tree, i, j)
            if (ancestor <= 0) then
               info = 2
               return
            end if
            covariance(i, j) = depth(ancestor)
            covariance(j, i) = covariance(i, j)
         end do
      end do
   end subroutine full_brownian_covariance

   pure function lowercase_model(model) result(key)
      !! Normalizes supported ACE GLS model names.
      character(len=*), intent(in) :: model !! User model name matched case-insensitively.
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
      case default
         key = 'invalid'
      end select
   end function lowercase_model

   pure real(dp) function quiet_nan() result(value)
      !! Returns an IEEE quiet NaN for undefined ancestral-state standard errors.
      use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
      value = ieee_value(0.0_dp, ieee_quiet_nan)
   end function quiet_nan

   subroutine prepare_continuous_problem(tree, phenotype, work, laplacian, rhs, info, kappa)
      !! Validates a Brownian ACE problem and constructs its internal-node weighted Laplacian system.
      type(phylo_tree), intent(in) :: tree !! Candidate rooted binary tree with branch lengths.
      real(dp), intent(in) :: phenotype(:) !! Tip trait values in numeric tip order.
      type(phylo_tree), intent(out) :: work !! Validated tree copy with optional kappa-transformed branch lengths.
      real(dp), allocatable, intent(out) :: laplacian(:, :) !! Internal-node weighted graph Laplacian.
      real(dp), allocatable, intent(out) :: rhs(:) !! Weighted tip-data right-hand side for ancestral states.
      integer, intent(out) :: info !! Zero on success or a validation status code.
      real(dp), intent(in), optional :: kappa !! Optional nonnegative branch-length power.
      integer, allocatable :: nchild(:)
      real(dp) :: exponent
      real(dp) :: weight
      integer :: child
      integer :: child_index
      integer :: e
      integer :: m
      integer :: parent
      integer :: parent_index

      info = 0
      work = tree
      if (.not. tree%valid() .or. .not. tree%has_lengths()) then
         info = 1
         return
      end if
      if (size(phenotype) /= tree%n_tip .or. tree%n_node /= tree%n_tip - 1) then
         info = 2
         return
      end if
      if (.not. all(ieee_is_finite(phenotype)) .or. .not. all(ieee_is_finite(tree%edge_length))) then
         info = 6
         return
      end if
      nchild = child_counts(tree)
      do parent = tree%n_tip + 1, tree%total_nodes()
         if (nchild(parent) /= 2) then
            info = 3
            return
         end if
      end do
      exponent = 1.0_dp
      if (present(kappa)) exponent = kappa
      if (.not. ieee_is_finite(exponent) .or. exponent < 0.0_dp) then
         info = 4
         return
      end if
      work%edge_length = work%edge_length**exponent
      if (any(work%edge_length <= 0.0_dp)) then
         info = 5
         return
      end if

      m = work%n_node
      allocate(laplacian(m, m), rhs(m))
      laplacian = 0.0_dp
      rhs = 0.0_dp
      do e = 1, work%nedge()
         parent = work%edge(e, 1)
         child = work%edge(e, 2)
         parent_index = parent - work%n_tip
         weight = 1.0_dp / work%edge_length(e)
         if (child <= work%n_tip) then
            laplacian(parent_index, parent_index) = laplacian(parent_index, parent_index) + weight
            rhs(parent_index) = rhs(parent_index) + weight * phenotype(child)
         else
            child_index = child - work%n_tip
            laplacian(parent_index, parent_index) = laplacian(parent_index, parent_index) + weight
            laplacian(child_index, child_index) = laplacian(child_index, child_index) + weight
            laplacian(parent_index, child_index) = laplacian(parent_index, child_index) - weight
            laplacian(child_index, parent_index) = laplacian(child_index, parent_index) - weight
         end if
      end do
   end subroutine prepare_continuous_problem

   pure real(dp) function edge_residual_sum(tree, phenotype, estimates) result(value)
      !! Returns `sum((parent-child)**2 / edge.length)` for fitted ancestral states.
      type(phylo_tree), intent(in) :: tree !! Rooted tree used for the Brownian edge objective.
      real(dp), intent(in) :: phenotype(:) !! Observed tip values in numeric tip order.
      real(dp), intent(in) :: estimates(:) !! Internal-node state estimates in ape node-number order.
      real(dp) :: child_value
      real(dp) :: parent_value
      integer :: child
      integer :: e
      integer :: parent

      value = 0.0_dp
      do e = 1, tree%nedge()
         parent = tree%edge(e, 1)
         child = tree%edge(e, 2)
         parent_value = estimates(parent - tree%n_tip)
         if (child <= tree%n_tip) then
            child_value = phenotype(child)
         else
            child_value = estimates(child - tree%n_tip)
         end if
         value = value + (parent_value - child_value)**2 / tree%edge_length(e)
      end do
   end function edge_residual_sum

end module ape_ace
