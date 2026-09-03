! SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only
!
! Ornstein-Uhlenbeck comparative model translated from ape R/compar.ou.R
! (Copyright 2005-2010 Emmanuel Paradis). The numerical likelihood is kept;
! formula/data-frame plumbing is omitted.
module ape_compar_ou
   use r_kinds, only : dp
   use ape_types, only : phylo_tree, parent_vector
   use ape_tree_algorithms, only : branching_times, dist_nodes
   use ape_pgls, only : pgls_result, pgls_fit
   use ape_optimize, only : bounded_problem, bounded_bfgs, finite_difference_hessian
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   implicit none
   private

   type, public :: compar_ou_result
      real(dp) :: alpha = 0.0_dp
      real(dp) :: sigma2 = 0.0_dp
      real(dp), allocatable :: theta(:)
      real(dp), allocatable :: standard_error(:)
      real(dp) :: deviance = huge(1.0_dp)
      real(dp) :: log_likelihood = -huge(1.0_dp)
      integer :: iterations = 0
      logical :: converged = .false.
   end type compar_ou_result

   type, extends(bounded_problem) :: ou_profile_problem
      type(phylo_tree) :: tree
      real(dp), allocatable :: phenotype(:)
      integer, allocatable :: shift_nodes(:)
   contains
      procedure :: value => ou_profile_value
      procedure :: gradient => ou_profile_gradient
   end type ou_profile_problem

   public :: compar_ou_fit
   public :: compar_ou_likelihood

contains

   subroutine compar_ou_fit(tree, phenotype, result, info, shift_nodes, alpha, lower_alpha, upper_alpha, &
      max_iter, tolerance)
      !! Fits ape's `compar.ou` model, profiling theta and sigma2 for a fixed or estimated alpha.
      type(phylo_tree), intent(in) :: tree !! Rooted ultrametric tree with positive branch lengths.
      real(dp), intent(in) :: phenotype(:) !! Continuous tip character in numeric tip order.
      type(compar_ou_result), intent(out) :: result !! Fitted alpha, sigma2, regime optima, likelihood, and conditional SEs.
      integer, intent(out) :: info !! Zero on success or nonzero for invalid data, optimization, or linear algebra failure.
      integer, intent(in), optional :: shift_nodes(:) !! Internal nodes where the OU optimum changes; root is ignored.
      real(dp), intent(in), optional :: alpha !! Fixed positive alpha; omitted to profile alpha deterministically.
      real(dp), intent(in), optional :: lower_alpha !! Lower alpha bound when estimated; default `1e-8`.
      real(dp), intent(in), optional :: upper_alpha !! Upper alpha bound when estimated; default 100.
      integer, intent(in), optional :: max_iter !! Maximum projected-BFGS iterations for alpha profiling.
      real(dp), intent(in), optional :: tolerance !! Optimizer tolerance; default `1e-8`.
      type(ou_profile_problem) :: problem
      type(pgls_result) :: gls
      real(dp), allocatable :: design(:, :)
      real(dp), allocatable :: covariance(:, :)
      real(dp), allocatable :: hessian(:, :)
      real(dp) :: lower(1)
      real(dp) :: objective
      real(dp) :: theta_alpha(1)
      real(dp) :: upper(1)
      integer, allocatable :: shifts(:)
      integer :: hinfo
      integer :: iterations
      integer :: opt_info
      integer :: status

      result = compar_ou_result()
      info = 0
      call validate_ou_input(tree, phenotype, shift_nodes, shifts, status)
      if (status /= 0) then
         info = status
         return
      end if
      if (present(alpha)) then
         if (.not. ieee_is_finite(alpha) .or. alpha <= 0.0_dp) then
            info = 4
            return
         end if
         theta_alpha(1) = alpha
         iterations = 0
         opt_info = 0
      else
         lower(1) = 1.0e-8_dp
         upper(1) = 100.0_dp
         if (present(lower_alpha)) lower(1) = lower_alpha
         if (present(upper_alpha)) upper(1) = upper_alpha
         if (lower(1) <= 0.0_dp .or. lower(1) >= upper(1) .or. .not. all(ieee_is_finite([lower(1), upper(1)]))) then
            info = 5
            return
         end if
         theta_alpha(1) = max(lower(1), min(upper(1), 0.1_dp))
         problem%tree = tree
         problem%phenotype = phenotype
         problem%shift_nodes = shifts
         call bounded_bfgs(problem, theta_alpha, lower, upper, objective, opt_info, iterations, max_iter, tolerance)
         if (opt_info /= 0 .and. opt_info /= 7) then
            info = 10 + opt_info
            return
         end if
      end if

      call ou_matrices(tree, shifts, theta_alpha(1), design, covariance, status)
      if (status /= 0) then
         info = 20 + status
         return
      end if
      call pgls_fit(phenotype, design, covariance, gls, status)
      if (status /= 0) then
         info = 30 + status
         return
      end if
      result%alpha = theta_alpha(1)
      result%sigma2 = gls%sigma2
      result%theta = gls%coefficients
      result%deviance = -2.0_dp * gls%log_likelihood
      result%log_likelihood = gls%log_likelihood
      result%iterations = iterations
      result%converged = opt_info == 0
      allocate(result%standard_error(2 + size(result%theta)))
      result%standard_error = 0.0_dp
      result%standard_error(2) = result%sigma2 * sqrt(2.0_dp / real(tree%n_tip, dp))
      result%standard_error(3:) = gls%standard_error
      if (present(alpha)) then
         result%standard_error(1) = 0.0_dp
      else
         problem%tree = tree
         problem%phenotype = phenotype
         problem%shift_nodes = shifts
         call finite_difference_hessian(problem, theta_alpha, hessian, hinfo, lower=lower, upper=upper)
         if (hinfo == 0 .and. hessian(1, 1) > 0.0_dp .and. ieee_is_finite(hessian(1, 1))) then
            result%standard_error(1) = sqrt(1.0_dp / hessian(1, 1))
         else
            result%standard_error(1) = quiet_nan()
         end if
      end if
   end subroutine compar_ou_fit

   function compar_ou_likelihood(tree, phenotype, alpha, sigma2, theta, info, shift_nodes) result(deviance)
      !! Evaluates ape's OU deviance at fixed alpha, sigma2, and regime optimum values.
      type(phylo_tree), intent(in) :: tree !! Rooted ultrametric tree with positive branch lengths.
      real(dp), intent(in) :: phenotype(:) !! Continuous tip character in numeric tip order.
      real(dp), intent(in) :: alpha !! Positive OU attraction parameter.
      real(dp), intent(in) :: sigma2 !! Positive residual/evolutionary variance multiplier.
      real(dp), intent(in) :: theta(:) !! One optimum per requested shift regime plus the ancestral/root regime.
      integer, intent(out) :: info !! Zero on success or nonzero for invalid dimensions or linear algebra failure.
      integer, intent(in), optional :: shift_nodes(:) !! Internal optimum-shift nodes; root is ignored.
      real(dp) :: deviance
      real(dp), allocatable :: covariance(:, :)
      real(dp), allocatable :: covariance_inverse(:, :)
      real(dp), allocatable :: design(:, :)
      real(dp), allocatable :: residual(:)
      real(dp) :: logdet
      integer, allocatable :: shifts(:)
      integer :: status

      deviance = huge(1.0_dp) / 100.0_dp
      call validate_ou_input(tree, phenotype, shift_nodes, shifts, status)
      if (status /= 0) then
         info = status
         return
      end if
      if (.not. ieee_is_finite(alpha) .or. alpha <= 0.0_dp .or. .not. ieee_is_finite(sigma2) .or. sigma2 <= 0.0_dp) then
         info = 4
         return
      end if
      call ou_matrices(tree, shifts, alpha, design, covariance, status)
      if (status /= 0 .or. size(theta) /= size(design, 2)) then
         info = 5 + status
         return
      end if
      call inverse_spd(covariance, covariance_inverse, logdet, status)
      if (status /= 0) then
         info = 20 + status
         return
      end if
      allocate(residual(tree%n_tip))
      residual = phenotype - matmul(design, theta)
      deviance = real(tree%n_tip, dp) * log(2.0_dp * acos(-1.0_dp) * sigma2) + logdet &
         + dot_product(residual, matmul(covariance_inverse, residual)) / sigma2
      info = 0
   end function compar_ou_likelihood

   function ou_profile_value(self, x) result(value)
      class(ou_profile_problem), intent(inout) :: self !! OU profile problem holding tree, traits, and shift regimes.
      real(dp), intent(in) :: x(:) !! One-element vector containing the candidate positive alpha.
      real(dp) :: value
      type(pgls_result) :: gls
      real(dp), allocatable :: covariance(:, :)
      real(dp), allocatable :: design(:, :)
      integer :: status

      value = huge(1.0_dp) / 100.0_dp
      if (size(x) /= 1 .or. x(1) <= 0.0_dp) return
      call ou_matrices(self%tree, self%shift_nodes, x(1), design, covariance, status)
      if (status /= 0) return
      call pgls_fit(self%phenotype, design, covariance, gls, status)
      if (status /= 0) return
      value = -2.0_dp * gls%log_likelihood
   end function ou_profile_value

   subroutine ou_profile_gradient(self, x, gradient)
      class(ou_profile_problem), intent(inout) :: self !! OU profile problem used for numerical alpha differentiation.
      real(dp), intent(in) :: x(:) !! One-element candidate alpha vector.
      real(dp), intent(out) :: gradient(:) !! Central finite-difference gradient of the profile deviance.
      real(dp) :: fm
      real(dp) :: fp
      real(dp) :: h
      real(dp) :: xm(1)
      real(dp) :: xp(1)

      h = epsilon(1.0_dp)**(1.0_dp / 3.0_dp) * max(1.0_dp, abs(x(1)))
      xp = x
      xm = x
      xp(1) = x(1) + h
      xm(1) = max(tiny(1.0_dp), x(1) - h)
      fp = self%value(xp)
      fm = self%value(xm)
      gradient(1) = (fp - fm) / (xp(1) - xm(1))
   end subroutine ou_profile_gradient

   subroutine ou_matrices(tree, shifts, alpha, design, covariance, info)
      !! Builds the regime design and OU covariance used by ape's `compar.ou` likelihood.
      type(phylo_tree), intent(in) :: tree !! Rooted ultrametric tree with positive branch lengths.
      integer, intent(in) :: shifts(:) !! Internal shift nodes excluding the root.
      real(dp), intent(in) :: alpha !! Positive OU attraction parameter.
      real(dp), allocatable, intent(out) :: design(:, :) !! Tip-by-regime mean-design matrix.
      real(dp), allocatable, intent(out) :: covariance(:, :) !! Tip OU covariance matrix before multiplying by sigma2.
      integer, intent(out) :: info !! Zero on success or nonzero for invalid tree ages/distances.
      real(dp), allocatable :: bt(:)
      real(dp), allocatable :: distance(:, :)
      real(dp), allocatable :: node_age(:)
      real(dp), allocatable :: wend(:, :)
      real(dp), allocatable :: wstart(:, :)
      integer, allocatable :: parent(:)
      integer, allocatable :: ordered_shift(:)
      integer, allocatable :: path(:)
      real(dp) :: tmax
      integer :: column
      integer :: i
      integer :: j
      integer :: root
      integer :: current_column

      info = 0
      root = tree%root()
      call branching_times(tree, bt, info)
      if (info /= 0) then
         allocate(design(0, 0), covariance(0, 0))
         return
      end if
      allocate(node_age(tree%total_nodes()))
      node_age = 0.0_dp
      node_age(tree%n_tip + 1:) = bt
      tmax = node_age(root)
      if (tmax <= 0.0_dp) then
         allocate(design(0, 0), covariance(0, 0))
         info = 2
         return
      end if
      ordered_shift = shifts
      call sort_nodes_by_age(ordered_shift, node_age)
      allocate(wend(tree%n_tip, size(shifts) + 1), wstart(tree%n_tip, size(shifts) + 1))
      wend = 0.0_dp
      wstart = 0.0_dp
      wstart(:, size(shifts) + 1) = tmax
      parent = parent_vector(tree)
      do i = 1, tree%n_tip
         call root_to_tip(parent, root, i, path)
         current_column = size(shifts) + 1
         do j = 2, size(path) - 1
            column = find_node(ordered_shift, path(j))
            if (column == 0) cycle
            wend(i, current_column) = node_age(path(j))
            wstart(i, column) = node_age(path(j))
            current_column = column
         end do
      end do
      allocate(design(tree%n_tip, size(shifts) + 1))
      design = exp(-alpha * wend) - exp(-alpha * wstart)
      call dist_nodes(tree, distance, info)
      if (info /= 0) then
         allocate(covariance(0, 0))
         return
      end if
      allocate(covariance(tree%n_tip, tree%n_tip))
      do j = 1, tree%n_tip
         do i = 1, tree%n_tip
            covariance(i, j) = exp(-alpha * distance(i, j)) &
               * (1.0_dp - exp(-2.0_dp * alpha * (tmax - 0.5_dp * distance(i, j))))
         end do
      end do
      if (.not. all(ieee_is_finite(covariance))) info = 3
   end subroutine ou_matrices

   subroutine validate_ou_input(tree, phenotype, shift_nodes, shifts, info)
      !! Validates compar.ou tree/data input and removes any root entry from the shift-node list.
      type(phylo_tree), intent(in) :: tree !! Candidate rooted tree.
      real(dp), intent(in) :: phenotype(:) !! Candidate numeric tip character.
      integer, intent(in), optional :: shift_nodes(:) !! Candidate internal shift nodes.
      integer, allocatable, intent(out) :: shifts(:) !! Validated shift nodes excluding the root.
      integer, intent(out) :: info !! Zero on success or nonzero for invalid dimensions/nodes/branch lengths.
      integer :: i
      integer :: k
      integer :: root

      info = 0
      if (.not. tree%valid() .or. .not. tree%has_lengths() .or. any(tree%edge_length <= 0.0_dp)) then
         allocate(shifts(0))
         info = 1
         return
      end if
      if (size(phenotype) /= tree%n_tip .or. .not. all(ieee_is_finite(phenotype))) then
         allocate(shifts(0))
         info = 2
         return
      end if
      root = tree%root()
      if (.not. present(shift_nodes)) then
         allocate(shifts(0))
         return
      end if
      if (any(shift_nodes <= tree%n_tip) .or. any(shift_nodes > tree%total_nodes())) then
         allocate(shifts(0))
         info = 3
         return
      end if
      allocate(shifts(count(shift_nodes /= root)))
      k = 0
      do i = 1, size(shift_nodes)
         if (shift_nodes(i) == root) cycle
         if (k > 0) then
            if (any(shifts(1:k) == shift_nodes(i))) cycle
         end if
         k = k + 1
         shifts(k) = shift_nodes(i)
      end do
      if (k < size(shifts)) shifts = shifts(1:k)
   end subroutine validate_ou_input

   subroutine inverse_spd(matrix, inverse, logdet, info)
      !! Computes an SPD inverse and log determinant through the shared linear-algebra package.
      use r_linalg, only : spd_inverse_logdet
      real(dp), intent(in) :: matrix(:, :) !! Symmetric positive-definite matrix.
      real(dp), allocatable, intent(out) :: inverse(:, :) !! Matrix inverse on success.
      real(dp), intent(out) :: logdet !! Natural logarithm of the determinant.
      integer, intent(out) :: info !! Shared SPD factorization status.

      call spd_inverse_logdet(matrix, inverse, logdet, info)
   end subroutine inverse_spd

   pure subroutine sort_nodes_by_age(nodes, age)
      !! Stable insertion-sorts node numbers by increasing branching age, matching `sort(bt[node])`.
      integer, intent(inout) :: nodes(:) !! Internal node numbers reordered in place.
      real(dp), intent(in) :: age(:) !! All-node age vector.
      integer :: i
      integer :: j
      integer :: key

      do i = 2, size(nodes)
         key = nodes(i)
         j = i - 1
         do while (j >= 1)
            if (age(nodes(j)) <= age(key)) exit
            nodes(j + 1) = nodes(j)
            j = j - 1
         end do
         nodes(j + 1) = key
      end do
   end subroutine sort_nodes_by_age

   pure integer function find_node(nodes, node) result(position)
      !! Finds a node's one-based regime-column position or zero when absent.
      integer, intent(in) :: nodes(:) !! Ordered shift-node list.
      integer, intent(in) :: node !! Ape node number to locate.
      integer :: i

      position = 0
      do i = 1, size(nodes)
         if (nodes(i) == node) then
            position = i
            return
         end if
      end do
   end function find_node

   pure subroutine root_to_tip(parent, root, tip, path)
      !! Builds one root-to-tip path from a node parent vector.
      integer, intent(in) :: parent(:) !! Parent node for every ape node number.
      integer, intent(in) :: root !! Root ape node number.
      integer, intent(in) :: tip !! Tip node number.
      integer, allocatable, intent(out) :: path(:) !! Root-to-tip node sequence.
      integer, allocatable :: reverse(:)
      integer :: current
      integer :: length

      allocate(reverse(size(parent)))
      current = tip
      length = 1
      reverse(1) = tip
      do while (current /= root)
         current = parent(current)
         if (current == 0) exit
         length = length + 1
         reverse(length) = current
      end do
      allocate(path(length))
      path = reverse(length:1:-1)
   end subroutine root_to_tip

   pure real(dp) function quiet_nan() result(value)
      !! Returns an IEEE quiet NaN for an undefined profile standard error.
      use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
      value = ieee_value(0.0_dp, ieee_quiet_nan)
   end function quiet_nan

end module ape_compar_ou
