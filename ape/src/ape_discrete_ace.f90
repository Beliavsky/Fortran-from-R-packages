! SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only
!
! Discrete ancestral-character likelihoods and rate fitting translated from
! ape R/ace.R (Copyright 2005-2024 Emmanuel Paradis and 2005 Ben Bolker).
module ape_discrete_ace
   use r_kinds, only : dp
   use r_linalg, only : general_complex_eigen, solve_system, inverse_matrix
   use ape_types, only : phylo_tree, child_counts
   use ape_optimize, only : bounded_problem, bounded_bfgs, finite_difference_hessian
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   implicit none
   private

   type, public :: ace_discrete_result
      real(dp), allocatable :: rates(:)
      real(dp), allocatable :: standard_error(:)
      real(dp), allocatable :: ancestral_likelihood(:, :)
      integer, allocatable :: index_matrix(:, :)
      real(dp) :: log_likelihood = -huge(1.0_dp)
      integer :: iterations = 0
      logical :: converged = .false.
   end type ace_discrete_result

   type, extends(bounded_problem) :: discrete_ace_problem
      type(phylo_tree) :: tree
      integer, allocatable :: states(:)
      integer, allocatable :: index_matrix(:, :)
   contains
      procedure :: value => discrete_problem_value
      procedure :: gradient => discrete_problem_gradient
   end type discrete_ace_problem

   interface ace_discrete_fit
      module procedure ace_discrete_fit_named
      module procedure ace_discrete_fit_index
   end interface ace_discrete_fit

   public :: ace_discrete_fit
   public :: ace_discrete_likelihood
   public :: ace_rate_index_matrix

contains

   subroutine ace_rate_index_matrix(n_states, model, index_matrix, info)
      !! Builds ape's ER, SYM, or ARD transition-rate index matrix.
      integer, intent(in) :: n_states !! Number of discrete character states; must be at least two.
      character(len=*), intent(in) :: model !! Named ape rate model: `ER`, `SYM`, or `ARD`, case-insensitive.
      integer, allocatable, intent(out) :: index_matrix(:, :) !! Off-diagonal parameter indices; diagonal entries are zero.
      integer, intent(out) :: info !! Zero on success or nonzero for an invalid state count or model name.
      character(len=:), allocatable :: key
      integer :: i
      integer :: j
      integer :: parameter

      info = 0
      if (n_states < 2) then
         allocate(index_matrix(0, 0))
         info = 1
         return
      end if
      allocate(index_matrix(n_states, n_states))
      index_matrix = 0
      key = uppercase(trim(model))
      select case (key)
      case ('ER')
         do j = 1, n_states
            do i = 1, n_states
               if (i /= j) index_matrix(i, j) = 1
            end do
         end do
      case ('ARD')
         parameter = 0
         do j = 1, n_states
            do i = 1, n_states
               if (i == j) cycle
               parameter = parameter + 1
               index_matrix(i, j) = parameter
            end do
         end do
      case ('SYM')
         parameter = 0
         do j = 1, n_states - 1
            do i = j + 1, n_states
               parameter = parameter + 1
               index_matrix(i, j) = parameter
               index_matrix(j, i) = parameter
            end do
         end do
      case default
         index_matrix = 0
         info = 2
      end select
   end subroutine ace_rate_index_matrix

   subroutine ace_discrete_fit_named(tree, states, n_states, model, result, info, initial_rate, marginal, max_iter, tolerance)
      !! Fits a named ape discrete-character continuous-time Markov model by maximum likelihood.
      type(phylo_tree), intent(in) :: tree !! Rooted fully dichotomous tree with nonnegative branch lengths.
      integer, intent(in) :: states(:) !! Tip states in numeric tip order; values 1..`n_states`, with zero denoting missing.
      integer, intent(in) :: n_states !! Number of observed/allowed character states.
      character(len=*), intent(in) :: model !! Named rate model: `ER`, `SYM`, or `ARD`.
      type(ace_discrete_result), intent(out) :: result !! Fitted rates, SEs, log likelihood, and internal-state likelihoods.
      integer, intent(out) :: info !! Zero on success or nonzero for validation, eigensystem, or optimization failure.
      real(dp), intent(in), optional :: initial_rate !! Common positive starting rate; default matches ape's `ip=0.1`.
      logical, intent(in), optional :: marginal !! Return upward marginal likelihoods when true; default ape joint smoothing.
      integer, intent(in), optional :: max_iter !! Maximum projected-BFGS iterations; default 500.
      real(dp), intent(in), optional :: tolerance !! Optimization tolerance; default `1e-8`.
      integer, allocatable :: index_matrix(:, :)
      integer :: status

      call ace_rate_index_matrix(n_states, model, index_matrix, status)
      if (status /= 0) then
         result = ace_discrete_result()
         info = status
         return
      end if
      call ace_discrete_fit_index(tree, states, index_matrix, result, info, initial_rate=initial_rate, &
         marginal=marginal, max_iter=max_iter, tolerance=tolerance)
   end subroutine ace_discrete_fit_named

   subroutine ace_discrete_fit_index(tree, states, index_matrix, result, info, initial_rate, initial_rates, marginal, &
      max_iter, tolerance)
      !! Fits a custom ape-style transition-rate index matrix by maximum likelihood.
      type(phylo_tree), intent(in) :: tree !! Rooted fully dichotomous tree with nonnegative branch lengths.
      integer, intent(in) :: states(:) !! Tip states in numeric order; zero denotes a missing/unknown state.
      integer, intent(in) :: index_matrix(:, :) !! Square off-diagonal map from transitions to 1-based rate parameters.
      type(ace_discrete_result), intent(out) :: result !! Fitted rates, SEs, log likelihood, and internal-state likelihoods.
      integer, intent(out) :: info !! Zero on success or nonzero for validation, eigensystem, or optimization failure.
      real(dp), intent(in), optional :: initial_rate !! Common starting rate used when `initial_rates` is absent; default 0.1.
      real(dp), intent(in), optional :: initial_rates(:) !! Optional separate starting value for every indexed rate parameter.
      logical, intent(in), optional :: marginal !! Return upward marginal likelihoods when true; default ape joint smoothing.
      integer, intent(in), optional :: max_iter !! Maximum projected-BFGS iterations; default 500.
      real(dp), intent(in), optional :: tolerance !! Optimization tolerance; default `1e-8`.
      type(discrete_ace_problem) :: problem
      real(dp), allocatable :: hessian(:, :)
      real(dp), allocatable :: hessian_inverse(:, :)
      real(dp), allocatable :: lower(:)
      real(dp), allocatable :: upper(:)
      real(dp), allocatable :: parameters(:)
      real(dp), allocatable :: ancestral(:, :)
      real(dp) :: deviance
      real(dp) :: start
      integer :: hinfo
      integer :: iterations
      integer :: np
      integer :: opt_info
      integer :: status
      logical :: use_marginal

      result = ace_discrete_result()
      info = 0
      call validate_discrete_input(tree, states, index_matrix, status)
      if (status /= 0) then
         info = status
         return
      end if
      np = maxval(index_matrix)
      if (np < 1) then
         info = 10
         return
      end if
      allocate(parameters(np), lower(np), upper(np))
      start = 0.1_dp
      if (present(initial_rate)) start = initial_rate
      if (.not. ieee_is_finite(start) .or. start < 0.0_dp) then
         info = 11
         return
      end if
      parameters = start
      if (present(initial_rates)) then
         if (size(initial_rates) /= np .or. any(initial_rates < 0.0_dp) .or. .not. all(ieee_is_finite(initial_rates))) then
            info = 12
            return
         end if
         parameters = initial_rates
      end if
      lower = 0.0_dp
      upper = 1.0e12_dp
      problem%tree = tree
      problem%states = states
      problem%index_matrix = index_matrix
      call bounded_bfgs(problem, parameters, lower, upper, deviance, opt_info, iterations, max_iter, tolerance)
      if (opt_info /= 0 .and. opt_info /= 7) then
         info = 20 + opt_info
         return
      end if

      use_marginal = .false.
      if (present(marginal)) use_marginal = marginal
      call ace_discrete_likelihood(tree, states, index_matrix, parameters, deviance, ancestral, status, use_marginal)
      if (status /= 0) then
         info = 40 + status
         return
      end if
      result%rates = parameters
      result%index_matrix = index_matrix
      result%ancestral_likelihood = ancestral
      result%log_likelihood = -0.5_dp * deviance
      result%iterations = iterations
      result%converged = opt_info == 0

      allocate(result%standard_error(np))
      result%standard_error = 0.0_dp
      call finite_difference_hessian(problem, parameters, hessian, hinfo, lower=lower, upper=upper)
      if (hinfo == 0 .and. all(ieee_is_finite(hessian))) then
         call inverse_matrix(hessian, hessian_inverse, hinfo)
         if (hinfo == 0) then
            do status = 1, np
               if (hessian_inverse(status, status) >= 0.0_dp) then
                  result%standard_error(status) = sqrt(hessian_inverse(status, status))
               else
                  result%standard_error(status) = ieee_nan()
               end if
            end do
         else
            result%standard_error = ieee_nan()
         end if
      else
         result%standard_error = ieee_nan()
      end if
   end subroutine ace_discrete_fit_index

   subroutine ace_discrete_likelihood(tree, states, index_matrix, rates, deviance, ancestral_likelihood, info, marginal)
      !! Evaluates ape's normalized pruning likelihood and optional internal-state smoothing at fixed rates.
      type(phylo_tree), intent(in) :: tree !! Rooted fully dichotomous tree with nonnegative branch lengths.
      integer, intent(in) :: states(:) !! Tip states in numeric order; zero denotes an unknown state.
      integer, intent(in) :: index_matrix(:, :) !! Square off-diagonal map from transitions to fitted rate indices.
      real(dp), intent(in) :: rates(:) !! Nonnegative transition-rate parameters referenced by `index_matrix`.
      real(dp), intent(out) :: deviance !! Minus twice the normalized pruning log likelihood, matching ape's `dev`.
      real(dp), allocatable, intent(out) :: ancestral_likelihood(:, :) !! Internal-node state likelihoods by node-number order.
      integer, intent(out) :: info !! Zero on success or nonzero for invalid input or failed eigensystem propagation.
      logical, intent(in), optional :: marginal !! Keep upward conditional likelihoods when true; otherwise smooth as ape does.
      real(dp), allocatable :: comp(:)
      real(dp), allocatable :: left(:)
      real(dp), allocatable :: likelihood(:, :)
      real(dp), allocatable :: right(:)
      real(dp), allocatable :: smoothed(:, :)
      complex(dp), allocatable :: eigenvalues(:)
      complex(dp), allocatable :: eigenvectors(:, :)
      real(dp), allocatable :: q(:, :)
      integer, allocatable :: child(:, :)
      integer, allocatable :: child_edge(:, :)
      logical, allocatable :: done(:)
      logical, allocatable :: smooth_done(:)
      integer :: anc
      integer :: child_node
      integer :: e
      integer :: i
      integer :: internal_count
      integer :: j
      integer :: n_states
      integer :: n_total
      integer :: root
      integer :: status
      logical :: progressed
      logical :: use_marginal

      status = 0
      deviance = huge(1.0_dp)
      call validate_discrete_input(tree, states, index_matrix, info)
      if (info /= 0) return
      if (size(rates) < maxval(index_matrix) .or. any(rates < 0.0_dp) .or. .not. all(ieee_is_finite(rates))) then
         info = 6
         return
      end if
      n_states = size(index_matrix, 1)
      n_total = tree%total_nodes()
      root = tree%root()
      allocate(q(n_states, n_states))
      call build_rate_matrix(index_matrix, rates, q)
      call general_complex_eigen(cmplx(q, 0.0_dp, kind=dp), eigenvalues, eigenvectors, status)
      if (status /= 0) then
         info = 10 + status
         return
      end if
      allocate(likelihood(n_total, n_states), comp(n_total), done(n_total))
      allocate(child(n_total, 2), child_edge(n_total, 2))
      likelihood = 0.0_dp
      comp = 0.0_dp
      done = .false.
      child = 0
      child_edge = 0
      do e = 1, tree%nedge()
         anc = tree%edge(e, 1)
         if (child(anc, 1) == 0) then
            i = 1
         else
            i = 2
         end if
         child(anc, i) = tree%edge(e, 2)
         child_edge(anc, i) = e
      end do
      do i = 1, tree%n_tip
         if (states(i) == 0) then
            likelihood(i, :) = 1.0_dp
         else
            likelihood(i, states(i)) = 1.0_dp
         end if
         done(i) = .true.
      end do
      internal_count = 0
      allocate(left(n_states), right(n_states))
      do while (internal_count < tree%n_node)
         progressed = .false.
         do anc = tree%n_tip + 1, n_total
            if (done(anc)) cycle
            if (child(anc, 1) == 0 .or. child(anc, 2) == 0) cycle
            if (.not. done(child(anc, 1)) .or. .not. done(child(anc, 2))) cycle
            call transition_vector(eigenvalues, eigenvectors, tree%edge_length(child_edge(anc, 1)), &
               likelihood(child(anc, 1), :), left, status)
            if (status /= 0) then
               info = 20 + status
               return
            end if
            call transition_vector(eigenvalues, eigenvectors, tree%edge_length(child_edge(anc, 2)), &
               likelihood(child(anc, 2), :), right, status)
            if (status /= 0) then
               info = 30 + status
               return
            end if
            likelihood(anc, :) = left * right
            comp(anc) = sum(likelihood(anc, :))
            if (.not. ieee_is_finite(comp(anc)) .or. comp(anc) <= tiny(1.0_dp)) then
               info = 40
               return
            end if
            likelihood(anc, :) = likelihood(anc, :) / comp(anc)
            done(anc) = .true.
            internal_count = internal_count + 1
            progressed = .true.
         end do
         if (.not. progressed) then
            info = 41
            return
         end if
      end do
      deviance = -2.0_dp * sum(log(comp(tree%n_tip + 1:n_total)))
      if (.not. ieee_is_finite(deviance)) then
         info = 42
         return
      end if

      allocate(ancestral_likelihood(tree%n_node, n_states))
      ancestral_likelihood = likelihood(tree%n_tip + 1:n_total, :)
      use_marginal = .false.
      if (present(marginal)) use_marginal = marginal
      if (use_marginal) then
         info = 0
         return
      end if

      allocate(smoothed(tree%n_node, n_states), smooth_done(tree%n_node))
      smoothed = ancestral_likelihood
      smooth_done = .false.
      smooth_done(root - tree%n_tip) = .true.
      internal_count = 1
      do while (internal_count < tree%n_node)
         progressed = .false.
         do anc = tree%n_tip + 1, n_total
            if (.not. smooth_done(anc - tree%n_tip)) cycle
            do j = 1, 2
               child_node = child(anc, j)
               if (child_node <= tree%n_tip) cycle
               if (smooth_done(child_node - tree%n_tip)) cycle
               call smooth_child(smoothed(anc - tree%n_tip, :), &
                  ancestral_likelihood(child_node - tree%n_tip, :), eigenvalues, eigenvectors, &
                  tree%edge_length(child_edge(anc, j)), smoothed(child_node - tree%n_tip, :), status)
               if (status /= 0) then
                  info = 50 + status
                  return
               end if
               smooth_done(child_node - tree%n_tip) = .true.
               internal_count = internal_count + 1
               progressed = .true.
            end do
         end do
         if (.not. progressed) then
            info = 51
            return
         end if
      end do
      ancestral_likelihood = smoothed
      info = 0
   end subroutine ace_discrete_likelihood

   function discrete_problem_value(self, x) result(value)
      class(discrete_ace_problem), intent(inout) :: self !! Discrete ACE data attached to the bounded optimizer.
      real(dp), intent(in) :: x(:) !! Candidate nonnegative transition rates.
      real(dp) :: value
      real(dp), allocatable :: ancestral(:, :)
      integer :: status

      call ace_discrete_likelihood(self%tree, self%states, self%index_matrix, x, value, ancestral, status, marginal=.true.)
      if (status /= 0 .or. .not. ieee_is_finite(value)) value = huge(1.0_dp) / 100.0_dp
   end function discrete_problem_value

   subroutine discrete_problem_gradient(self, x, gradient)
      class(discrete_ace_problem), intent(inout) :: self !! Discrete ACE data attached to the bounded optimizer.
      real(dp), intent(in) :: x(:) !! Candidate nonnegative transition rates.
      real(dp), intent(out) :: gradient(:) !! Finite-difference gradient of minus twice the log likelihood.
      real(dp), allocatable :: xp(:)
      real(dp), allocatable :: xm(:)
      real(dp) :: f0
      real(dp) :: fm
      real(dp) :: fp
      real(dp) :: h
      integer :: i

      allocate(xp(size(x)), xm(size(x)))
      f0 = self%value(x)
      do i = 1, size(x)
         h = epsilon(1.0_dp)**(1.0_dp / 3.0_dp) * max(1.0_dp, abs(x(i)))
         xp = x
         xp(i) = x(i) + h
         fp = self%value(xp)
         if (x(i) > h) then
            xm = x
            xm(i) = x(i) - h
            fm = self%value(xm)
            gradient(i) = (fp - fm) / (2.0_dp * h)
         else
            gradient(i) = (fp - f0) / h
         end if
      end do
   end subroutine discrete_problem_gradient

   subroutine validate_discrete_input(tree, states, index_matrix, info)
      !! Checks the rooted-binary tree, tip states, and custom transition-index matrix.
      type(phylo_tree), intent(in) :: tree !! Tree to validate for discrete pruning.
      integer, intent(in) :: states(:) !! Candidate tip states, where zero is missing.
      integer, intent(in) :: index_matrix(:, :) !! Candidate square transition-rate index matrix.
      integer, intent(out) :: info !! Zero when the discrete likelihood inputs are valid.
      integer, allocatable :: counts(:)
      integer :: i
      integer :: n_states

      info = 0
      if (.not. tree%valid() .or. .not. tree%has_lengths()) then
         info = 1
         return
      end if
      if (tree%n_node /= tree%n_tip - 1 .or. size(states) /= tree%n_tip) then
         info = 2
         return
      end if
      if (any(tree%edge_length < 0.0_dp) .or. .not. all(ieee_is_finite(tree%edge_length))) then
         info = 3
         return
      end if
      n_states = size(index_matrix, 1)
      if (n_states < 2 .or. size(index_matrix, 2) /= n_states) then
         info = 4
         return
      end if
      if (any(states < 0) .or. any(states > n_states) .or. any(index_matrix < 0)) then
         info = 5
         return
      end if
      do i = 1, n_states
         if (index_matrix(i, i) /= 0) then
            info = 5
            return
         end if
      end do
      counts = child_counts(tree)
      do i = tree%n_tip + 1, tree%total_nodes()
         if (counts(i) /= 2) then
            info = 2
            return
         end if
      end do
   end subroutine validate_discrete_input

   pure subroutine build_rate_matrix(index_matrix, rates, q)
      !! Expands an ape-style rate-index matrix into a continuous-time Markov generator.
      integer, intent(in) :: index_matrix(:, :) !! Off-diagonal map from transitions to 1-based rate parameters.
      real(dp), intent(in) :: rates(:) !! Nonnegative rate values referenced by the index matrix.
      real(dp), intent(out) :: q(:, :) !! Generator matrix with off-diagonal rates and rows summing to zero.
      integer :: i
      integer :: j
      integer :: parameter

      q = 0.0_dp
      do j = 1, size(q, 2)
         do i = 1, size(q, 1)
            if (i == j) cycle
            parameter = index_matrix(i, j)
            if (parameter > 0 .and. parameter <= size(rates)) q(i, j) = rates(parameter)
         end do
      end do
      do i = 1, size(q, 1)
         q(i, i) = -sum(q(i, :))
      end do
   end subroutine build_rate_matrix

   subroutine transition_vector(eigenvalues, eigenvectors, branch_length, input, output, info)
      !! Applies `exp(Q*t)` to one state-likelihood vector from a cached eigendecomposition of `Q`.
      complex(dp), intent(in) :: eigenvalues(:) !! Eigenvalues of the continuous-time Markov generator.
      complex(dp), intent(in) :: eigenvectors(:, :) !! Right eigenvectors of the generator as matrix columns.
      real(dp), intent(in) :: branch_length !! Nonnegative branch duration multiplying the generator.
      real(dp), intent(in) :: input(:) !! Child-state conditional likelihood vector.
      real(dp), intent(out) :: output(:) !! Parent-state conditional likelihood vector after transition propagation.
      integer, intent(out) :: info !! Zero on success or the shared linear-solve status code.
      complex(dp), allocatable :: coefficients(:)
      complex(dp), allocatable :: rhs(:)
      complex(dp), allocatable :: propagated(:)
      integer :: i

      info = 0
      allocate(coefficients(size(input)), rhs(size(input)), propagated(size(input)))
      rhs = cmplx(input, 0.0_dp, kind=dp)
      call solve_system(eigenvectors, rhs, coefficients, info)
      if (info /= 0) return
      propagated = matmul(eigenvectors, exp(eigenvalues * branch_length) * coefficients)
      output = real(propagated, dp)
      do i = 1, size(output)
         if (output(i) < 0.0_dp .and. output(i) > -1.0e-12_dp) output(i) = 0.0_dp
      end do
      if (any(output < -1.0e-10_dp) .or. .not. all(ieee_is_finite(output))) info = 1
   end subroutine transition_vector

   subroutine transition_matrix(eigenvalues, eigenvectors, branch_length, matrix, info)
      !! Forms `exp(Q*t)` from a cached generator eigendecomposition.
      complex(dp), intent(in) :: eigenvalues(:) !! Eigenvalues of the continuous-time Markov generator.
      complex(dp), intent(in) :: eigenvectors(:, :) !! Right eigenvectors of the generator as matrix columns.
      real(dp), intent(in) :: branch_length !! Nonnegative branch duration multiplying the generator.
      real(dp), intent(out) :: matrix(:, :) !! Transition matrix with the same order as the eigensystem.
      integer, intent(out) :: info !! Zero on success or nonzero when propagation fails.
      real(dp), allocatable :: basis(:)
      real(dp), allocatable :: column(:)
      integer :: j

      info = 0
      allocate(basis(size(matrix, 1)), column(size(matrix, 1)))
      matrix = 0.0_dp
      do j = 1, size(matrix, 2)
         basis = 0.0_dp
         basis(j) = 1.0_dp
         call transition_vector(eigenvalues, eigenvectors, branch_length, basis, column, info)
         if (info /= 0) return
         matrix(:, j) = column
      end do
   end subroutine transition_matrix

   subroutine smooth_child(parent_likelihood, upward_child, eigenvalues, eigenvectors, branch_length, child_likelihood, info)
      !! Performs ape's downward internal-node likelihood update for one internal child.
      real(dp), intent(in) :: parent_likelihood(:) !! Already-smoothed likelihood vector at the parent node.
      real(dp), intent(in) :: upward_child(:) !! Upward pruning likelihood vector at the child node.
      complex(dp), intent(in) :: eigenvalues(:) !! Eigenvalues of the fitted Markov generator.
      complex(dp), intent(in) :: eigenvectors(:, :) !! Right eigenvectors of the fitted Markov generator.
      real(dp), intent(in) :: branch_length !! Branch duration connecting parent and child.
      real(dp), intent(out) :: child_likelihood(:) !! Smoothed state likelihood vector at the internal child.
      integer, intent(out) :: info !! Zero on success or nonzero for a failed transition or zero denominator.
      real(dp), allocatable :: denominator(:)
      real(dp), allocatable :: matrix(:, :)
      real(dp), allocatable :: ratio(:)
      real(dp) :: total

      info = 0
      allocate(matrix(size(upward_child), size(upward_child)), denominator(size(upward_child)), ratio(size(upward_child)))
      call transition_matrix(eigenvalues, eigenvectors, branch_length, matrix, info)
      if (info /= 0) return
      denominator = matmul(upward_child, matrix)
      if (any(denominator <= tiny(1.0_dp))) then
         info = 2
         return
      end if
      ratio = parent_likelihood / denominator
      child_likelihood = matmul(ratio, matrix) * upward_child
      total = sum(child_likelihood)
      if (.not. ieee_is_finite(total) .or. total <= tiny(1.0_dp)) then
         info = 3
         return
      end if
      child_likelihood = child_likelihood / total
      info = 0
   end subroutine smooth_child

   pure function uppercase(text) result(value)
      !! Converts ASCII letters to uppercase for case-insensitive model-name matching.
      character(len=*), intent(in) :: text !! Input character string.
      character(len=len(text)) :: value
      integer :: code
      integer :: i

      value = text
      do i = 1, len(text)
         code = iachar(value(i:i))
         if (code >= iachar('a') .and. code <= iachar('z')) value(i:i) = achar(code - 32)
      end do
   end function uppercase

   pure real(dp) function ieee_nan() result(value)
      !! Returns an IEEE quiet NaN without self-comparison idioms.
      use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
      value = ieee_value(0.0_dp, ieee_quiet_nan)
   end function ieee_nan

end module ape_discrete_ace
