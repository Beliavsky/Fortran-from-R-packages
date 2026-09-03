! SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only
!
! Penalized-likelihood molecular dating translated from ape R/chronopl.R
! (Copyright 2005-2012 Emmanuel Paradis). Cross-validation and R object
! attributes are intentionally omitted; the numerical fit is retained.
module ape_chronopl
   use r_kinds, only : dp
   use ape_types, only : phylo_tree, parent_vector
   use ape_optimize, only : bounded_problem, bounded_bfgs
   use ape_misc_statistics, only : regularized_gamma_q
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   implicit none
   private

   type, public :: chronopl_result
      type(phylo_tree) :: tree
      real(dp), allocatable :: rates(:)
      real(dp), allocatable :: node_age(:)
      real(dp) :: penalized_log_likelihood = -huge(1.0_dp)
      integer :: iterations = 0
      logical :: converged = .false.
   end type chronopl_result

   type, extends(bounded_problem) :: chronopl_problem
      type(phylo_tree) :: tree
      real(dp), allocatable :: scaled_length(:)
      real(dp), allocatable :: fixed_age(:)
      integer, allocatable :: unknown_node(:)
      integer, allocatable :: parent_edge(:)
      integer, allocatable :: basal_edge(:)
      real(dp) :: lambda = 1.0_dp
   contains
      procedure :: value => chronopl_problem_value
      procedure :: gradient => chronopl_problem_gradient
   end type chronopl_problem

   type, public :: chronos_clock_result
      type(phylo_tree) :: tree
      real(dp), allocatable :: node_age(:)
      real(dp) :: rate = 0.0_dp
      real(dp) :: log_likelihood = -huge(1.0_dp)
      real(dp) :: phiic = huge(1.0_dp)
      integer :: n_parameters = 0
      integer :: iterations = 0
      logical :: converged = .false.
   end type chronos_clock_result

   type, extends(bounded_problem) :: chronos_clock_problem
      type(phylo_tree) :: tree
      real(dp), allocatable :: observed_length(:)
      real(dp), allocatable :: fixed_age(:)
      integer, allocatable :: unknown_node(:)
   contains
      procedure :: value => chronos_clock_problem_value
      procedure :: gradient => chronos_clock_problem_gradient
   end type chronos_clock_problem

   type, public :: chronos_result
      type(phylo_tree) :: tree
      real(dp), allocatable :: rates(:)
      real(dp), allocatable :: frequencies(:)
      real(dp), allocatable :: node_age(:)
      real(dp) :: log_likelihood = -huge(1.0_dp)
      real(dp) :: penalized_log_likelihood = -huge(1.0_dp)
      real(dp) :: phiic = huge(1.0_dp)
      real(dp) :: lambda = 1.0_dp
      integer :: n_parameters = 0
      integer :: iterations = 0
      logical :: converged = .false.
      character(len=10) :: model = ""
   end type chronos_result

   type, extends(bounded_problem) :: chronos_problem
      type(phylo_tree) :: tree
      real(dp), allocatable :: observed_length(:)
      real(dp), allocatable :: fixed_age(:)
      integer, allocatable :: unknown_node(:)
      integer, allocatable :: parent_edge(:)
      integer, allocatable :: basal_edge(:)
      real(dp) :: lambda = 1.0_dp
      integer :: model_code = 1
      integer :: n_rate = 0
      integer :: n_frequency = 0
   contains
      procedure :: value => chronos_problem_value
      procedure :: gradient => chronos_problem_gradient
   end type chronos_problem

   public :: chronopl_fit
   public :: chronopl_objective
   public :: chronos_clock_fit
   public :: chronos_fit
   public :: chronos_objective

contains

   subroutine chronopl_fit(tree, lambda, result, info, age_min, age_max, calibration_nodes, scale, tolerance, max_iter)
      !! Fits ape's deterministic penalized-likelihood dating objective without cross-validation.
      type(phylo_tree), intent(in) :: tree !! Rooted tree with strictly positive observed branch lengths.
      real(dp), intent(in) :: lambda !! Nonnegative smoothing penalty multiplying adjacent-rate differences.
      type(chronopl_result), intent(out) :: result !! Dated tree, per-edge rates, node ages, objective, and convergence status.
      integer, intent(out) :: info !! Zero on success; nonzero for invalid data, calibration, or optimizer failure.
      real(dp), intent(in), optional :: age_min(:) !! Calibration lower/fixed ages; scalar or one value per calibration node.
      real(dp), intent(in), optional :: age_max(:) !! Optional calibration upper ages; scalar or one value per calibration node.
      integer, intent(in), optional :: calibration_nodes(:) !! Internal ape node numbers; default is the root only.
      real(dp), intent(in), optional :: scale !! Positive ape `S` branch-count scale; default one.
      real(dp), intent(in), optional :: tolerance !! Positive optimizer/bound tolerance; default `1e-8`.
      integer, intent(in), optional :: max_iter !! Maximum projected-BFGS iterations; default 500.
      type(chronopl_problem) :: problem
      real(dp), allocatable :: age(:)
      real(dp), allocatable :: calibration_lower(:)
      real(dp), allocatable :: calibration_upper(:)
      real(dp), allocatable :: lower(:)
      real(dp), allocatable :: parameters(:)
      real(dp), allocatable :: upper(:)
      integer, allocatable :: calibration(:)
      logical, allocatable :: fixed(:)
      real(dp) :: objective
      real(dp) :: s
      real(dp) :: tol
      integer :: e
      integer :: i
      integer :: iterations
      integer :: n_unknown
      integer :: opt_info
      integer :: status

      result = chronopl_result()
      info = 0
      if (.not. tree%valid() .or. .not. tree%has_lengths() .or. any(tree%edge_length <= 0.0_dp)) then
         info = 1
         return
      end if
      if (.not. all(ieee_is_finite(tree%edge_length)) .or. .not. ieee_is_finite(lambda) .or. lambda < 0.0_dp) then
         info = 2
         return
      end if
      s = 1.0_dp
      if (present(scale)) s = scale
      if (.not. ieee_is_finite(s) .or. s <= 0.0_dp) then
         info = 3
         return
      end if
      tol = 1.0e-8_dp
      if (present(tolerance)) tol = tolerance
      if (.not. ieee_is_finite(tol) .or. tol <= 0.0_dp .or. tol >= 0.5_dp) then
         info = 4
         return
      end if

      call prepare_calibrations(tree, age_min, age_max, calibration_nodes, calibration, calibration_lower, &
         calibration_upper, age, fixed, status)
      if (status /= 0) then
         info = 10 + status
         return
      end if
      call initialize_ages(tree, calibration, calibration_lower, calibration_upper, age, status)
      if (status /= 0) then
         info = 20 + status
         return
      end if
      call mark_fixed_calibrations(tree, calibration, calibration_lower, calibration_upper, age, fixed)
      call collect_unknown_nodes(tree, fixed, problem%unknown_node)
      n_unknown = size(problem%unknown_node)

      problem%tree = tree
      problem%lambda = lambda
      problem%scaled_length = tree%edge_length / s
      problem%fixed_age = age
      call build_edge_relations(tree, problem%parent_edge, problem%basal_edge)
      allocate(parameters(tree%nedge() + n_unknown))
      allocate(lower(size(parameters)), upper(size(parameters)))
      parameters(1:tree%nedge()) = tree%edge_length
      lower(1:tree%nedge()) = tol
      upper(1:tree%nedge()) = 1.0_dp - tol
      if (n_unknown > 0) then
         do i = 1, n_unknown
            parameters(tree%nedge() + i) = age(problem%unknown_node(i))
            lower(tree%nedge() + i) = tol
            upper(tree%nedge() + i) = 1.0_dp / tol
            call apply_calibration_bound(problem%unknown_node(i), calibration, calibration_lower, calibration_upper, &
               lower(tree%nedge() + i), upper(tree%nedge() + i))
         end do
      end if
      if (any(lower > upper)) then
         info = 30
         return
      end if
      parameters = max(lower, min(upper, parameters))
      call bounded_bfgs(problem, parameters, lower, upper, objective, opt_info, iterations, max_iter, tol)
      if (opt_info /= 0 .and. opt_info /= 7) then
         info = 40 + opt_info
         return
      end if

      result%tree = tree
      result%rates = parameters(1:tree%nedge())
      age = problem%fixed_age
      do i = 1, n_unknown
         age(problem%unknown_node(i)) = parameters(tree%nedge() + i)
      end do
      allocate(result%node_age(tree%n_node))
      result%node_age = age(tree%n_tip + 1:tree%total_nodes())
      do e = 1, tree%nedge()
         result%tree%edge_length(e) = age(tree%edge(e, 1)) - age(tree%edge(e, 2))
      end do
      if (any(result%tree%edge_length <= 0.0_dp)) then
         info = 50
         return
      end if
      result%penalized_log_likelihood = -objective
      result%iterations = iterations
      result%converged = opt_info == 0
   end subroutine chronopl_fit

   subroutine chronos_clock_fit(tree, result, info, age_min, age_max, calibration_nodes, tolerance, max_iter)
      !! Fits the single-rate `chronos(..., model="clock")` Poisson dating branch deterministically.
      type(phylo_tree), intent(in) :: tree !! Rooted tree with nonnegative observed branch lengths/counts.
      type(chronos_clock_result), intent(out) :: result !! Dated tree, common rate, ages, log likelihood, and PHIIC.
      integer, intent(out) :: info !! Zero on success or nonzero for invalid data, calibration, or optimization failure.
      real(dp), intent(in), optional :: age_min(:) !! Calibration lower/fixed ages; scalar or one per calibration node.
      real(dp), intent(in), optional :: age_max(:) !! Optional calibration upper ages; scalar or one per calibration node.
      integer, intent(in), optional :: calibration_nodes(:) !! Internal calibrated nodes; default root only.
      real(dp), intent(in), optional :: tolerance !! Positive optimizer/bound tolerance; default `1e-8`.
      integer, intent(in), optional :: max_iter !! Maximum projected-BFGS iterations; default 500.
      type(chronos_clock_problem) :: problem
      real(dp), allocatable :: age(:)
      real(dp), allocatable :: calibration_lower(:)
      real(dp), allocatable :: calibration_upper(:)
      real(dp), allocatable :: lower(:)
      real(dp), allocatable :: parameters(:)
      real(dp), allocatable :: real_length(:)
      real(dp), allocatable :: upper(:)
      integer, allocatable :: calibration(:)
      logical, allocatable :: fixed(:)
      real(dp) :: objective
      real(dp) :: tol
      integer :: e
      integer :: i
      integer :: iterations
      integer :: n_unknown
      integer :: opt_info
      integer :: status

      result = chronos_clock_result()
      info = 0
      if (.not. tree%valid() .or. .not. tree%has_lengths() .or. any(tree%edge_length < 0.0_dp)) then
         info = 1
         return
      end if
      if (.not. all(ieee_is_finite(tree%edge_length))) then
         info = 2
         return
      end if
      tol = 1.0e-8_dp
      if (present(tolerance)) tol = tolerance
      if (.not. ieee_is_finite(tol) .or. tol <= 0.0_dp .or. tol >= 0.5_dp) then
         info = 3
         return
      end if
      call prepare_calibrations(tree, age_min, age_max, calibration_nodes, calibration, calibration_lower, &
         calibration_upper, age, fixed, status)
      if (status /= 0) then
         info = 10 + status
         return
      end if
      call initialize_ages(tree, calibration, calibration_lower, calibration_upper, age, status)
      if (status /= 0) then
         info = 20 + status
         return
      end if
      call mark_fixed_calibrations(tree, calibration, calibration_lower, calibration_upper, age, fixed)
      call collect_unknown_nodes(tree, fixed, problem%unknown_node)
      n_unknown = size(problem%unknown_node)
      problem%tree = tree
      problem%observed_length = tree%edge_length
      problem%fixed_age = age
      allocate(real_length(tree%nedge()))
      do e = 1, tree%nedge()
         real_length(e) = age(tree%edge(e, 1)) - age(tree%edge(e, 2))
      end do
      if (any(real_length <= 0.0_dp)) then
         info = 30
         return
      end if
      allocate(parameters(1 + n_unknown), lower(1 + n_unknown), upper(1 + n_unknown))
      parameters(1) = 0.5_dp * (minval(tree%edge_length / real_length) + maxval(tree%edge_length / real_length))
      lower(1) = tol
      upper(1) = 1.0e5_dp - tol
      do i = 1, n_unknown
         parameters(1 + i) = age(problem%unknown_node(i))
         lower(1 + i) = tol
         upper(1 + i) = 1.0_dp / tol
         call apply_calibration_bound(problem%unknown_node(i), calibration, calibration_lower, calibration_upper, &
            lower(1 + i), upper(1 + i))
      end do
      parameters = max(lower, min(upper, parameters))
      call bounded_bfgs(problem, parameters, lower, upper, objective, opt_info, iterations, max_iter, tol)
      if (opt_info /= 0 .and. opt_info /= 7) then
         info = 40 + opt_info
         return
      end if
      result%tree = tree
      result%rate = parameters(1)
      age = problem%fixed_age
      do i = 1, n_unknown
         age(problem%unknown_node(i)) = parameters(1 + i)
      end do
      allocate(result%node_age(tree%n_node))
      result%node_age = age(tree%n_tip + 1:tree%total_nodes())
      do e = 1, tree%nedge()
         result%tree%edge_length(e) = age(tree%edge(e, 1)) - age(tree%edge(e, 2))
      end do
      result%log_likelihood = -objective
      result%n_parameters = 1 + n_unknown
      result%phiic = -2.0_dp * result%log_likelihood + 2.0_dp * real(result%n_parameters, dp)
      result%iterations = iterations
      result%converged = opt_info == 0
   end subroutine chronos_clock_fit

   recursive subroutine chronos_fit(tree, result, info, lambda, model, nb_rate_cat, age_min, age_max, calibration_nodes, &
      tolerance, max_iter)
      !! Fits ape chronos correlated, relaxed, discrete, or clock Poisson molecular-dating models.
      type(phylo_tree), intent(in) :: tree !! Rooted tree with nonnegative observed substitution/count branch lengths.
      type(chronos_result), intent(out) :: result !! Dated tree, fitted rates/frequencies, likelihoods, and PHIIC summary.
      integer, intent(out) :: info !! Zero on success; nonzero for invalid input, calibration, or optimization failure.
      real(dp), intent(in), optional :: lambda !! Nonnegative smoothing multiplier; default one.
      character(len=*), intent(in), optional :: model !! Model name: correlated, relaxed, discrete, or clock; default correlated.
      integer, intent(in), optional :: nb_rate_cat !! Number of discrete rate categories; default ten and capped at edge count.
      real(dp), intent(in), optional :: age_min(:) !! Calibration lower/fixed ages; scalar or one per calibration node.
      real(dp), intent(in), optional :: age_max(:) !! Optional calibration upper ages; scalar or one per calibration node.
      integer, intent(in), optional :: calibration_nodes(:) !! Internal calibrated nodes; default root only.
      real(dp), intent(in), optional :: tolerance !! Positive optimizer/bound tolerance; default `1e-8`.
      integer, intent(in), optional :: max_iter !! Maximum projected-BFGS iterations; default 1000.
      type(chronos_clock_result) :: clock_result
      type(chronos_problem) :: problem
      real(dp), allocatable :: age(:)
      real(dp), allocatable :: calibration_lower(:)
      real(dp), allocatable :: calibration_upper(:)
      real(dp), allocatable :: initial_edge_rate(:)
      real(dp), allocatable :: lower(:)
      real(dp), allocatable :: parameters(:)
      real(dp), allocatable :: real_length(:)
      real(dp), allocatable :: sorted_rate(:)
      real(dp), allocatable :: upper(:)
      integer, allocatable :: calibration(:)
      logical, allocatable :: fixed(:)
      character(len=10) :: model_name
      real(dp) :: lam
      real(dp) :: objective
      real(dp) :: tol
      integer :: e
      integer :: i
      integer :: iterations
      integer :: k
      integer :: n_edge
      integer :: n_frequency
      integer :: n_unknown
      integer :: opt_info
      integer :: status

      result = chronos_result()
      info = 0
      if (.not. tree%valid() .or. .not. tree%has_lengths() .or. any(tree%edge_length < 0.0_dp)) then
         info = 1
         return
      end if
      if (.not. all(ieee_is_finite(tree%edge_length))) then
         info = 2
         return
      end if
      lam = 1.0_dp
      if (present(lambda)) lam = lambda
      if (.not. ieee_is_finite(lam) .or. lam < 0.0_dp) then
         info = 3
         return
      end if
      model_name = "correlated"
      if (present(model)) model_name = lowercase_model(model)
      select case (trim(model_name))
      case ("clock")
         call chronos_clock_fit(tree, clock_result, status, age_min, age_max, calibration_nodes, tolerance, max_iter)
         if (status /= 0) then
            info = status
            return
         end if
         result%tree = clock_result%tree
         allocate(result%rates(1), result%frequencies(1))
         result%rates(1) = clock_result%rate
         result%frequencies(1) = 1.0_dp
         result%node_age = clock_result%node_age
         result%log_likelihood = clock_result%log_likelihood
         result%penalized_log_likelihood = clock_result%log_likelihood
         result%phiic = clock_result%phiic
         result%lambda = lam
         result%n_parameters = clock_result%n_parameters
         result%iterations = clock_result%iterations
         result%converged = clock_result%converged
         result%model = "clock"
         return
      case ("correlated")
         problem%model_code = 1
      case ("relaxed")
         problem%model_code = 2
      case ("discrete")
         problem%model_code = 3
      case default
         info = 4
         return
      end select
      tol = 1.0e-8_dp
      if (present(tolerance)) tol = tolerance
      if (.not. ieee_is_finite(tol) .or. tol <= 0.0_dp .or. tol >= 0.5_dp) then
         info = 5
         return
      end if

      call prepare_calibrations(tree, age_min, age_max, calibration_nodes, calibration, calibration_lower, &
         calibration_upper, age, fixed, status)
      if (status /= 0) then
         info = 10 + status
         return
      end if
      call initialize_ages(tree, calibration, calibration_lower, calibration_upper, age, status)
      if (status /= 0) then
         info = 20 + status
         return
      end if
      call mark_fixed_calibrations(tree, calibration, calibration_lower, calibration_upper, age, fixed)
      call collect_unknown_nodes(tree, fixed, problem%unknown_node)
      n_unknown = size(problem%unknown_node)
      n_edge = tree%nedge()
      allocate(real_length(n_edge), initial_edge_rate(n_edge))
      do e = 1, n_edge
         real_length(e) = age(tree%edge(e, 1)) - age(tree%edge(e, 2))
      end do
      if (any(real_length <= 0.0_dp)) then
         info = 30
         return
      end if
      initial_edge_rate = tree%edge_length / real_length

      problem%tree = tree
      problem%observed_length = tree%edge_length
      problem%fixed_age = age
      problem%lambda = lam
      call build_edge_relations(tree, problem%parent_edge, problem%basal_edge)
      if (problem%model_code == 3) then
         k = 10
         if (present(nb_rate_cat)) k = nb_rate_cat
         if (k < 1) then
            info = 31
            return
         end if
         k = min(k, n_edge)
         if (k == 1) then
            call chronos_fit(tree, result, info, lam, "clock", 1, age_min, age_max, calibration_nodes, tol, max_iter)
            if (info == 0) result%model = "discrete"
            return
         end if
         problem%n_rate = k
         problem%n_frequency = k - 1
      else
         problem%n_rate = n_edge
         problem%n_frequency = 0
      end if
      n_frequency = problem%n_frequency
      allocate(parameters(problem%n_rate + n_unknown + n_frequency))
      allocate(lower(size(parameters)), upper(size(parameters)))
      lower = tol
      upper = 1.0_dp / tol
      if (problem%model_code == 3) then
         sorted_rate = initial_edge_rate
         call sort_real_in_place(sorted_rate)
         do i = 1, problem%n_rate
            parameters(i) = quantile_type7_sorted(sorted_rate, (real(i, dp) - 0.5_dp) / real(problem%n_rate, dp))
         end do
      else
         parameters(1:n_edge) = initial_edge_rate
      end if
      lower(1:problem%n_rate) = tol
      upper(1:problem%n_rate) = 1.0e5_dp - tol
      do i = 1, n_unknown
         parameters(problem%n_rate + i) = age(problem%unknown_node(i))
         lower(problem%n_rate + i) = tol
         upper(problem%n_rate + i) = 1.0_dp / tol
         call apply_calibration_bound(problem%unknown_node(i), calibration, calibration_lower, calibration_upper, &
            lower(problem%n_rate + i), upper(problem%n_rate + i))
      end do
      if (n_frequency > 0) then
         parameters(problem%n_rate + n_unknown + 1:) = 1.0_dp / real(problem%n_rate, dp)
         lower(problem%n_rate + n_unknown + 1:) = 0.0_dp
         upper(problem%n_rate + n_unknown + 1:) = 1.0_dp
      end if
      parameters = max(lower, min(upper, parameters))
      call bounded_bfgs(problem, parameters, lower, upper, objective, opt_info, iterations, max_iter, tol)
      if (opt_info /= 0 .and. opt_info /= 7) then
         info = 40 + opt_info
         return
      end if

      result%tree = tree
      allocate(result%rates(problem%n_rate))
      result%rates = parameters(1:problem%n_rate)
      if (n_frequency > 0) then
         allocate(result%frequencies(problem%n_rate))
         result%frequencies(1:n_frequency) = parameters(problem%n_rate + n_unknown + 1:)
         result%frequencies(problem%n_rate) = 1.0_dp - sum(result%frequencies(1:n_frequency))
      else
         allocate(result%frequencies(0))
      end if
      age = problem%fixed_age
      do i = 1, n_unknown
         age(problem%unknown_node(i)) = parameters(problem%n_rate + i)
      end do
      allocate(result%node_age(tree%n_node))
      result%node_age = age(tree%n_tip + 1:tree%total_nodes())
      do e = 1, n_edge
         result%tree%edge_length(e) = age(tree%edge(e, 1)) - age(tree%edge(e, 2))
      end do
      result%log_likelihood = chronos_unpenalized_loglik(problem, parameters)
      result%penalized_log_likelihood = -objective
      result%n_parameters = size(parameters)
      result%lambda = lam
      result%phiic = chronos_phiic(problem, parameters, result%log_likelihood)
      result%iterations = iterations
      result%converged = opt_info == 0
      result%model = trim(model_name)
   end subroutine chronos_fit

   function chronos_objective(tree, observed_length, rates, ages, lambda, model, info, frequencies) result(value)
      !! Evaluates the negative chronos penalized likelihood at fixed rates, ages, and optional discrete frequencies.
      type(phylo_tree), intent(in) :: tree !! Rooted topology defining edge chronology and adjacency.
      real(dp), intent(in) :: observed_length(:) !! Nonnegative observed substitution/count length for each edge.
      real(dp), intent(in) :: rates(:) !! Positive edge rates or discrete rate-category values.
      real(dp), intent(in) :: ages(:) !! Age for every ape node number, with tip ages normally zero.
      real(dp), intent(in) :: lambda !! Nonnegative smoothing multiplier for correlated or relaxed models.
      character(len=*), intent(in) :: model !! One of correlated, relaxed, discrete, or clock.
      integer, intent(out) :: info !! Zero on success; nonzero for dimension, model, or numerical invalidity.
      real(dp), intent(in), optional :: frequencies(:) !! Full discrete-category probabilities summing to one.
      real(dp) :: value
      type(chronos_problem) :: problem
      real(dp), allocatable :: parameters(:)
      character(len=10) :: model_name
      integer :: e
      integer :: k

      info = 0
      value = huge(1.0_dp) / 100.0_dp
      if (.not. tree%valid() .or. size(observed_length) /= tree%nedge() .or. size(ages) /= tree%total_nodes()) then
         info = 1
         return
      end if
      if (any(observed_length < 0.0_dp) .or. any(rates <= 0.0_dp) .or. lambda < 0.0_dp) then
         info = 2
         return
      end if
      model_name = lowercase_model(model)
      problem%tree = tree
      problem%observed_length = observed_length
      problem%fixed_age = ages
      allocate(problem%unknown_node(0))
      problem%lambda = lambda
      call build_edge_relations(tree, problem%parent_edge, problem%basal_edge)
      select case (trim(model_name))
      case ("correlated")
         if (size(rates) /= tree%nedge()) then
            info = 3
            return
         end if
         problem%model_code = 1
         problem%n_rate = size(rates)
         problem%n_frequency = 0
         parameters = rates
      case ("relaxed")
         if (size(rates) /= tree%nedge()) then
            info = 3
            return
         end if
         problem%model_code = 2
         problem%n_rate = size(rates)
         problem%n_frequency = 0
         parameters = rates
      case ("clock")
         if (size(rates) /= 1) then
            info = 3
            return
         end if
         problem%model_code = 3
         problem%n_rate = 1
         problem%n_frequency = 0
         parameters = rates
      case ("discrete")
         if (.not. present(frequencies)) then
            info = 4
            return
         end if
         k = size(rates)
         if (size(frequencies) /= k .or. k < 1 .or. any(frequencies < 0.0_dp) &
            .or. abs(sum(frequencies) - 1.0_dp) > 1.0e-10_dp) then
            info = 5
            return
         end if
         problem%model_code = 3
         problem%n_rate = k
         problem%n_frequency = max(0, k - 1)
         allocate(parameters(k + problem%n_frequency))
         parameters(1:k) = rates
         if (k > 1) parameters(k + 1:) = frequencies(1:k - 1)
      case default
         info = 6
         return
      end select
      do e = 1, tree%nedge()
         if (ages(tree%edge(e, 1)) - ages(tree%edge(e, 2)) <= 0.0_dp) then
            info = 7
            return
         end if
      end do
      value = chronos_problem_value(problem, parameters)
      if (.not. ieee_is_finite(value)) info = 8
   end function chronos_objective

   function chronos_problem_value(self, x) result(value)
      class(chronos_problem), intent(inout) :: self !! Chronos data and model definition attached to the bounded optimizer.
      real(dp), intent(in) :: x(:) !! Rates followed by unknown ages and, for multi-rate discrete models, free frequencies.
      real(dp) :: value
      real(dp), allocatable :: age(:)
      real(dp), allocatable :: branch_time(:)
      real(dp), allocatable :: sorted_rate(:)
      real(dp) :: loglik
      real(dp) :: mean_basal
      real(dp) :: mean_rate
      real(dp) :: penalty
      real(dp) :: probability
      real(dp) :: target
      integer :: e
      integer :: i
      integer :: n_edge
      integer :: n_unknown

      value = huge(1.0_dp) / 100.0_dp
      n_edge = self%tree%nedge()
      n_unknown = size(self%unknown_node)
      if (size(x) /= self%n_rate + n_unknown + self%n_frequency .or. .not. all(ieee_is_finite(x))) return
      if (any(x(1:self%n_rate) <= 0.0_dp)) return
      age = self%fixed_age
      do i = 1, n_unknown
         age(self%unknown_node(i)) = x(self%n_rate + i)
      end do
      allocate(branch_time(n_edge))
      do e = 1, n_edge
         branch_time(e) = age(self%tree%edge(e, 1)) - age(self%tree%edge(e, 2))
      end do
      if (any(branch_time <= tiny(1.0_dp))) return

      select case (self%model_code)
      case (1)
         loglik = edge_poisson_loglik(self%observed_length, branch_time, x(1:self%n_rate))
         if (.not. ieee_is_finite(loglik)) return
         penalty = 0.0_dp
         do e = 1, n_edge
            if (self%parent_edge(e) > 0) penalty = penalty + (x(e) - x(self%parent_edge(e)))**2
         end do
         value = -loglik + self%lambda * penalty
         if (size(self%basal_edge) > 1) then
            mean_basal = sum(x(self%basal_edge)) / real(size(self%basal_edge), dp)
            value = value - self%lambda * sum((x(self%basal_edge) - mean_basal)**2) / &
               real(size(self%basal_edge) - 1, dp)
         end if
      case (2)
         loglik = edge_poisson_loglik(self%observed_length, branch_time, x(1:self%n_rate))
         if (.not. ieee_is_finite(loglik)) return
         sorted_rate = x(1:self%n_rate)
         call sort_real_in_place(sorted_rate)
         mean_rate = sum(sorted_rate) / real(self%n_rate, dp)
         penalty = 0.0_dp
         do i = 1, self%n_rate
            target = real(i, dp) / real(self%n_rate, dp)
            probability = 1.0_dp - regularized_gamma_q(mean_rate, sorted_rate(i))
            penalty = penalty + (target - probability)**2
         end do
         value = -loglik + self%lambda * penalty
      case (3)
         value = -discrete_poisson_loglik(self%observed_length, branch_time, x(1:self%n_rate), &
            x(self%n_rate + n_unknown + 1:))
      end select
   end function chronos_problem_value

   subroutine chronos_problem_gradient(self, x, gradient)
      class(chronos_problem), intent(inout) :: self !! Chronos model whose exact objective is differentiated numerically.
      real(dp), intent(in) :: x(:) !! Rates, unknown ages, and optional discrete frequencies at the evaluation point.
      real(dp), intent(out) :: gradient(:) !! Central/one-sided finite-difference gradient of the exact chronos objective.
      real(dp), allocatable :: xm(:)
      real(dp), allocatable :: xp(:)
      real(dp) :: f0
      real(dp) :: fm
      real(dp) :: fp
      real(dp) :: h
      real(dp), parameter :: bad = huge(1.0_dp) / 1000.0_dp
      integer :: i

      gradient = 0.0_dp
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
         if (abs(fm) < bad .and. abs(fp) < bad) then
            gradient(i) = (fp - fm) / (2.0_dp * h)
         else if (abs(fp) < bad .and. abs(f0) < bad) then
            gradient(i) = (fp - f0) / h
         else if (abs(fm) < bad .and. abs(f0) < bad) then
            gradient(i) = (f0 - fm) / h
         else
            gradient(i) = 0.0_dp
         end if
      end do
   end subroutine chronos_problem_gradient

   function chronos_unpenalized_loglik(problem, x) result(loglik)
      !! Returns the unpenalized Poisson or discrete-mixture log likelihood for a fitted chronos parameter vector.
      type(chronos_problem), intent(in) :: problem !! Fitted chronos model definition and fixed data.
      real(dp), intent(in) :: x(:) !! Fitted rates, unknown ages, and optional free frequencies.
      real(dp) :: loglik
      real(dp), allocatable :: age(:)
      real(dp), allocatable :: branch_time(:)
      integer :: e
      integer :: i
      integer :: n_unknown

      n_unknown = size(problem%unknown_node)
      age = problem%fixed_age
      do i = 1, n_unknown
         age(problem%unknown_node(i)) = x(problem%n_rate + i)
      end do
      allocate(branch_time(problem%tree%nedge()))
      do e = 1, problem%tree%nedge()
         branch_time(e) = age(problem%tree%edge(e, 1)) - age(problem%tree%edge(e, 2))
      end do
      if (problem%model_code == 3) then
         loglik = discrete_poisson_loglik(problem%observed_length, branch_time, x(1:problem%n_rate), &
            x(problem%n_rate + n_unknown + 1:))
      else
         loglik = edge_poisson_loglik(problem%observed_length, branch_time, x(1:problem%n_rate))
      end if
   end function chronos_unpenalized_loglik

   function chronos_phiic(problem, x, loglik) result(phiic)
      !! Computes ape's PHIIC expression for correlated, relaxed, and discrete chronos models.
      type(chronos_problem), intent(in) :: problem !! Fitted model topology, penalty definition, and parameter dimensions.
      real(dp), intent(in) :: x(:) !! Fitted rates, ages, and optional free frequencies.
      real(dp), intent(in) :: loglik !! Unpenalized fitted Poisson log likelihood.
      real(dp) :: phiic
      real(dp), allocatable :: phi(:)
      real(dp), allocatable :: sorted_rate(:)
      real(dp) :: mean_basal
      real(dp) :: mean_rate
      real(dp) :: probability
      real(dp) :: basal_variance
      integer :: e
      integer :: i
      integer :: k

      k = size(x)
      if (problem%model_code == 3) then
         phiic = -2.0_dp * loglik + 2.0_dp * real(k, dp)
         return
      end if
      if (problem%model_code == 1) then
         allocate(phi(count(problem%parent_edge > 0)))
         i = 0
         basal_variance = 0.0_dp
         if (size(problem%basal_edge) > 1) then
            mean_basal = sum(x(problem%basal_edge)) / real(size(problem%basal_edge), dp)
            basal_variance = sum((x(problem%basal_edge) - mean_basal)**2) / real(size(problem%basal_edge) - 1, dp)
         end if
         do e = 1, problem%tree%nedge()
            if (problem%parent_edge(e) <= 0) cycle
            i = i + 1
            phi(i) = (x(e) - x(problem%parent_edge(e)))**2 + basal_variance
         end do
      else
         allocate(phi(problem%n_rate))
         sorted_rate = x(1:problem%n_rate)
         call sort_real_in_place(sorted_rate)
         mean_rate = sum(sorted_rate) / real(problem%n_rate, dp)
         do i = 1, problem%n_rate
            probability = 1.0_dp - regularized_gamma_q(mean_rate, sorted_rate(i))
            phi(i) = (real(i, dp) / real(problem%n_rate, dp) - probability)**2
         end do
      end if
      phiic = -2.0_dp * loglik + 2.0_dp * real(k, dp) + problem%lambda * sqrt(sum(phi**2))
   end function chronos_phiic

   pure function edge_poisson_loglik(observed_length, branch_time, rates) result(loglik)
      !! Returns the independent-edge Poisson log likelihood for one fitted rate per edge.
      real(dp), intent(in) :: observed_length(:) !! Observed substitution/count branch lengths.
      real(dp), intent(in) :: branch_time(:) !! Positive chronological branch durations.
      real(dp), intent(in) :: rates(:) !! Positive per-edge rates.
      real(dp) :: loglik
      real(dp) :: mean_count
      integer :: e

      loglik = 0.0_dp
      if (size(observed_length) /= size(branch_time) .or. size(rates) /= size(branch_time)) then
         loglik = -huge(1.0_dp) / 100.0_dp
         return
      end if
      do e = 1, size(rates)
         mean_count = rates(e) * branch_time(e)
         if (mean_count <= 0.0_dp) then
            loglik = -huge(1.0_dp) / 100.0_dp
            return
         end if
         loglik = loglik + observed_length(e) * log(mean_count) - mean_count - log_gamma(observed_length(e) + 1.0_dp)
      end do
   end function edge_poisson_loglik

   function discrete_poisson_loglik(observed_length, branch_time, rates, free_frequency) result(loglik)
      !! Returns ape's branchwise finite-mixture Poisson log likelihood for discrete chronos rate categories.
      real(dp), intent(in) :: observed_length(:) !! Observed substitution/count branch lengths.
      real(dp), intent(in) :: branch_time(:) !! Positive chronological branch durations.
      real(dp), intent(in) :: rates(:) !! Positive discrete rate-category values.
      real(dp), intent(in) :: free_frequency(:) !! First `k-1` probabilities; final probability is one minus their sum.
      real(dp) :: loglik
      real(dp), allocatable :: frequency(:)
      real(dp), allocatable :: terms(:)
      real(dp) :: maximum_term
      real(dp) :: mean_count
      integer :: e
      integer :: j

      loglik = -huge(1.0_dp) / 100.0_dp
      if (size(free_frequency) /= max(0, size(rates) - 1)) return
      allocate(frequency(size(rates)), terms(size(rates)))
      if (size(rates) > 1) frequency(1:size(rates) - 1) = free_frequency
      frequency(size(rates)) = 1.0_dp - sum(free_frequency)
      if (any(frequency < 0.0_dp) .or. any(frequency > 1.0_dp)) return
      loglik = 0.0_dp
      do e = 1, size(branch_time)
         terms = -huge(1.0_dp) / 100.0_dp
         do j = 1, size(rates)
            if (frequency(j) <= 0.0_dp) cycle
            mean_count = rates(j) * branch_time(e)
            if (mean_count <= 0.0_dp) return
            terms(j) = log(frequency(j)) + observed_length(e) * log(mean_count) - mean_count - &
               log_gamma(observed_length(e) + 1.0_dp)
         end do
         maximum_term = maxval(terms)
         if (maximum_term < -huge(1.0_dp) / 1000.0_dp) return
         loglik = loglik + maximum_term + log(sum(exp(terms - maximum_term)))
      end do
   end function discrete_poisson_loglik

   pure function lowercase_model(text) result(output)
      !! Converts an ASCII model name to lower case for ape-style case-insensitive dispatch.
      character(len=*), intent(in) :: text !! Model name to normalize.
      character(len=10) :: output
      integer :: code
      integer :: i
      integer :: n

      output = ""
      n = min(len_trim(text), len(output))
      do i = 1, n
         code = iachar(text(i:i))
         if (code >= iachar('A') .and. code <= iachar('Z')) then
            output(i:i) = achar(code + iachar('a') - iachar('A'))
         else
            output(i:i) = text(i:i)
         end if
      end do
   end function lowercase_model

   pure subroutine sort_real_in_place(values)
      !! Stable insertion-sorts a real vector in ascending order.
      real(dp), intent(inout) :: values(:) !! Values sorted ascending in place.
      real(dp) :: key
      integer :: i
      integer :: j

      do i = 2, size(values)
         key = values(i)
         j = i - 1
         do while (j >= 1)
            if (values(j) <= key) exit
            values(j + 1) = values(j)
            j = j - 1
         end do
         values(j + 1) = key
      end do
   end subroutine sort_real_in_place

   pure real(dp) function quantile_type7_sorted(values, probability) result(value)
      !! Evaluates R's default type-7 sample quantile on an ascending sorted vector.
      real(dp), intent(in) :: values(:) !! Ascending sample values.
      real(dp), intent(in) :: probability !! Probability in [0,1].
      real(dp) :: fraction
      real(dp) :: position
      integer :: lower_index

      if (size(values) == 1) then
         value = values(1)
         return
      end if
      position = 1.0_dp + real(size(values) - 1, dp) * max(0.0_dp, min(1.0_dp, probability))
      lower_index = min(size(values) - 1, max(1, int(floor(position))))
      fraction = position - real(lower_index, dp)
      value = (1.0_dp - fraction) * values(lower_index) + fraction * values(lower_index + 1)
   end function quantile_type7_sorted

   function chronopl_objective(tree, observed_length, rates, ages, lambda, info) result(value)
      !! Evaluates the negative ape penalized log likelihood at supplied edge rates and all-node ages.
      type(phylo_tree), intent(in) :: tree !! Tree topology defining edge adjacency.
      real(dp), intent(in) :: observed_length(:) !! Scaled observed substitutions/counts for each tree edge.
      real(dp), intent(in) :: rates(:) !! Positive rate associated with each edge.
      real(dp), intent(in) :: ages(:) !! Age for every ape node number, including zero-age tips when appropriate.
      real(dp), intent(in) :: lambda !! Nonnegative rate-smoothing penalty.
      integer, intent(out) :: info !! Zero on success; nonzero for invalid dimensions or parameters.
      real(dp) :: value
      integer, allocatable :: basal(:)
      integer, allocatable :: parent_edge(:)
      real(dp), allocatable :: real_length(:)
      integer :: e

      info = 0
      value = huge(1.0_dp) / 100.0_dp
      if (size(observed_length) /= tree%nedge() .or. size(rates) /= tree%nedge() &
         .or. size(ages) /= tree%total_nodes()) then
         info = 1
         return
      end if
      if (any(observed_length < 0.0_dp) .or. any(rates <= 0.0_dp) .or. lambda < 0.0_dp) then
         info = 2
         return
      end if
      if (.not. all(ieee_is_finite(observed_length)) .or. .not. all(ieee_is_finite(rates)) &
         .or. .not. all(ieee_is_finite(ages)) .or. .not. ieee_is_finite(lambda)) then
         info = 3
         return
      end if
      allocate(real_length(tree%nedge()))
      do e = 1, tree%nedge()
         real_length(e) = ages(tree%edge(e, 1)) - ages(tree%edge(e, 2))
      end do
      if (any(real_length <= 0.0_dp)) then
         info = 4
         return
      end if
      call build_edge_relations(tree, parent_edge, basal)
      value = negative_penalized_loglik(observed_length, real_length, rates, parent_edge, basal, lambda)
   end function chronopl_objective

   function chronos_clock_problem_value(self, x) result(value)
      class(chronos_clock_problem), intent(inout) :: self !! Clock-model tree/data attached to the bounded optimizer.
      real(dp), intent(in) :: x(:) !! Common rate followed by unknown internal-node ages.
      real(dp) :: value
      real(dp), allocatable :: age(:)
      real(dp) :: branch_time
      real(dp) :: b
      integer :: e
      integer :: i

      value = huge(1.0_dp) / 100.0_dp
      if (size(x) /= 1 + size(self%unknown_node) .or. x(1) <= 0.0_dp .or. .not. all(ieee_is_finite(x))) return
      age = self%fixed_age
      do i = 1, size(self%unknown_node)
         age(self%unknown_node(i)) = x(1 + i)
      end do
      value = 0.0_dp
      do e = 1, self%tree%nedge()
         branch_time = age(self%tree%edge(e, 1)) - age(self%tree%edge(e, 2))
         if (branch_time <= tiny(1.0_dp)) then
            value = huge(1.0_dp) / 100.0_dp
            return
         end if
         b = x(1) * branch_time
         value = value - self%observed_length(e) * log(b) + b + log_gamma(self%observed_length(e) + 1.0_dp)
      end do
   end function chronos_clock_problem_value

   subroutine chronos_clock_problem_gradient(self, x, gradient)
      class(chronos_clock_problem), intent(inout) :: self !! Clock-model tree/data attached to the bounded optimizer.
      real(dp), intent(in) :: x(:) !! Common rate followed by unknown internal-node ages.
      real(dp), intent(out) :: gradient(:) !! Analytic gradient of the negative Poisson log likelihood.
      real(dp), allocatable :: age(:)
      real(dp), allocatable :: real_length(:)
      integer :: e
      integer :: i
      integer :: incoming
      integer :: node

      gradient = 0.0_dp
      age = self%fixed_age
      do i = 1, size(self%unknown_node)
         age(self%unknown_node(i)) = x(1 + i)
      end do
      allocate(real_length(self%tree%nedge()))
      do e = 1, self%tree%nedge()
         real_length(e) = age(self%tree%edge(e, 1)) - age(self%tree%edge(e, 2))
      end do
      if (any(real_length <= tiny(1.0_dp))) return
      gradient(1) = sum(real_length - self%observed_length / x(1))
      do i = 1, size(self%unknown_node)
         node = self%unknown_node(i)
         do e = 1, self%tree%nedge()
            if (self%tree%edge(e, 1) == node) then
               gradient(1 + i) = gradient(1 + i) + x(1) - self%observed_length(e) / real_length(e)
            end if
         end do
         incoming = 0
         do e = 1, self%tree%nedge()
            if (self%tree%edge(e, 2) == node) then
               incoming = e
               exit
            end if
         end do
         if (incoming > 0) then
            gradient(1 + i) = gradient(1 + i) - x(1) + self%observed_length(incoming) / real_length(incoming)
         end if
      end do
   end subroutine chronos_clock_problem_gradient

   function chronopl_problem_value(self, x) result(value)
      class(chronopl_problem), intent(inout) :: self !! Penalized-likelihood dating data attached to the optimizer.
      real(dp), intent(in) :: x(:) !! Edge rates followed by ages of currently unknown internal nodes.
      real(dp) :: value
      real(dp), allocatable :: age(:)
      real(dp), allocatable :: real_length(:)
      integer :: e
      integer :: i
      integer :: n_edge

      value = huge(1.0_dp) / 100.0_dp
      n_edge = self%tree%nedge()
      if (size(x) /= n_edge + size(self%unknown_node)) return
      if (any(x(1:n_edge) <= 0.0_dp) .or. .not. all(ieee_is_finite(x))) return
      age = self%fixed_age
      do i = 1, size(self%unknown_node)
         age(self%unknown_node(i)) = x(n_edge + i)
      end do
      allocate(real_length(n_edge))
      do e = 1, n_edge
         real_length(e) = age(self%tree%edge(e, 1)) - age(self%tree%edge(e, 2))
      end do
      if (any(real_length <= tiny(1.0_dp))) return
      value = negative_penalized_loglik(self%scaled_length, real_length, x(1:n_edge), self%parent_edge, &
         self%basal_edge, self%lambda)
   end function chronopl_problem_value

   subroutine chronopl_problem_gradient(self, x, gradient)
      class(chronopl_problem), intent(inout) :: self !! Penalized-likelihood dating data attached to the optimizer.
      real(dp), intent(in) :: x(:) !! Edge rates followed by ages of unknown internal nodes.
      real(dp), intent(out) :: gradient(:) !! Analytic gradient of the negative penalized log likelihood.
      real(dp), allocatable :: age(:)
      real(dp), allocatable :: real_length(:)
      real(dp) :: difference
      real(dp) :: mean_basal
      integer :: child_edge
      integer :: e
      integer :: i
      integer :: incoming
      integer :: n_basal
      integer :: n_edge
      integer :: node
      integer :: parent

      gradient = 0.0_dp
      n_edge = self%tree%nedge()
      age = self%fixed_age
      do i = 1, size(self%unknown_node)
         age(self%unknown_node(i)) = x(n_edge + i)
      end do
      allocate(real_length(n_edge))
      do e = 1, n_edge
         real_length(e) = age(self%tree%edge(e, 1)) - age(self%tree%edge(e, 2))
      end do
      if (any(real_length <= tiny(1.0_dp))) return

      gradient(1:n_edge) = real_length - self%scaled_length / x(1:n_edge)
      do e = 1, n_edge
         parent = self%parent_edge(e)
         if (parent <= 0) cycle
         difference = x(e) - x(parent)
         gradient(e) = gradient(e) + 2.0_dp * self%lambda * difference
         gradient(parent) = gradient(parent) - 2.0_dp * self%lambda * difference
      end do
      n_basal = size(self%basal_edge)
      if (n_basal > 1) then
         mean_basal = sum(x(self%basal_edge)) / real(n_basal, dp)
         do i = 1, n_basal
            e = self%basal_edge(i)
            gradient(e) = gradient(e) + 2.0_dp * self%lambda * (x(e) - mean_basal) / real(n_basal - 1, dp)
         end do
      end if

      do i = 1, size(self%unknown_node)
         node = self%unknown_node(i)
         do child_edge = 1, n_edge
            if (self%tree%edge(child_edge, 1) /= node) cycle
            gradient(n_edge + i) = gradient(n_edge + i) + x(child_edge) &
               - self%scaled_length(child_edge) / real_length(child_edge)
         end do
         incoming = 0
         do e = 1, n_edge
            if (self%tree%edge(e, 2) == node) then
               incoming = e
               exit
            end if
         end do
         if (incoming > 0) then
            gradient(n_edge + i) = gradient(n_edge + i) - x(incoming) &
               + self%scaled_length(incoming) / real_length(incoming)
         end if
      end do
   end subroutine chronopl_problem_gradient

   pure real(dp) function negative_penalized_loglik(observed_length, real_length, rates, parent_edge, basal, lambda) &
      result(value)
      !! Returns `-(loglik - lambda * roughness)` for fixed rates and chronological branch lengths.
      real(dp), intent(in) :: observed_length(:) !! Scaled observed substitutions/counts on edges.
      real(dp), intent(in) :: real_length(:) !! Positive chronological duration of every edge.
      real(dp), intent(in) :: rates(:) !! Positive rate for every edge.
      integer, intent(in) :: parent_edge(:) !! Incoming parent-edge index for each edge, zero for basal edges.
      integer, intent(in) :: basal(:) !! Edge indices descending directly from the root.
      real(dp), intent(in) :: lambda !! Nonnegative roughness multiplier.
      real(dp) :: b
      real(dp) :: loglik
      real(dp) :: mean_basal
      real(dp) :: penalty
      integer :: e
      integer :: n_basal

      loglik = 0.0_dp
      do e = 1, size(rates)
         b = rates(e) * real_length(e)
         if (b <= 0.0_dp) then
            value = huge(1.0_dp) / 100.0_dp
            return
         end if
         loglik = loglik - b + observed_length(e) * log(b) - log_gamma(observed_length(e) + 1.0_dp)
      end do
      penalty = 0.0_dp
      do e = 1, size(rates)
         if (parent_edge(e) > 0) penalty = penalty + (rates(e) - rates(parent_edge(e)))**2
      end do
      n_basal = size(basal)
      if (n_basal > 1) then
         mean_basal = sum(rates(basal)) / real(n_basal, dp)
         penalty = penalty + sum((rates(basal) - mean_basal)**2) / real(n_basal - 1, dp)
      end if
      value = -(loglik - lambda * penalty)
   end function negative_penalized_loglik

   subroutine prepare_calibrations(tree, age_min, age_max, calibration_nodes, calibration, lower, upper, age, fixed, info)
      !! Validates/broadcasts chronopl calibration input and creates initial age/fixed-node arrays.
      type(phylo_tree), intent(in) :: tree !! Tree whose internal nodes may be calibrated.
      real(dp), intent(in), optional :: age_min(:) !! Lower/fixed calibration ages; scalar or one per calibrated node.
      real(dp), intent(in), optional :: age_max(:) !! Optional upper calibration ages; scalar or one per calibrated node.
      integer, intent(in), optional :: calibration_nodes(:) !! Internal node numbers; default root.
      integer, allocatable, intent(out) :: calibration(:) !! Validated calibration node numbers.
      real(dp), allocatable, intent(out) :: lower(:) !! Broadcast lower ages matching `calibration`.
      real(dp), allocatable, intent(out) :: upper(:) !! Broadcast upper ages matching `calibration`.
      real(dp), allocatable, intent(out) :: age(:) !! All-node initial age array, tips initialized to zero.
      logical, allocatable, intent(out) :: fixed(:) !! Mask of nodes fixed before optimization.
      integer, intent(out) :: info !! Zero on success or nonzero for invalid calibration dimensions/ranges.
      integer :: ncal

      info = 0
      if (present(calibration_nodes)) then
         ncal = size(calibration_nodes)
         if (ncal < 1) then
            info = 1
            return
         end if
         calibration = calibration_nodes
      else
         ncal = 1
         allocate(calibration(1))
         calibration(1) = tree%root()
      end if
      if (any(calibration <= tree%n_tip) .or. any(calibration > tree%total_nodes())) then
         info = 2
         return
      end if
      call broadcast_real(age_min, ncal, 1.0_dp, lower, info)
      if (info /= 0) then
         info = 3
         return
      end if
      if (present(age_max)) then
         call broadcast_real(age_max, ncal, 0.0_dp, upper, info)
         if (info /= 0) then
            info = 4
            return
         end if
      else
         upper = lower
      end if
      if (any(lower < 0.0_dp) .or. any(upper < lower) .or. .not. all(ieee_is_finite(lower)) &
         .or. .not. all(ieee_is_finite(upper))) then
         info = 5
         return
      end if
      allocate(age(tree%total_nodes()), fixed(tree%total_nodes()))
      age = 0.0_dp
      fixed = .false.
      fixed(1:tree%n_tip) = .true.
   end subroutine prepare_calibrations

   subroutine initialize_ages(tree, calibration, lower, upper, age, info)
      !! Reproduces chronopl's root-to-tip interpolation used for positive starting chronological edge lengths.
      type(phylo_tree), intent(in) :: tree !! Rooted tree topology used to form root-to-tip paths.
      integer, intent(in) :: calibration(:) !! Internal calibration nodes.
      real(dp), intent(in) :: lower(:) !! Calibration lower/fixed ages.
      real(dp), intent(in) :: upper(:) !! Calibration upper/fixed ages.
      real(dp), intent(inout) :: age(:) !! All-node ages; tips are zero and internal ages are filled on output.
      integer, intent(out) :: info !! Zero on success or nonzero if interpolation cannot create positive edge lengths.
      integer, allocatable :: order(:)
      integer, allocatable :: parent(:)
      integer, allocatable :: path(:)
      logical, allocatable :: known(:)
      integer, allocatable :: known_count(:)
      real(dp) :: step
      integer :: i
      integer :: j
      integer :: k
      integer :: length
      integer :: root
      integer :: tip

      info = 0
      parent = parent_vector(tree)
      allocate(known(tree%total_nodes()))
      known = .false.
      known(1:tree%n_tip) = .true.
      do i = 1, size(calibration)
         age(calibration(i)) = 0.5_dp * (lower(i) + upper(i))
         known(calibration(i)) = .true.
      end do
      root = tree%root()
      if (.not. known(root)) then
         age(root) = 3.0_dp * maxval(upper)
         known(root) = .true.
      end if
      allocate(known_count(tree%n_tip), order(tree%n_tip))
      do tip = 1, tree%n_tip
         call root_to_tip_path(parent, root, tip, path)
         known_count(tip) = count(known(path))
         order(tip) = tip
      end do
      call sort_tips_by_known(order, known_count)
      do k = 1, tree%n_tip
         tip = order(k)
         call root_to_tip_path(parent, root, tip, path)
         length = size(path)
         i = 2
         do while (i <= length)
            if (known(path(i))) then
               i = i + 1
               cycle
            end if
            j = i + 1
            do while (j <= length .and. .not. known(path(j)))
               j = j + 1
            end do
            if (j > length) then
               info = 1
               return
            end if
            step = (age(path(i - 1)) - age(path(j))) / real(j - i + 1, dp)
            do while (i < j)
               age(path(i)) = age(path(i - 1)) - step
               known(path(i)) = .true.
               i = i + 1
            end do
            i = j + 1
         end do
      end do
      do i = 1, tree%nedge()
         if (age(tree%edge(i, 1)) - age(tree%edge(i, 2)) <= 0.0_dp) then
            info = 2
            return
         end if
      end do
   end subroutine initialize_ages

   subroutine mark_fixed_calibrations(tree, calibration, lower, upper, age, fixed)
      !! Marks fixed calibration nodes and resets their exact ages after initialization.
      type(phylo_tree), intent(in) :: tree !! Tree defining valid node numbers.
      integer, intent(in) :: calibration(:) !! Internal calibrated nodes.
      real(dp), intent(in) :: lower(:) !! Calibration lower ages.
      real(dp), intent(in) :: upper(:) !! Calibration upper ages.
      real(dp), intent(inout) :: age(:) !! All-node initial ages updated for exact calibrations.
      logical, intent(inout) :: fixed(:) !! All-node fixed mask with tips already fixed.
      integer :: i

      fixed(1:tree%n_tip) = .true.
      do i = 1, size(calibration)
         if (abs(lower(i) - upper(i)) <= 0.0_dp) then
            age(calibration(i)) = lower(i)
            fixed(calibration(i)) = .true.
         end if
      end do
   end subroutine mark_fixed_calibrations

   pure subroutine collect_unknown_nodes(tree, fixed, unknown)
      !! Collects internal nodes whose ages remain optimization variables.
      type(phylo_tree), intent(in) :: tree !! Tree defining the internal-node number range.
      logical, intent(in) :: fixed(:) !! All-node mask of tips and exactly calibrated nodes.
      integer, allocatable, intent(out) :: unknown(:) !! Internal ape node numbers not fixed by calibration.
      integer :: i
      integer :: k

      allocate(unknown(count(.not. fixed(tree%n_tip + 1:tree%total_nodes()))))
      k = 0
      do i = tree%n_tip + 1, tree%total_nodes()
         if (fixed(i)) cycle
         k = k + 1
         unknown(k) = i
      end do
   end subroutine collect_unknown_nodes

   pure subroutine build_edge_relations(tree, parent_edge, basal)
      !! Builds the parent-edge map and list of edges descending directly from the root.
      type(phylo_tree), intent(in) :: tree !! Rooted tree whose edge adjacency is required.
      integer, allocatable, intent(out) :: parent_edge(:) !! Parent edge for each edge, or zero for basal edges.
      integer, allocatable, intent(out) :: basal(:) !! Indices of edges whose parent node is the root.
      integer :: e
      integer :: j
      integer :: k
      integer :: root

      root = tree%root()
      allocate(parent_edge(tree%nedge()), basal(count(tree%edge(:, 1) == root)))
      parent_edge = 0
      k = 0
      do e = 1, tree%nedge()
         if (tree%edge(e, 1) == root) then
            k = k + 1
            basal(k) = e
         end if
         do j = 1, tree%nedge()
            if (tree%edge(j, 2) == tree%edge(e, 1)) then
               parent_edge(e) = j
               exit
            end if
         end do
      end do
   end subroutine build_edge_relations

   pure subroutine apply_calibration_bound(node, calibration, lower_calibration, upper_calibration, lower, upper)
      !! Replaces a free node's generic optimizer bounds with any interval calibration on that node.
      integer, intent(in) :: node !! Internal ape node number whose optimization bounds are requested.
      integer, intent(in) :: calibration(:) !! Calibrated node numbers.
      real(dp), intent(in) :: lower_calibration(:) !! Calibration lower ages.
      real(dp), intent(in) :: upper_calibration(:) !! Calibration upper ages.
      real(dp), intent(inout) :: lower !! Current lower age bound, replaced on a matching calibration.
      real(dp), intent(inout) :: upper !! Current upper age bound, replaced on a matching calibration.
      integer :: i

      do i = 1, size(calibration)
         if (calibration(i) /= node) cycle
         if (abs(lower_calibration(i) - upper_calibration(i)) > 0.0_dp) then
            lower = lower_calibration(i)
            upper = upper_calibration(i)
         end if
         return
      end do
   end subroutine apply_calibration_bound

   subroutine broadcast_real(values, n, default_value, output, info)
      !! Broadcasts an optional scalar/vector real argument to a required calibration length.
      real(dp), intent(in), optional :: values(:) !! Optional scalar or length-`n` values to broadcast.
      integer, intent(in) :: n !! Required output length.
      real(dp), intent(in) :: default_value !! Value used when the optional input is absent.
      real(dp), allocatable, intent(out) :: output(:) !! Broadcast vector of length `n`.
      integer, intent(out) :: info !! Zero on success or one for an incompatible input length.

      info = 0
      allocate(output(n))
      if (.not. present(values)) then
         output = default_value
      else if (size(values) == 1) then
         output = values(1)
      else if (size(values) == n) then
         output = values
      else
         output = default_value
         info = 1
      end if
   end subroutine broadcast_real

   pure subroutine root_to_tip_path(parent, root, tip, path)
      !! Constructs one root-to-tip node path from a parent vector.
      integer, intent(in) :: parent(:) !! Parent node number for every node, with zero at the root.
      integer, intent(in) :: root !! Root ape node number.
      integer, intent(in) :: tip !! Terminal ape node number.
      integer, allocatable, intent(out) :: path(:) !! Node sequence from root through the requested tip.
      integer, allocatable :: reverse_path(:)
      integer :: current
      integer :: length

      allocate(reverse_path(size(parent)))
      length = 1
      reverse_path(1) = tip
      current = tip
      do while (current /= root)
         current = parent(current)
         if (current == 0) exit
         length = length + 1
         reverse_path(length) = current
      end do
      allocate(path(length))
      path = reverse_path(length:1:-1)
   end subroutine root_to_tip_path

   pure subroutine sort_tips_by_known(order, known_count)
      !! Stable insertion-sorts tip indices by decreasing initial calibrated-node count.
      integer, intent(inout) :: order(:) !! Tip indices reordered in place.
      integer, intent(in) :: known_count(:) !! Initial number of known ages along each corresponding tip path.
      integer :: i
      integer :: j
      integer :: key

      do i = 2, size(order)
         key = order(i)
         j = i - 1
         do while (j >= 1)
            if (known_count(order(j)) >= known_count(key)) exit
            order(j + 1) = order(j)
            j = j - 1
         end do
         order(j + 1) = key
      end do
   end subroutine sort_tips_by_known

end module ape_chronopl
