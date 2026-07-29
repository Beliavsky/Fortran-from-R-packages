! RiskPortfolios Fortran, derived from RiskPortfolios 2.1.7.
! Original code Copyright (C) 2013-2021 David Ardia.
! Original authors: David Ardia, Kris Boudt, Jean-Philippe Gagnon-Fleury.
! SPDX-License-Identifier: GPL-2.0-or-later
module riskportfolios_portfolio
   use riskportfolios_kinds, only : dp
   use riskportfolios_linalg, only : solve_linear
   use riskportfolios_stats, only : covariance_to_correlation, quantile_type7, median_value
   use riskportfolios_optimization, only : projected_gradient, project_feasible, &
      OBJECTIVE_MEAN_VARIANCE, OBJECTIVE_MINIMUM_VARIANCE, OBJECTIVE_ERC, &
      OBJECTIVE_MAXIMUM_DIVERSIFICATION, OBJECTIVE_RISK_EFFICIENT, &
      OBJECTIVE_MAXIMUM_DECORRELATION
   implicit none
   private

   integer, parameter, public :: PORT_MEAN_VARIANCE = 1
   integer, parameter, public :: PORT_MINIMUM_VARIANCE = 2
   integer, parameter, public :: PORT_INVERSE_VOLATILITY = 3
   integer, parameter, public :: PORT_EQUAL_RISK_CONTRIBUTION = 4
   integer, parameter, public :: PORT_MAXIMUM_DIVERSIFICATION = 5
   integer, parameter, public :: PORT_RISK_EFFICIENT = 6
   integer, parameter, public :: PORT_MAXIMUM_DECORRELATION = 7

   integer, parameter, public :: CONSTRAINT_NONE = 1
   integer, parameter, public :: CONSTRAINT_LONG_ONLY = 2
   integer, parameter, public :: CONSTRAINT_GROSS = 3
   integer, parameter, public :: CONSTRAINT_USER = 4

   type, public :: portfolio_control
      integer :: constraint = CONSTRAINT_NONE
      real(dp) :: gross_limit = 1.6_dp
      real(dp) :: gamma = 0.8773_dp
      integer :: max_iterations = 4000
      real(dp) :: tolerance = 1.0e-10_dp
      real(dp), allocatable :: lower_bounds(:)
      real(dp), allocatable :: upper_bounds(:)
      real(dp), allocatable :: initial_weights(:)
   end type portfolio_control

   public :: optimal_portfolio, portfolio_risk_contributions

contains

   subroutine optimal_portfolio(sigma, weights, portfolio_type, mu, semidev, control, info)
      real(dp), intent(in) :: sigma(:, :)
      real(dp), allocatable, intent(out) :: weights(:)
      integer, intent(in), optional :: portfolio_type
      real(dp), intent(in), optional :: mu(:), semidev(:)
      type(portfolio_control), intent(in), optional :: control
      integer, intent(out), optional :: info
      type(portfolio_control) :: ctr
      real(dp), allocatable :: lower(:), upper(:), initial(:), solution(:), projected(:)
      real(dp), allocatable :: rhs_solution(:), target(:)
      real(dp) :: one(size(sigma, 1)), sd(size(sigma, 1))
      real(dp) :: rho(size(sigma, 1), size(sigma, 2)), denominator
      integer :: n, ptype, stat
      logical :: use_gross

      n = size(sigma, 1)
      allocate(weights(n))
      weights = 0.0_dp
      if (size(sigma, 2) /= n .or. n == 0) then
         if (present(info)) info = -1
         return
      end if
      ctr = portfolio_control()
      if (present(control)) ctr = control
      ptype = PORT_MEAN_VARIANCE
      if (present(portfolio_type)) ptype = portfolio_type

      call prepare_bounds(n, ptype, ctr, lower, upper, initial, use_gross, stat)
      if (stat /= 0) then
         if (present(info)) info = stat
         return
      end if
      one = 1.0_dp

      select case (ptype)
      case (PORT_MEAN_VARIANCE)
         if (.not. present(mu) .or. size(mu) /= n) then
            stat = -10
         else if (ctr%constraint == CONSTRAINT_NONE) then
            call solve_linear(sigma, mu, rhs_solution, stat)
            if (stat == 0) then
               denominator = sum(rhs_solution)
               if (abs(denominator) > tiny(1.0_dp)) then
                  weights = rhs_solution / denominator
               else
                  stat = -11
               end if
            end if
         else
            allocate(solution(n))
            call projected_gradient(OBJECTIVE_MEAN_VARIANCE, sigma, mu, ctr%gamma, &
               initial, lower, upper, use_gross, ctr%gross_limit, &
               ctr%max_iterations, ctr%tolerance, solution, stat)
            weights = solution
         end if

      case (PORT_MINIMUM_VARIANCE)
         if (ctr%constraint == CONSTRAINT_NONE) then
            call solve_linear(sigma, one, rhs_solution, stat)
            if (stat == 0) then
               denominator = sum(rhs_solution)
               if (abs(denominator) > tiny(1.0_dp)) then
                  weights = rhs_solution / denominator
               else
                  stat = -12
               end if
            end if
         else
            allocate(solution(n))
            call projected_gradient(OBJECTIVE_MINIMUM_VARIANCE, sigma, one, ctr%gamma, &
               initial, lower, upper, use_gross, ctr%gross_limit, &
               ctr%max_iterations, ctr%tolerance, solution, stat)
            weights = solution
         end if

      case (PORT_INVERSE_VOLATILITY)
         sd = 0.0_dp
         block
            integer :: i
            do i = 1, n
               sd(i) = sqrt(max(sigma(i, i), tiny(1.0_dp)))
            end do
         end block
         weights = 1.0_dp / sd
         weights = weights / sum(weights)
         stat = 0

      case (PORT_EQUAL_RISK_CONTRIBUTION)
         allocate(solution(n))
         call projected_gradient(OBJECTIVE_ERC, sigma, one, ctr%gamma, initial, lower, &
            upper, use_gross, ctr%gross_limit, ctr%max_iterations, ctr%tolerance, &
            solution, stat)
         weights = solution

      case (PORT_MAXIMUM_DIVERSIFICATION)
         allocate(solution(n))
         call projected_gradient(OBJECTIVE_MAXIMUM_DIVERSIFICATION, sigma, one, &
            ctr%gamma, initial, lower, upper, use_gross, ctr%gross_limit, &
            ctr%max_iterations, ctr%tolerance, solution, stat)
         weights = solution

      case (PORT_RISK_EFFICIENT)
         if (.not. present(semidev) .or. size(semidev) /= n) then
            stat = -13
         else
            allocate(target(n), solution(n))
            call risk_efficiency_scores(semidev, target)
            call projected_gradient(OBJECTIVE_RISK_EFFICIENT, sigma, target, ctr%gamma, &
               initial, lower, upper, use_gross, ctr%gross_limit, &
               ctr%max_iterations, ctr%tolerance, solution, stat)
            weights = solution
         end if

      case (PORT_MAXIMUM_DECORRELATION)
         rho = covariance_to_correlation(sigma)
         if (ctr%constraint == CONSTRAINT_NONE) then
            call solve_linear(rho, one, rhs_solution, stat)
            if (stat == 0) then
               denominator = sum(rhs_solution)
               if (abs(denominator) > tiny(1.0_dp)) then
                  weights = rhs_solution / denominator
               else
                  stat = -14
               end if
            end if
         else
            allocate(solution(n))
            call projected_gradient(OBJECTIVE_MAXIMUM_DECORRELATION, rho, one, ctr%gamma, &
               initial, lower, upper, use_gross, ctr%gross_limit, &
               ctr%max_iterations, ctr%tolerance, solution, stat)
            weights = solution
         end if

      case default
         stat = -2
      end select

      if (stat == 0 .and. ptype /= PORT_INVERSE_VOLATILITY) then
         allocate(projected(n))
         call project_feasible(weights, lower, upper, use_gross, ctr%gross_limit, &
            projected, stat)
         weights = projected
      end if
      if (present(info)) info = stat


   end subroutine optimal_portfolio

   subroutine prepare_bounds(n, portfolio_type, ctr, lower, upper, initial, &
      use_gross, info)
      integer, intent(in) :: n, portfolio_type
      type(portfolio_control), intent(in) :: ctr
      real(dp), allocatable, intent(out) :: lower(:), upper(:), initial(:)
      logical, intent(out) :: use_gross
      integer, intent(out) :: info
      real(dp), parameter :: broad_bound = 1.0e6_dp
      integer :: stat
      real(dp) :: projected_initial(n)

      allocate(lower(n), upper(n), initial(n))
      use_gross = ctr%constraint == CONSTRAINT_GROSS
      select case (ctr%constraint)
      case (CONSTRAINT_NONE)
         lower = -broad_bound
         upper = broad_bound
      case (CONSTRAINT_LONG_ONLY)
         lower = 0.0_dp
         upper = 1.0_dp
      case (CONSTRAINT_GROSS)
         lower = -ctr%gross_limit
         upper = ctr%gross_limit
      case (CONSTRAINT_USER)
         if (.not. allocated(ctr%lower_bounds) .or. &
             .not. allocated(ctr%upper_bounds)) then
            info = -20
            return
         end if
         if (size(ctr%lower_bounds) /= n .or. size(ctr%upper_bounds) /= n) then
            info = -21
            return
         end if
         lower = ctr%lower_bounds
         upper = ctr%upper_bounds
      case default
         info = -22
         return
      end select

      if (allocated(ctr%lower_bounds) .and. ctr%constraint /= CONSTRAINT_USER) then
         if (size(ctr%lower_bounds) == n) lower = ctr%lower_bounds
      end if
      if (allocated(ctr%upper_bounds) .and. ctr%constraint /= CONSTRAINT_USER) then
         if (size(ctr%upper_bounds) == n) upper = ctr%upper_bounds
      end if

      if (portfolio_type == PORT_RISK_EFFICIENT) then
         if (.not. allocated(ctr%lower_bounds)) lower = 1.0_dp / real(2 * n, dp)
         if (.not. allocated(ctr%upper_bounds)) upper = 2.0_dp / real(n, dp)
      end if

      if (allocated(ctr%initial_weights)) then
         if (size(ctr%initial_weights) /= n) then
            info = -23
            return
         end if
         initial = ctr%initial_weights
      else if (all(lower > -0.5_dp * broad_bound) .and. &
               all(upper < 0.5_dp * broad_bound)) then
         initial = 0.5_dp * (lower + upper)
      else
         initial = 1.0_dp / real(n, dp)
      end if
      call project_feasible(initial, lower, upper, use_gross, ctr%gross_limit, &
         projected_initial, stat)
      initial = projected_initial
      info = stat
   end subroutine prepare_bounds

   subroutine risk_efficiency_scores(semidev, scores)
      real(dp), intent(in) :: semidev(:)
      real(dp), intent(out) :: scores(:)
      real(dp) :: cutoffs(11), values(size(semidev))
      logical :: selected(size(semidev))
      integer :: i, n_selected

      cutoffs(1) = 0.0_dp
      do i = 1, 10
         cutoffs(i + 1) = quantile_type7(semidev, real(i, dp) / 10.0_dp)
      end do
      scores = semidev
      do i = 1, 10
         selected = semidev > cutoffs(i) .and. semidev <= cutoffs(i + 1)
         n_selected = count(selected)
         if (n_selected > 0) then
            values = 0.0_dp
            values(1:n_selected) = pack(semidev, selected)
            where (selected)
               scores = median_value(values(1:n_selected))
            end where
         end if
      end do
   end subroutine risk_efficiency_scores

   pure function portfolio_risk_contributions(sigma, weights, proportional) result(rc)
      real(dp), intent(in) :: sigma(:, :), weights(:)
      logical, intent(in), optional :: proportional
      real(dp) :: rc(size(weights))
      real(dp) :: sigma_w(size(weights)), variance, portfolio_sd
      logical :: use_proportional

      sigma_w = matmul(sigma, weights)
      variance = max(dot_product(weights, sigma_w), tiny(1.0_dp))
      use_proportional = .true.
      if (present(proportional)) use_proportional = proportional
      if (use_proportional) then
         rc = weights * sigma_w / variance
      else
         portfolio_sd = sqrt(variance)
         rc = weights * sigma_w / portfolio_sd
      end if
   end function portfolio_risk_contributions

end module riskportfolios_portfolio
