! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2012-2014 Alexios Galanos and Bernhard Pfaff.
! Modern Fortran translation Copyright (C) 2026 OpenAI.
! This derivative work is distributed under GPL-3.0-or-later.
! Derived from parma 1.7, Copyright (C) 2012-2014 Alexios Galanos and Bernhard Pfaff.
module parma_constraints
   use parma_kinds, only: dp
   use parma_types, only: parma_spec
   use parma_risk, only: quadratic_variance
   implicit none
   private
   public :: budget_residual, target_residual, turnover_value
   public :: buy_turnover_value, sell_turnover_value, linear_constraint_penalty
   public :: constraint_penalty, constraints_feasible
   public :: inequality_turnover, inequality_buy_sell_turnover
   public :: inequality_variance

contains

   function budget_residual(weights, budget) result(value)
      real(dp), intent(in) :: weights(:), budget
      real(dp) :: value
      value = sum(weights) - budget
   end function budget_residual

   function target_residual(weights, mu, target) result(value)
      real(dp), intent(in) :: weights(:), mu(:), target
      real(dp) :: value
      value = dot_product(weights,mu) - target
   end function target_residual

   function turnover_value(weights, initial) result(value)
      real(dp), intent(in) :: weights(:), initial(:)
      real(dp) :: value
      value = sum(abs(weights-initial))
   end function turnover_value

   function buy_turnover_value(weights, initial) result(value)
      real(dp), intent(in) :: weights(:), initial(:)
      real(dp) :: value
      value = sum(max(weights-initial,0.0_dp))
   end function buy_turnover_value

   function sell_turnover_value(weights, initial) result(value)
      real(dp), intent(in) :: weights(:), initial(:)
      real(dp) :: value
      value = sum(max(initial-weights,0.0_dp))
   end function sell_turnover_value

   function inequality_turnover(weights, initial, limit) result(value)
      real(dp), intent(in) :: weights(:), initial(:), limit
      real(dp) :: value
      value = turnover_value(weights,initial)-limit
   end function inequality_turnover

   function inequality_buy_sell_turnover(weights, initial, buy_limit, sell_limit) result(values)
      real(dp), intent(in) :: weights(:), initial(:), buy_limit, sell_limit
      real(dp) :: values(2)
      values(1) = buy_turnover_value(weights,initial)-buy_limit
      values(2) = sell_turnover_value(weights,initial)-sell_limit
   end function inequality_buy_sell_turnover

   function inequality_variance(weights, covariance, limit) result(value)
      real(dp), intent(in) :: weights(:), covariance(:,:), limit
      real(dp) :: value
      value = quadratic_variance(weights,covariance)-limit
   end function inequality_variance

   function linear_constraint_penalty(weights, eq_a, eq_b, ineq_a, ineq_lb, ineq_ub) result(value)
      real(dp), intent(in) :: weights(:)
      real(dp), intent(in), optional :: eq_a(:,:), eq_b(:)
      real(dp), intent(in), optional :: ineq_a(:,:), ineq_lb(:), ineq_ub(:)
      real(dp) :: value
      real(dp), allocatable :: ax(:)

      value = 0.0_dp
      if (present(eq_a) .and. present(eq_b)) then
         if (size(eq_a,1) > 0) then
            ax = matmul(eq_a,weights)-eq_b
            value = value + sum(ax*ax)
         end if
      end if
      if (present(ineq_a)) then
         if (size(ineq_a,1) > 0) then
            ax = matmul(ineq_a,weights)
            if (present(ineq_lb)) value = value + sum(max(ineq_lb-ax,0.0_dp)**2)
            if (present(ineq_ub)) value = value + sum(max(ax-ineq_ub,0.0_dp)**2)
         end if
      end if
   end function linear_constraint_penalty

   function constraint_penalty(weights, spec) result(value)
      real(dp), intent(in) :: weights(:)
      type(parma_spec), intent(in) :: spec
      real(dp) :: value, reward

      if (spec%leverage > 0.0_dp) then
         value = (sum(sqrt(weights*weights+1.0e-20_dp))-spec%leverage)**2
      else
         value = budget_residual(weights,spec%budget)**2
      end if
      if (allocated(spec%lb)) value = value + sum(max(spec%lb-weights,0.0_dp)**2)
      if (allocated(spec%ub)) value = value + sum(max(weights-spec%ub,0.0_dp)**2)
      if (allocated(spec%mu)) then
         reward = dot_product(weights,spec%mu)
         if (spec%target_is_equality) then
            value = value + (reward-spec%target)**2
         else if (abs(spec%target) > tiny(1.0_dp)) then
            value = value + max(spec%target-reward,0.0_dp)**2
         end if
      end if
      if (allocated(spec%eq_a) .and. allocated(spec%eq_b)) then
         value = value + linear_constraint_penalty(weights,eq_a=spec%eq_a,eq_b=spec%eq_b)
      end if
      if (allocated(spec%ineq_a)) then
         if (allocated(spec%ineq_lb) .and. allocated(spec%ineq_ub)) then
            value = value + linear_constraint_penalty(weights,ineq_a=spec%ineq_a, &
               ineq_lb=spec%ineq_lb,ineq_ub=spec%ineq_ub)
         else if (allocated(spec%ineq_lb)) then
            value = value + linear_constraint_penalty(weights,ineq_a=spec%ineq_a, &
               ineq_lb=spec%ineq_lb)
         else if (allocated(spec%ineq_ub)) then
            value = value + linear_constraint_penalty(weights,ineq_a=spec%ineq_a, &
               ineq_ub=spec%ineq_ub)
         end if
      end if
      if (allocated(spec%initial)) then
         if (spec%turnover_limit < huge(1.0_dp)/2.0_dp) then
            value = value + max(turnover_value(weights,spec%initial)-spec%turnover_limit,0.0_dp)**2
         end if
         if (spec%buy_turnover_limit < huge(1.0_dp)/2.0_dp) then
            value = value + max(buy_turnover_value(weights,spec%initial)-spec%buy_turnover_limit,0.0_dp)**2
         end if
         if (spec%sell_turnover_limit < huge(1.0_dp)/2.0_dp) then
            value = value + max(sell_turnover_value(weights,spec%initial)-spec%sell_turnover_limit,0.0_dp)**2
         end if
      end if
      if (allocated(spec%cov) .and. spec%variance_limit < huge(1.0_dp)/2.0_dp) then
         value = value + max(quadratic_variance(weights,spec%cov)-spec%variance_limit,0.0_dp)**2
      end if
      if (spec%max_positions > 0) then
         value = value+real(max(count(abs(weights)>1.0e-6_dp)-spec%max_positions,0),dp)**2
      end if
   end function constraint_penalty

   function constraints_feasible(weights, spec, tol) result(feasible)
      real(dp), intent(in) :: weights(:)
      type(parma_spec), intent(in) :: spec
      real(dp), intent(in), optional :: tol
      logical :: feasible
      real(dp) :: epsilonx, reward
      real(dp), allocatable :: ax(:)

      epsilonx = 1.0e-7_dp
      if (present(tol)) epsilonx = tol
      if (spec%leverage > 0.0_dp) then
         feasible = abs(sum(abs(weights))-spec%leverage) <= epsilonx
      else
         feasible = abs(sum(weights)-spec%budget) <= epsilonx
      end if
      if (.not. feasible) return
      if (allocated(spec%lb)) feasible = feasible .and. all(weights >= spec%lb-epsilonx)
      if (allocated(spec%ub)) feasible = feasible .and. all(weights <= spec%ub+epsilonx)
      if (.not. feasible) return
      if (allocated(spec%mu)) then
         reward = dot_product(weights,spec%mu)
         if (spec%target_is_equality) then
            feasible = abs(reward-spec%target) <= epsilonx
         else if (abs(spec%target) > tiny(1.0_dp)) then
            feasible = reward >= spec%target-epsilonx
         end if
      end if
      if (.not. feasible) return
      if (allocated(spec%eq_a) .and. allocated(spec%eq_b)) then
         feasible = all(abs(matmul(spec%eq_a,weights)-spec%eq_b) <= epsilonx)
      end if
      if (.not. feasible) return
      if (allocated(spec%ineq_a)) then
         ax = matmul(spec%ineq_a,weights)
         if (allocated(spec%ineq_lb)) feasible = feasible .and. all(ax >= spec%ineq_lb-epsilonx)
         if (allocated(spec%ineq_ub)) feasible = feasible .and. all(ax <= spec%ineq_ub+epsilonx)
      end if
      if (.not. feasible) return
      if (allocated(spec%initial)) then
         feasible = feasible .and. turnover_value(weights,spec%initial) <= spec%turnover_limit+epsilonx
         feasible = feasible .and. buy_turnover_value(weights,spec%initial) <= spec%buy_turnover_limit+epsilonx
         feasible = feasible .and. sell_turnover_value(weights,spec%initial) <= spec%sell_turnover_limit+epsilonx
      end if
      if (.not. feasible) return
      if (allocated(spec%cov)) then
         feasible = feasible .and. quadratic_variance(weights,spec%cov) <= spec%variance_limit+epsilonx
      end if
      if (spec%max_positions > 0) then
         feasible = feasible .and. count(abs(weights)>epsilonx) <= spec%max_positions
      end if
   end function constraints_feasible

end module parma_constraints
