! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2012-2014 Alexios Galanos and Bernhard Pfaff.
! Modern Fortran translation Copyright (C) 2026 OpenAI.
! This derivative work is distributed under GPL-3.0-or-later.
! Native convex QP support replacing the upstream quadprog binding.
module parma_qp
   use parma_kinds, only: dp
   use parma_types, only: qp_result
   use parma_linalg, only: project_box_budget, vector_norm
   implicit none
   private
   public :: qp_box_budget

contains

   subroutine qp_box_budget(h,g,lb,ub,budget,result,x0,target_a,target_b, &
      target_equality,max_iter,tol)
      real(dp), intent(in) :: h(:,:),g(:),lb(:),ub(:),budget
      type(qp_result), intent(out) :: result
      real(dp), intent(in), optional :: x0(:),target_a(:),target_b
      logical, intent(in), optional :: target_equality
      integer, intent(in), optional :: max_iter
      real(dp), intent(in), optional :: tol
      real(dp), allocatable :: x(:),trial(:),projected(:),gradient(:)
      real(dp) :: f,ftrial,step,epsilonx,reward,denom,delta
      integer :: n,iter,maxit,info
      logical :: target_eq

      n = size(g)
      maxit = 5000
      if (present(max_iter)) maxit = max_iter
      epsilonx = 1.0e-10_dp
      if (present(tol)) epsilonx = tol
      target_eq = .false.
      if (present(target_equality)) target_eq = target_equality
      allocate(x(n),trial(n),projected(n),gradient(n),result%x(n))
      x = budget/real(n,dp)
      if (present(x0)) x = x0
      call project_box_budget(x,lb,ub,budget,trial,info)
      if (info /= 0) then
         result%status = 2
         result%message = 'infeasible box and budget constraints'
         result%x = trial
         return
      end if
      x = trial
      step = 1.0_dp/max(1.0_dp,maxval(abs(h)))
      f = objective(x)
      do iter = 1,maxit
         gradient = matmul(h,x)-g
         trial = x-step*gradient
         call project_box_budget(trial,lb,ub,budget,projected,info)
         trial = projected
         if (present(target_a) .and. present(target_b)) then
            reward = dot_product(target_a,trial)
            denom = dot_product(target_a-sum(target_a)/real(n,dp),target_a)
            if (abs(denom) > tiny(1.0_dp)) then
               if (target_eq .or. reward < target_b) then
                  delta = (target_b-reward)/denom
                  trial = trial+delta*(target_a-sum(target_a)/real(n,dp))
                  call project_box_budget(trial,lb,ub,budget,projected,info)
                  trial = projected
               end if
            end if
         end if
         ftrial = objective(trial)
         if (ftrial <= f) then
            if (vector_norm(trial-x) <= epsilonx*(1.0_dp+vector_norm(x))) then
               x = trial
               exit
            end if
            x = trial
            f = ftrial
            step = min(1.1_dp*step,1.0_dp)
         else
            step = 0.5_dp*step
            if (step < 1.0e-15_dp) exit
         end if
      end do
      result%x = x
      result%objective = objective(x)
      result%iterations = min(iter,maxit)
      result%status = 0
      result%message = 'converged'

   contains
      function objective(z) result(value)
         real(dp), intent(in) :: z(:)
         real(dp) :: value
         value = 0.5_dp*dot_product(z,matmul(h,z))-dot_product(g,z)
      end function objective
   end subroutine qp_box_budget

end module parma_qp
