! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2012-2014 Alexios Galanos and Bernhard Pfaff.
! Modern Fortran translation Copyright (C) 2026 OpenAI.
! This derivative work is distributed under GPL-3.0-or-later.
! The original package wraps the Lobo-Vandenberghe-Boyd SOCP code.
! This module supplies a native Fortran logarithmic-barrier implementation.
module parma_socp
   use parma_kinds, only: dp
   use parma_types, only: socp_result
   use parma_linalg, only: solve_linear, vector_norm
   implicit none
   private
   public :: socp_solve, socp_max_violation, socp_feasible

contains

   subroutine socp_solve(f,a,b,c,d,cone_sizes,result,x0,max_iter,tol)
      real(dp), intent(in) :: f(:),a(:,:),b(:),c(:,:),d(:)
      integer, intent(in) :: cone_sizes(:)
      type(socp_result), intent(out) :: result
      real(dp), intent(in), optional :: x0(:)
      integer, intent(in), optional :: max_iter
      real(dp), intent(in), optional :: tol
      real(dp), allocatable :: x(:),g(:),h(:,:),dx(:),trial(:)
      real(dp) :: t, mu, epsilonx, decrement, alpha, phi, phi_trial
      real(dp) :: violation
      integer :: n, iter, outer, maxit, info

      n = size(f)
      maxit = 200
      if (present(max_iter)) maxit = max_iter
      epsilonx = 1.0e-8_dp
      if (present(tol)) epsilonx = tol
      allocate(x(n),g(n),h(n,n),dx(n),trial(n),result%x(n))
      x = 0.0_dp
      if (present(x0)) x = x0
      call phase1_socp(a,b,c,d,cone_sizes,x,maxit,epsilonx,info)
      if (info /= 0) then
         result%x = x
         result%max_violation = socp_max_violation(x,a,b,c,d,cone_sizes)
         result%status = 2
         result%message = 'failed to obtain a strictly feasible point'
         return
      end if
      t = 1.0_dp
      mu = 10.0_dp
      iter = 0
      do outer = 1, 40
         do
            iter = iter+1
            call barrier_derivatives(x,t,f,a,b,c,d,cone_sizes,g,h,phi,info)
            if (info /= 0) exit
            call solve_linear(h,-g,dx,info)
            if (info /= 0) exit
            decrement = -dot_product(g,dx)
            if (decrement/2.0_dp <= epsilonx) exit
            alpha = 1.0_dp
            do
               trial = x+alpha*dx
               if (socp_feasible(trial,a,b,c,d,cone_sizes,strict=.true.)) then
                  call barrier_value(trial,t,f,a,b,c,d,cone_sizes,phi_trial,info)
                  if (info == 0 .and. phi_trial <= phi+0.01_dp*alpha*dot_product(g,dx)) exit
               end if
               alpha = 0.5_dp*alpha
               if (alpha < 1.0e-14_dp) exit
            end do
            if (alpha < 1.0e-14_dp) exit
            x = trial
            if (iter >= maxit) exit
         end do
         if (real(size(cone_sizes),dp)/t <= epsilonx) exit
         t = t*mu
         if (iter >= maxit) exit
      end do
      violation = socp_max_violation(x,a,b,c,d,cone_sizes)
      result%x = x
      result%objective = dot_product(f,x)
      result%max_violation = violation
      result%iterations = iter
      if (violation <= sqrt(epsilonx)) then
         result%status = 0
         result%message = 'converged'
      else
         result%status = 1
         result%message = 'iteration limit or residual violation'
      end if
   end subroutine socp_solve

   subroutine phase1_socp(a,b,c,d,cone_sizes,x,max_iter,tol,info)
      real(dp), intent(in) :: a(:,:),b(:),c(:,:),d(:),tol
      integer, intent(in) :: cone_sizes(:),max_iter
      real(dp), intent(inout) :: x(:)
      integer, intent(out) :: info
      real(dp), allocatable :: gradient(:),trial(:)
      real(dp) :: objective,trial_obj,step,violation,v,u_norm
      integer :: iter,cone,row0,row1

      allocate(gradient(size(x)),trial(size(x)))
      info = 1
      do iter = 1, max_iter
         objective = 0.0_dp
         gradient = 0.0_dp
         row0 = 1
         do cone = 1, size(cone_sizes)
            row1 = row0+cone_sizes(cone)-1
            u_norm = vector_norm(matmul(a(row0:row1,:),x)+b(row0:row1))
            v = dot_product(c(cone,:),x)+d(cone)
            violation = u_norm-v+1.0e-6_dp
            if (violation > 0.0_dp) then
               objective = objective+0.5_dp*violation*violation
               if (u_norm > tiny(1.0_dp)) then
                  gradient = gradient+violation*(matmul(transpose(a(row0:row1,:)), &
                     matmul(a(row0:row1,:),x)+b(row0:row1))/u_norm-c(cone,:))
               else
                  gradient = gradient-violation*c(cone,:)
               end if
            end if
            row0 = row1+1
         end do
         if (socp_feasible(x,a,b,c,d,cone_sizes,strict=.true.)) then
            info = 0
            return
         end if
         step = 1.0_dp/max(1.0_dp,vector_norm(gradient))
         do
            trial = x-step*gradient
            trial_obj = phase1_objective(trial,a,b,c,d,cone_sizes)
            if (trial_obj <= objective-1.0e-4_dp*step*dot_product(gradient,gradient)) exit
            step = step*0.5_dp
            if (step < 1.0e-14_dp) exit
         end do
         x = trial
         if (objective < tol*tol) then
            if (socp_feasible(x,a,b,c,d,cone_sizes,strict=.true.)) info = 0
            return
         end if
      end do
   end subroutine phase1_socp

   function phase1_objective(x,a,b,c,d,cone_sizes) result(value)
      real(dp), intent(in) :: x(:),a(:,:),b(:),c(:,:),d(:)
      integer, intent(in) :: cone_sizes(:)
      real(dp) :: value,v,violation
      integer :: cone,row0,row1
      value = 0.0_dp
      row0 = 1
      do cone = 1, size(cone_sizes)
         row1 = row0+cone_sizes(cone)-1
         v = dot_product(c(cone,:),x)+d(cone)
         violation = vector_norm(matmul(a(row0:row1,:),x)+b(row0:row1))-v+1.0e-6_dp
         value = value+0.5_dp*max(violation,0.0_dp)**2
         row0 = row1+1
      end do
   end function phase1_objective

   subroutine barrier_value(x,t,f,a,b,c,d,cone_sizes,value,info)
      real(dp), intent(in) :: x(:),t,f(:),a(:,:),b(:),c(:,:),d(:)
      integer, intent(in) :: cone_sizes(:)
      real(dp), intent(out) :: value
      integer, intent(out) :: info
      real(dp) :: v,s
      real(dp), allocatable :: u(:)
      integer :: cone,row0,row1

      value = t*dot_product(f,x)
      info = 0
      row0 = 1
      do cone = 1, size(cone_sizes)
         row1 = row0+cone_sizes(cone)-1
         u = matmul(a(row0:row1,:),x)+b(row0:row1)
         v = dot_product(c(cone,:),x)+d(cone)
         s = v*v-dot_product(u,u)
         if (s <= 0.0_dp .or. v <= 0.0_dp) then
            info = 1
            value = huge(1.0_dp)
            return
         end if
         value = value-log(s)
         row0 = row1+1
      end do
   end subroutine barrier_value

   subroutine barrier_derivatives(x,t,f,a,b,c,d,cone_sizes,g,h,value,info)
      real(dp), intent(in) :: x(:),t,f(:),a(:,:),b(:),c(:,:),d(:)
      integer, intent(in) :: cone_sizes(:)
      real(dp), intent(out) :: g(:),h(:,:),value
      integer, intent(out) :: info
      real(dp), allocatable :: u(:),gs(:),hs(:,:)
      real(dp) :: v,s
      integer :: cone,row0,row1,n

      n = size(x)
      allocate(gs(n),hs(n,n))
      g = t*f
      h = 0.0_dp
      value = t*dot_product(f,x)
      info = 0
      row0 = 1
      do cone = 1, size(cone_sizes)
         row1 = row0+cone_sizes(cone)-1
         u = matmul(a(row0:row1,:),x)+b(row0:row1)
         v = dot_product(c(cone,:),x)+d(cone)
         s = v*v-dot_product(u,u)
         if (s <= 0.0_dp .or. v <= 0.0_dp) then
            info = 1
            return
         end if
         gs = 2.0_dp*(v*c(cone,:)-matmul(transpose(a(row0:row1,:)),u))
         hs = 2.0_dp*(outer(c(cone,:),c(cone,:))- &
            matmul(transpose(a(row0:row1,:)),a(row0:row1,:)))
         g = g-gs/s
         h = h+outer(gs,gs)/(s*s)-hs/s
         value = value-log(s)
         row0 = row1+1
      end do
      h = h+1.0e-12_dp*identity(n)
   end subroutine barrier_derivatives

   function socp_max_violation(x,a,b,c,d,cone_sizes) result(value)
      real(dp), intent(in) :: x(:),a(:,:),b(:),c(:,:),d(:)
      integer, intent(in) :: cone_sizes(:)
      real(dp) :: value,v
      integer :: cone,row0,row1
      value = -huge(1.0_dp)
      row0 = 1
      do cone = 1, size(cone_sizes)
         row1 = row0+cone_sizes(cone)-1
         v = vector_norm(matmul(a(row0:row1,:),x)+b(row0:row1))- &
            (dot_product(c(cone,:),x)+d(cone))
         value = max(value,v)
         row0 = row1+1
      end do
      if (size(cone_sizes) == 0) value = 0.0_dp
   end function socp_max_violation

   function socp_feasible(x,a,b,c,d,cone_sizes,strict) result(feasible)
      real(dp), intent(in) :: x(:),a(:,:),b(:),c(:,:),d(:)
      integer, intent(in) :: cone_sizes(:)
      logical, intent(in), optional :: strict
      logical :: feasible,is_strict
      real(dp) :: rhs,lhs,margin
      integer :: cone,row0,row1

      is_strict = .false.
      if (present(strict)) is_strict = strict
      margin = 0.0_dp
      if (is_strict) margin = 1.0e-10_dp
      feasible = .true.
      row0 = 1
      do cone = 1, size(cone_sizes)
         row1 = row0+cone_sizes(cone)-1
         lhs = vector_norm(matmul(a(row0:row1,:),x)+b(row0:row1))
         rhs = dot_product(c(cone,:),x)+d(cone)
         if (lhs+margin > rhs) then
            feasible = .false.
            return
         end if
         row0 = row1+1
      end do
   end function socp_feasible

   function outer(x,y) result(a)
      real(dp), intent(in) :: x(:),y(:)
      real(dp) :: a(size(x),size(y))
      a = spread(x,2,size(y))*spread(y,1,size(x))
   end function outer

   function identity(n) result(a)
      integer, intent(in) :: n
      real(dp) :: a(n,n)
      integer :: i
      a = 0.0_dp
      do i = 1,n
         a(i,i) = 1.0_dp
      end do
   end function identity

end module parma_socp
