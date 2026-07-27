! SPDX-License-Identifier: GPL-2.0-or-later
!
! Computational translation of gogarch, copyright (C) 2008-2026 Bernhard Pfaff.
! Fortran translation copyright (C) 2026 translation contributors.
! Distributed under the GNU General Public License, version 2 or later.
module gogarch_optimizer
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   use gogarch_kinds, only : dp
   implicit none
   private

   type, public :: optimizer_result
      real(dp), allocatable :: x(:)
      real(dp) :: value = huge(1.0_dp)
      integer :: iterations = 0
      integer :: evaluations = 0
      integer :: status = 1
   end type optimizer_result

   abstract interface
      function objective_function(x, context) result(value)
         import :: dp
         real(dp), intent(in) :: x(:)
         class(*), intent(in) :: context
         real(dp) :: value
      end function objective_function
   end interface

   public :: nelder_mead

contains

   recursive function nelder_mead(objective, x0, context, step, tolerance, max_iterations) result(result)
      procedure(objective_function) :: objective
      real(dp), intent(in) :: x0(:)
      class(*), intent(in) :: context
      real(dp), intent(in), optional :: step, tolerance
      integer, intent(in), optional :: max_iterations
      type(optimizer_result) :: result
      real(dp), allocatable :: simplex(:,:), values(:), centroid(:), xr(:), xe(:), xc(:)
      real(dp) :: alpha, gamma, rho, sigma, initial_step, tol
      real(dp) :: fr, fe, fc, spread_x, spread_f
      integer :: n, maxit, iter, i, best, worst, second_worst

      n = size(x0)
      maxit = 500
      if (present(max_iterations)) maxit = max(1,max_iterations)
      initial_step = 0.15_dp
      if (present(step)) initial_step = max(step,1.0e-8_dp)
      tol = 1.0e-8_dp
      if (present(tolerance)) tol = max(tolerance,epsilon(1.0_dp))
      alpha = 1.0_dp
      gamma = 2.0_dp
      rho = 0.5_dp
      sigma = 0.5_dp

      allocate(simplex(n,n+1),values(n+1),centroid(n),xr(n),xe(n),xc(n),result%x(n))
      simplex(:,1) = x0
      do i = 1, n
         simplex(:,i+1) = x0
         simplex(i,i+1) = simplex(i,i+1)+initial_step*max(1.0_dp,abs(x0(i)))
      end do
      do i = 1, n+1
         values(i) = safe_objective(objective,simplex(:,i),context)
      end do
      result%evaluations = n+1

      do iter = 1, maxit
         call order_simplex(simplex,values)
         spread_f = maxval(abs(values-values(1)))
         spread_x = maxval(abs(simplex-spread(simplex(:,1),2,n+1)))
         if (spread_f <= tol*(1.0_dp+abs(values(1))) .and. spread_x <= sqrt(tol)*(1.0_dp+maxval(abs(simplex(:,1))))) then
            result%status = 0
            exit
         end if

         best = 1
         worst = n+1
         second_worst = n
         centroid = sum(simplex(:,1:n),dim=2)/real(n,dp)
         xr = centroid+alpha*(centroid-simplex(:,worst))
         fr = safe_objective(objective,xr,context)
         result%evaluations = result%evaluations+1

         if (fr < values(best)) then
            xe = centroid+gamma*(xr-centroid)
            fe = safe_objective(objective,xe,context)
            result%evaluations = result%evaluations+1
            if (fe < fr) then
               simplex(:,worst) = xe
               values(worst) = fe
            else
               simplex(:,worst) = xr
               values(worst) = fr
            end if
         else if (fr < values(second_worst)) then
            simplex(:,worst) = xr
            values(worst) = fr
         else
            if (fr < values(worst)) then
               xc = centroid+rho*(xr-centroid)
            else
               xc = centroid-rho*(centroid-simplex(:,worst))
            end if
            fc = safe_objective(objective,xc,context)
            result%evaluations = result%evaluations+1
            if (fc < min(fr,values(worst))) then
               simplex(:,worst) = xc
               values(worst) = fc
            else
               do i = 2, n+1
                  simplex(:,i) = simplex(:,best)+sigma*(simplex(:,i)-simplex(:,best))
                  values(i) = safe_objective(objective,simplex(:,i),context)
               end do
               result%evaluations = result%evaluations+n
            end if
         end if
      end do

      call order_simplex(simplex,values)
      result%x = simplex(:,1)
      result%value = values(1)
      result%iterations = min(iter,maxit)
   end function nelder_mead

   recursive function safe_objective(objective, x, context) result(value)
      procedure(objective_function) :: objective
      real(dp), intent(in) :: x(:)
      class(*), intent(in) :: context
      real(dp) :: value
      value = objective(x,context)
      if (.not. ieee_is_finite(value)) value = huge(1.0_dp)/100.0_dp
   end function safe_objective

   subroutine order_simplex(simplex, values)
      real(dp), intent(inout) :: simplex(:,:), values(:)
      real(dp) :: temp_value, temp_vector(size(simplex,1))
      integer :: i, j
      do i = 2, size(values)
         temp_value = values(i)
         temp_vector = simplex(:,i)
         j = i-1
         do while (j >= 1)
            if (values(j) <= temp_value) exit
            values(j+1) = values(j)
            simplex(:,j+1) = simplex(:,j)
            j = j-1
         end do
         values(j+1) = temp_value
         simplex(:,j+1) = temp_vector
      end do
   end subroutine order_simplex

end module gogarch_optimizer
