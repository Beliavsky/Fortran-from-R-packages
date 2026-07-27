! Part of the experimental modern Fortran translation of rugarch 1.5-6.
! This file was created or modified for the Fortran project on 2026-07-23.
! Original rugarch authors retain copyright; see NOTICE and ORIGIN.md.
! SPDX-License-Identifier: GPL-3.0-only

module rugarch_optimizer
   use rugarch_kinds, only : dp
   implicit none
   private

   type, public :: optimizer_result
      real(dp), allocatable :: x(:)
      real(dp) :: objective = huge(1.0_dp)
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

   function nelder_mead(fun, x0, context, step, tolerance, max_iterations) result(result)
      procedure(objective_function) :: fun
      real(dp), intent(in) :: x0(:)
      class(*), intent(in) :: context
      real(dp), intent(in), optional :: step, tolerance
      integer, intent(in), optional :: max_iterations
      type(optimizer_result) :: result
      real(dp), allocatable :: simplex(:,:), f(:), centroid(:), xr(:), xe(:), xc(:)
      real(dp) :: initial_step, tol, fr, fe, fc, spread_x, spread_f
      real(dp), parameter :: alpha = 1.0_dp
      real(dp), parameter :: gamma = 2.0_dp
      real(dp), parameter :: rho = 0.5_dp
      real(dp), parameter :: sigma = 0.5_dp
      integer :: n, maxit, i, j, iter

      n = size(x0)
      initial_step = 0.1_dp
      if (present(step)) initial_step = step
      tol = 1.0e-7_dp
      if (present(tolerance)) tol = tolerance
      maxit = 2000
      if (present(max_iterations)) maxit = max_iterations

      allocate(simplex(n,n+1),f(n+1),centroid(n),xr(n),xe(n),xc(n))
      simplex(:,1) = x0
      do j = 2, n+1
         simplex(:,j) = x0
         i = j-1
         simplex(i,j) = simplex(i,j) + initial_step*max(1.0_dp,abs(x0(i)))
      end do
      do j = 1, n+1
         f(j) = fun(simplex(:,j),context)
      end do
      result%evaluations = n+1

      do iter = 1, maxit
         call sort_simplex(simplex,f)
         spread_f = maxval(abs(f-f(1)))
         spread_x = maxval(abs(simplex-spread(simplex(:,1),dim=2,ncopies=n+1)))
         if (spread_f <= tol*(1.0_dp+abs(f(1))) .and. &
             spread_x <= sqrt(tol)*(1.0_dp+maxval(abs(simplex(:,1))))) then
            result%status = 0
            exit
         end if

         centroid = sum(simplex(:,1:n),dim=2)/real(n,dp)
         xr = centroid + alpha*(centroid-simplex(:,n+1))
         fr = fun(xr,context)
         result%evaluations = result%evaluations+1

         if (fr < f(1)) then
            xe = centroid + gamma*(xr-centroid)
            fe = fun(xe,context)
            result%evaluations = result%evaluations+1
            if (fe < fr) then
               simplex(:,n+1) = xe
               f(n+1) = fe
            else
               simplex(:,n+1) = xr
               f(n+1) = fr
            end if
         else if (fr < f(n)) then
            simplex(:,n+1) = xr
            f(n+1) = fr
         else
            if (fr < f(n+1)) then
               xc = centroid + rho*(xr-centroid)
            else
               xc = centroid + rho*(simplex(:,n+1)-centroid)
            end if
            fc = fun(xc,context)
            result%evaluations = result%evaluations+1
            if (fc < min(fr,f(n+1))) then
               simplex(:,n+1) = xc
               f(n+1) = fc
            else
               do j = 2, n+1
                  simplex(:,j) = simplex(:,1)+sigma*(simplex(:,j)-simplex(:,1))
                  f(j) = fun(simplex(:,j),context)
               end do
               result%evaluations = result%evaluations+n
            end if
         end if
      end do

      call sort_simplex(simplex,f)
      allocate(result%x(n))
      result%x = simplex(:,1)
      result%objective = f(1)
      result%iterations = min(iter,maxit)
   end function nelder_mead

   subroutine sort_simplex(simplex, f)
      real(dp), intent(inout) :: simplex(:,:), f(:)
      real(dp), allocatable :: tmpx(:)
      real(dp) :: tmpf
      integer :: i, j, k, n

      n = size(f)
      allocate(tmpx(size(simplex,1)))
      do i = 1, n-1
         k = i
         do j = i+1, n
            if (f(j) < f(k)) k = j
         end do
         if (k /= i) then
            tmpf = f(i); f(i) = f(k); f(k) = tmpf
            tmpx = simplex(:,i)
            simplex(:,i) = simplex(:,k)
            simplex(:,k) = tmpx
         end if
      end do
   end subroutine sort_simplex

end module rugarch_optimizer
