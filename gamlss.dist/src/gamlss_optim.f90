! Copyright (C) 1998-2025 T. W. Yee, University of Auckland.
! Modern Fortran computational translation, 2026.
! SPDX-License-Identifier: GPL-3.0-only
module gamlss_optim
   use gamlss_kinds, only : dp
   implicit none
   private

   abstract interface
      function objective_function(x) result(f)
         import :: dp
         real(dp), intent(in) :: x(:)
         real(dp) :: f
      end function objective_function
   end interface

   public :: bfgs_minimize, numerical_gradient, numerical_hessian

contains

   subroutine numerical_gradient(fun, x, grad)
      procedure(objective_function) :: fun
      real(dp), intent(in) :: x(:)
      real(dp), intent(out) :: grad(:)
      real(dp), allocatable :: xp(:), xm(:)
      real(dp) :: h
      integer :: j

      allocate(xp(size(x)), xm(size(x)))
      do j = 1, size(x)
         h = epsilon(1.0_dp)**(1.0_dp/3.0_dp) * max(1.0_dp, abs(x(j)))
         xp = x
         xm = x
         xp(j) = x(j) + h
         xm(j) = x(j) - h
         grad(j) = (fun(xp) - fun(xm)) / (2.0_dp*h)
      end do
   end subroutine numerical_gradient

   subroutine numerical_hessian(fun, x, hess)
      procedure(objective_function) :: fun
      real(dp), intent(in) :: x(:)
      real(dp), intent(out) :: hess(:, :)
      real(dp), allocatable :: xp(:), xm(:), xpp(:), xpm(:), xmp(:), xmm(:)
      real(dp) :: hi, hj, f0
      integer :: i, j, n

      n = size(x)
      if (size(hess, 1) /= n .or. size(hess, 2) /= n) error stop "numerical_hessian: shape mismatch"
      allocate(xp(n), xm(n), xpp(n), xpm(n), xmp(n), xmm(n))
      f0 = fun(x)
      hess = 0.0_dp
      do i = 1, n
         hi = epsilon(1.0_dp)**0.25_dp * max(1.0_dp, abs(x(i)))
         xp = x
         xm = x
         xp(i) = x(i) + hi
         xm(i) = x(i) - hi
         hess(i, i) = (fun(xp) - 2.0_dp*f0 + fun(xm)) / (hi*hi)
         do j = i + 1, n
            hj = epsilon(1.0_dp)**0.25_dp * max(1.0_dp, abs(x(j)))
            xpp = x
            xpm = x
            xmp = x
            xmm = x
            xpp(i) = x(i) + hi
            xpp(j) = x(j) + hj
            xpm(i) = x(i) + hi
            xpm(j) = x(j) - hj
            xmp(i) = x(i) - hi
            xmp(j) = x(j) + hj
            xmm(i) = x(i) - hi
            xmm(j) = x(j) - hj
            hess(i, j) = (fun(xpp) - fun(xpm) - fun(xmp) + fun(xmm)) / (4.0_dp*hi*hj)
            hess(j, i) = hess(i, j)
         end do
      end do
   end subroutine numerical_hessian

   subroutine bfgs_minimize(fun, x, fval, status, max_iter, tol, inverse_hessian)
      procedure(objective_function) :: fun
      real(dp), intent(inout) :: x(:)
      real(dp), intent(out) :: fval
      integer, intent(out) :: status
      integer, intent(in), optional :: max_iter
      real(dp), intent(in), optional :: tol
      real(dp), intent(out), optional :: inverse_hessian(:, :)
      real(dp), allocatable :: hmat(:, :), ident(:, :), grad(:), grad_new(:)
      real(dp), allocatable :: direction(:), xnew(:), s(:), y(:), a(:, :)
      real(dp) :: tolerance, step, fnew, slope, rho, ys, gnorm, xnorm
      integer :: i, iter, n, niter

      n = size(x)
      if (n == 0) then
         fval = fun(x)
         status = 0
         return
      end if
      niter = 200
      if (present(max_iter)) niter = max_iter
      tolerance = 1.0e-7_dp
      if (present(tol)) tolerance = tol
      allocate(hmat(n,n), ident(n,n), grad(n), grad_new(n), direction(n), xnew(n), s(n), y(n), a(n,n))
      ident = 0.0_dp
      do i = 1, n
         ident(i, i) = 1.0_dp
      end do
      hmat = ident
      fval = fun(x)
      if (.not. is_finite(fval)) then
         status = 2
         if (present(inverse_hessian)) inverse_hessian = hmat
         return
      end if
      call numerical_gradient(fun, x, grad)
      status = 1

      do iter = 1, niter
         gnorm = maxval(abs(grad))
         xnorm = max(1.0_dp, maxval(abs(x)))
         if (gnorm <= tolerance*(1.0_dp + abs(fval))/xnorm) then
            status = 0
            exit
         end if

         direction = -matmul(hmat, grad)
         slope = dot_product(grad, direction)
         if (slope >= -tiny(1.0_dp) .or. .not. is_finite(slope)) then
            hmat = ident
            direction = -grad
            slope = -dot_product(grad, grad)
         end if

         step = 1.0_dp
         do
            xnew = x + step*direction
            fnew = fun(xnew)
            if (is_finite(fnew)) then
               if (fnew <= fval + 1.0e-4_dp*step*slope) exit
            end if
            step = 0.5_dp*step
            if (step < 1.0e-12_dp) exit
         end do
         if (step < 1.0e-12_dp) then
            status = 3
            exit
         end if

         s = xnew - x
         call numerical_gradient(fun, xnew, grad_new)
         y = grad_new - grad
         ys = dot_product(y, s)
         if (ys > sqrt(epsilon(1.0_dp))*sqrt(dot_product(y,y)*dot_product(s,s))) then
            rho = 1.0_dp/ys
            a = ident - rho*outer_product(s, y)
            hmat = matmul(a, matmul(hmat, transpose(a))) + rho*outer_product(s, s)
         else
            hmat = ident
         end if

         x = xnew
         fval = fnew
         grad = grad_new
         if (maxval(abs(s)) <= tolerance*max(1.0_dp, maxval(abs(x)))) then
            status = 0
            exit
         end if
      end do
      if (present(inverse_hessian)) inverse_hessian = hmat
   end subroutine bfgs_minimize

   pure function outer_product(a, b) result(c)
      real(dp), intent(in) :: a(:), b(:)
      real(dp) :: c(size(a), size(b))
      integer :: j
      do j = 1, size(b)
         c(:, j) = a*b(j)
      end do
   end function outer_product

   elemental logical function is_finite(x) result(ok)
      real(dp), intent(in) :: x
      ok = abs(x) <= huge(x)
   end function is_finite

end module gamlss_optim
