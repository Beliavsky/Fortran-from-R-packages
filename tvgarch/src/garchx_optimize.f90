! SPDX-License-Identifier: GPL-2.0-or-later
! This file is part of a modern Fortran translation of garchx.
! Copyright (C) 2026 translation contributors.
! Original garchx package copyright (C) Genaro Sucarrat.
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 2 of the License, or
! (at your option) any later version.
module garchx_optimize
   use garchx_kinds, only : dp
   implicit none
   private
   public :: bounded_nelder_mead, numerical_hessian, objective_function

   abstract interface
      function objective_function(x) result(value)
         import dp
         real(dp), intent(in) :: x(:)
         real(dp) :: value
      end function objective_function
   end interface
contains
   subroutine bounded_nelder_mead(fun, x0, lower, upper, xbest, fbest, status, &
                                  iterations, max_iter, rel_tol)
      procedure(objective_function) :: fun
      real(dp), intent(in) :: x0(:), lower(:), upper(:)
      real(dp), allocatable, intent(out) :: xbest(:)
      real(dp), intent(out) :: fbest
      integer, intent(out) :: status, iterations
      integer, intent(in), optional :: max_iter
      real(dp), intent(in), optional :: rel_tol
      real(dp), parameter :: alpha = 1.0_dp, gamma = 2.0_dp
      real(dp), parameter :: rho = 0.5_dp, sigma = 0.5_dp
      integer :: n, maxit, j, worst, second_worst, best
      real(dp) :: tol, fspread, xspread, step
      real(dp), allocatable :: simplex(:, :), f(:), centroid(:), xr(:), xe(:), xc(:)

      n = size(x0)
      if (n < 1 .or. size(lower) /= n .or. size(upper) /= n .or. any(lower > upper)) then
         status = 2
         iterations = 0
         fbest = huge(1.0_dp)
         allocate(xbest(0))
         return
      end if
      maxit = 3000
      if (present(max_iter)) maxit = max_iter
      tol = 1.0e-8_dp
      if (present(rel_tol)) tol = rel_tol
      allocate(simplex(n, n+1), f(n+1), centroid(n), xr(n), xe(n), xc(n), xbest(n))
      simplex(:, 1) = min(max(x0, lower), upper)
      do j = 1, n
         simplex(:, j+1) = simplex(:, 1)
         step = 0.05_dp*max(1.0_dp, abs(simplex(j, 1)))
         if (simplex(j, j+1)+step <= upper(j)) then
            simplex(j, j+1) = simplex(j, j+1)+step
         else if (simplex(j, j+1)-step >= lower(j)) then
            simplex(j, j+1) = simplex(j, j+1)-step
         else
            simplex(j, j+1) = 0.5_dp*(lower(j)+upper(j))
         end if
      end do
      do j = 1, n+1
         f(j) = fun(simplex(:, j))
      end do

      status = 1
      do iterations = 1, maxit
         best = minloc(f, dim=1)
         worst = maxloc(f, dim=1)
         second_worst = best
         do j = 1, n+1
            if (j /= worst) then
               if (second_worst == best .or. f(j) > f(second_worst)) second_worst = j
            end if
         end do
         fspread = maxval(abs(f-f(best)))/max(1.0_dp, abs(f(best)))
         xspread = 0.0_dp
         do j = 1, n+1
            xspread = max(xspread, maxval(abs(simplex(:, j)-simplex(:, best))/ &
                      max(1.0_dp, abs(simplex(:, best)))))
         end do
         if (fspread <= tol .and. xspread <= sqrt(tol)) then
            status = 0
            exit
         end if

         centroid = 0.0_dp
         do j = 1, n+1
            if (j /= worst) centroid = centroid+simplex(:, j)
         end do
         centroid = centroid/real(n, dp)
         xr = min(max(centroid+alpha*(centroid-simplex(:, worst)), lower), upper)
         if (fun(xr) < f(best)) then
            xe = min(max(centroid+gamma*(xr-centroid), lower), upper)
            if (fun(xe) < fun(xr)) then
               simplex(:, worst) = xe
               f(worst) = fun(xe)
            else
               simplex(:, worst) = xr
               f(worst) = fun(xr)
            end if
         else if (fun(xr) < f(second_worst)) then
            simplex(:, worst) = xr
            f(worst) = fun(xr)
         else
            if (fun(xr) < f(worst)) then
               xc = min(max(centroid+rho*(xr-centroid), lower), upper)
            else
               xc = min(max(centroid-rho*(centroid-simplex(:, worst)), lower), upper)
            end if
            if (fun(xc) < min(fun(xr), f(worst))) then
               simplex(:, worst) = xc
               f(worst) = fun(xc)
            else
               do j = 1, n+1
                  if (j /= best) then
                     simplex(:, j) = min(max(simplex(:, best)+ &
                                      sigma*(simplex(:, j)-simplex(:, best)), lower), upper)
                     f(j) = fun(simplex(:, j))
                  end if
               end do
            end if
         end if
      end do
      if (iterations > maxit) iterations = maxit
      best = minloc(f, dim=1)
      xbest = simplex(:, best)
      fbest = f(best)
   end subroutine bounded_nelder_mead

   subroutine numerical_hessian(fun, x, hessian)
      procedure(objective_function) :: fun
      real(dp), intent(in) :: x(:)
      real(dp), allocatable, intent(out) :: hessian(:, :)
      integer :: n, i, j
      real(dp) :: f0, hi, hj
      real(dp), allocatable :: xp(:), xm(:), xpp(:), xpm(:), xmp(:), xmm(:), h(:)

      n = size(x)
      allocate(hessian(n, n), xp(n), xm(n), xpp(n), xpm(n), xmp(n), xmm(n), h(n))
      h = epsilon(1.0_dp)**0.25_dp*max(1.0_dp, abs(x))
      f0 = fun(x)
      hessian = 0.0_dp
      do i = 1, n
         hi = h(i)
         xp = x
         xm = x
         xp(i) = xp(i)+hi
         xm(i) = xm(i)-hi
         hessian(i, i) = (fun(xp)-2.0_dp*f0+fun(xm))/(hi*hi)
         do j = i+1, n
            hj = h(j)
            xpp = x
            xpm = x
            xmp = x
            xmm = x
            xpp(i) = xpp(i)+hi; xpp(j) = xpp(j)+hj
            xpm(i) = xpm(i)+hi; xpm(j) = xpm(j)-hj
            xmp(i) = xmp(i)-hi; xmp(j) = xmp(j)+hj
            xmm(i) = xmm(i)-hi; xmm(j) = xmm(j)-hj
            hessian(i, j) = (fun(xpp)-fun(xpm)-fun(xmp)+fun(xmm))/(4.0_dp*hi*hj)
            hessian(j, i) = hessian(i, j)
         end do
      end do
   end subroutine numerical_hessian
end module garchx_optimize
