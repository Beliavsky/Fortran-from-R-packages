! SPDX-License-Identifier: GPL-3.0-only
module yrnd_optimize
   use yrnd_kinds, only : dp
   implicit none
   private

   abstract interface
      function objective_function(x) result(value)
         import dp
         real(dp), intent(in) :: x(:)
         real(dp) :: value
      end function objective_function
   end interface

   public :: nelder_mead_bounded

contains

   subroutine nelder_mead_bounded(fun, start, lower, upper, solution, value, status, &
                                  max_iter, tolerance)
      procedure(objective_function) :: fun
      real(dp), intent(in) :: start(:), lower(:), upper(:)
      real(dp), intent(out) :: solution(size(start)), value
      integer, intent(out) :: status
      integer, intent(in), optional :: max_iter
      real(dp), intent(in), optional :: tolerance

      real(dp), parameter :: alpha = 1.0_dp
      real(dp), parameter :: gamma = 2.0_dp
      real(dp), parameter :: rho = 0.5_dp
      real(dp), parameter :: sigma = 0.5_dp
      integer :: n, niter, iter, i, j, ilo, ihi, inhi
      real(dp) :: tol, spread, fspread
      real(dp), allocatable :: simplex(:, :), f(:), centroid(:), xr(:), xe(:), xc(:)

      n = size(start)
      if (size(lower) /= n .or. size(upper) /= n) error stop "nelder_mead_bounded: size mismatch"
      niter = 1000
      if (present(max_iter)) niter = max_iter
      tol = 1.0e-8_dp
      if (present(tolerance)) tol = tolerance

      allocate(simplex(n, n + 1), f(n + 1), centroid(n), xr(n), xe(n), xc(n))
      simplex(:, 1) = clamp(start, lower, upper)
      do j = 2, n + 1
         simplex(:, j) = simplex(:, 1)
         i = j - 1
         if (upper(i) > lower(i)) then
            simplex(i, j) = simplex(i, j) + 0.05_dp * (upper(i) - lower(i))
            if (simplex(i, j) > upper(i)) then
               simplex(i, j) = simplex(i, 1) - 0.05_dp * (upper(i) - lower(i))
            end if
         end if
         simplex(:, j) = clamp(simplex(:, j), lower, upper)
      end do
      do j = 1, n + 1
         f(j) = fun(simplex(:, j))
      end do

      status = 1
      do iter = 1, niter
         call order_indices(f, ilo, ihi, inhi)
         spread = 0.0_dp
         do j = 1, n + 1
            spread = max(spread, maxval(abs(simplex(:, j) - simplex(:, ilo))))
         end do
         fspread = maxval(abs(f - f(ilo)))
         if (spread <= tol * (1.0_dp + maxval(abs(simplex(:, ilo)))) .and. &
             fspread <= tol * (1.0_dp + abs(f(ilo)))) then
            status = 0
            exit
         end if

         centroid = 0.0_dp
         do j = 1, n + 1
            if (j /= ihi) centroid = centroid + simplex(:, j)
         end do
         centroid = centroid / real(n, dp)

         xr = clamp(centroid + alpha * (centroid - simplex(:, ihi)), lower, upper)
         if (fun(xr) < f(ilo)) then
            xe = clamp(centroid + gamma * (xr - centroid), lower, upper)
            if (fun(xe) < fun(xr)) then
               simplex(:, ihi) = xe
               f(ihi) = fun(xe)
            else
               simplex(:, ihi) = xr
               f(ihi) = fun(xr)
            end if
         else if (fun(xr) < f(inhi)) then
            simplex(:, ihi) = xr
            f(ihi) = fun(xr)
         else
            if (fun(xr) < f(ihi)) then
               xc = clamp(centroid + rho * (xr - centroid), lower, upper)
            else
               xc = clamp(centroid + rho * (simplex(:, ihi) - centroid), lower, upper)
            end if
            if (fun(xc) < min(f(ihi), fun(xr))) then
               simplex(:, ihi) = xc
               f(ihi) = fun(xc)
            else
               do j = 1, n + 1
                  if (j /= ilo) then
                     simplex(:, j) = clamp(simplex(:, ilo) + &
                        sigma * (simplex(:, j) - simplex(:, ilo)), lower, upper)
                     f(j) = fun(simplex(:, j))
                  end if
               end do
            end if
         end if
      end do

      call order_indices(f, ilo, ihi, inhi)
      solution = simplex(:, ilo)
      value = f(ilo)

   contains

      pure function clamp(x, lo, hi) result(y)
         real(dp), intent(in) :: x(:), lo(:), hi(:)
         real(dp) :: y(size(x))
         y = min(max(x, lo), hi)
      end function clamp

      subroutine order_indices(values, low, high, next_high)
         real(dp), intent(in) :: values(:)
         integer, intent(out) :: low, high, next_high
         integer :: k
         low = minloc(values, dim=1)
         high = maxloc(values, dim=1)
         next_high = low
         do k = 1, size(values)
            if (k == high) cycle
            if (next_high == high .or. values(k) > values(next_high)) next_high = k
         end do
      end subroutine order_indices

   end subroutine nelder_mead_bounded

end module yrnd_optimize
