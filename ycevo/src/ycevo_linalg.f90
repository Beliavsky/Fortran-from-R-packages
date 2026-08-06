module ycevo_linalg
   use ycevo_kinds, only : dp
   use ycevo_status, only : ycevo_success, ycevo_err_singular
   implicit none
   private

   public :: solve_linear_system, weighted_quadratic_fit

contains

   subroutine solve_linear_system(a, b, x, status)
      real(dp), intent(in) :: a(:, :), b(:)
      real(dp), allocatable, intent(out) :: x(:)
      integer, intent(out) :: status
      real(dp), allocatable :: aa(:, :), bb(:), row(:)
      real(dp) :: factor, pivot_scale
      integer :: i, k, n, pivot

      n = size(b)
      status = ycevo_err_singular
      allocate(x(n))
      x = 0.0_dp
      if (size(a, 1) /= n .or. size(a, 2) /= n) return
      allocate(aa(n, n), bb(n), row(n))
      aa = a
      bb = b

      do k = 1, n - 1
         pivot = k - 1 + maxloc(abs(aa(k:n, k)), dim=1)
         pivot_scale = max(1.0_dp, maxval(abs(aa(k:n, k:n))))
         if (abs(aa(pivot, k)) <= 100.0_dp*epsilon(1.0_dp)*pivot_scale) return
         if (pivot /= k) then
            row = aa(k, :)
            aa(k, :) = aa(pivot, :)
            aa(pivot, :) = row
            factor = bb(k)
            bb(k) = bb(pivot)
            bb(pivot) = factor
         end if
         do i = k + 1, n
            factor = aa(i, k) / aa(k, k)
            aa(i, k) = 0.0_dp
            aa(i, k+1:n) = aa(i, k+1:n) - factor*aa(k, k+1:n)
            bb(i) = bb(i) - factor*bb(k)
         end do
      end do
      if (abs(aa(n, n)) <= 100.0_dp*epsilon(1.0_dp)*max(1.0_dp, abs(aa(n,n)))) return

      x(n) = bb(n) / aa(n, n)
      do i = n - 1, 1, -1
         x(i) = (bb(i) - dot_product(aa(i, i+1:n), x(i+1:n))) / aa(i, i)
      end do
      status = ycevo_success
   end subroutine solve_linear_system

   subroutine weighted_quadratic_fit(x, y, w, coef, status)
      real(dp), intent(in) :: x(:), y(:), w(:)
      real(dp), intent(out) :: coef(3)
      integer, intent(out) :: status
      real(dp) :: a(3,3), b(3), xp(0:4)
      real(dp), allocatable :: sol(:)
      integer :: i, j, k

      coef = 0.0_dp
      status = ycevo_err_singular
      if (size(y) /= size(x) .or. size(w) /= size(x)) return
      a = 0.0_dp
      b = 0.0_dp
      do i = 1, size(x)
         xp(0) = 1.0_dp
         do j = 1, 4
            xp(j) = xp(j-1)*x(i)
         end do
         do j = 1, 3
            b(j) = b(j) + w(i)*y(i)*xp(j-1)
            do k = 1, 3
               a(j,k) = a(j,k) + w(i)*xp(j+k-2)
            end do
         end do
      end do
      call solve_linear_system(a, b, sol, status)
      if (status == ycevo_success) coef = sol
   end subroutine weighted_quadratic_fit

end module ycevo_linalg
