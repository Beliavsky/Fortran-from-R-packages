! SPDX-License-Identifier: GPL-3.0-only
module smoots_linalg
   use smoots_kinds, only : dp
   use smoots_status, only : sm_ok, sm_invalid_input, sm_singular
   implicit none
   private
   public :: solve_linear_system, invert_matrix, least_squares_normal
contains
   subroutine solve_linear_system(a, b, x, status)
      real(dp), intent(in) :: a(:,:), b(:)
      real(dp), intent(out) :: x(:)
      integer, intent(out) :: status
      real(dp), allocatable :: aa(:,:), bb(:)
      real(dp) :: pivot_abs, factor, temp, scale
      integer :: n, i, j, k, pivot

      n = size(b)
      x = 0.0_dp
      status = sm_ok
      if (size(a,1) /= n .or. size(a,2) /= n .or. size(x) /= n) then
         status = sm_invalid_input
         return
      end if
      if (n == 0) return
      allocate(aa(n,n), bb(n))
      aa = a
      bb = b
      scale = max(1.0_dp, maxval(abs(aa)))
      do k = 1, n - 1
         pivot = k
         pivot_abs = abs(aa(k,k))
         do i = k + 1, n
            if (abs(aa(i,k)) > pivot_abs) then
               pivot = i
               pivot_abs = abs(aa(i,k))
            end if
         end do
         if (pivot_abs <= 100.0_dp*epsilon(1.0_dp)*scale) then
            status = sm_singular
            return
         end if
         if (pivot /= k) then
            do j = k, n
               temp = aa(k,j); aa(k,j) = aa(pivot,j); aa(pivot,j) = temp
            end do
            temp = bb(k); bb(k) = bb(pivot); bb(pivot) = temp
         end if
         do i = k + 1, n
            factor = aa(i,k)/aa(k,k)
            aa(i,k) = 0.0_dp
            aa(i,k+1:n) = aa(i,k+1:n) - factor*aa(k,k+1:n)
            bb(i) = bb(i) - factor*bb(k)
         end do
      end do
      if (abs(aa(n,n)) <= 100.0_dp*epsilon(1.0_dp)*scale) then
         status = sm_singular
         return
      end if
      x(n) = bb(n)/aa(n,n)
      do i = n - 1, 1, -1
         x(i) = (bb(i) - dot_product(aa(i,i+1:n), x(i+1:n)))/aa(i,i)
      end do
   end subroutine solve_linear_system

   subroutine invert_matrix(a, ainv, status)
      real(dp), intent(in) :: a(:,:)
      real(dp), intent(out) :: ainv(:,:)
      integer, intent(out) :: status
      real(dp), allocatable :: e(:), x(:)
      integer :: n, j, istat
      n = size(a,1)
      ainv = 0.0_dp
      if (size(a,2) /= n .or. size(ainv,1) /= n .or. size(ainv,2) /= n) then
         status = sm_invalid_input
         return
      end if
      allocate(e(n), x(n))
      do j = 1, n
         e = 0.0_dp
         e(j) = 1.0_dp
         call solve_linear_system(a, e, x, istat)
         if (istat /= sm_ok) then
            status = istat
            return
         end if
         ainv(:,j) = x
      end do
      status = sm_ok
   end subroutine invert_matrix

   subroutine least_squares_normal(x, y, beta, status, ridge)
      real(dp), intent(in) :: x(:,:), y(:)
      real(dp), intent(out) :: beta(:)
      integer, intent(out) :: status
      real(dp), intent(in), optional :: ridge
      real(dp), allocatable :: xtx(:,:), xty(:)
      real(dp) :: rr
      integer :: j, p
      p = size(x,2)
      if (size(x,1) /= size(y) .or. size(beta) /= p) then
         status = sm_invalid_input
         beta = 0.0_dp
         return
      end if
      allocate(xtx(p,p), xty(p))
      xtx = matmul(transpose(x), x)
      xty = matmul(transpose(x), y)
      rr = 0.0_dp
      if (present(ridge)) rr = max(0.0_dp, ridge)
      do j = 1, p
         xtx(j,j) = xtx(j,j) + rr
      end do
      call solve_linear_system(xtx, xty, beta, status)
   end subroutine least_squares_normal
end module smoots_linalg
