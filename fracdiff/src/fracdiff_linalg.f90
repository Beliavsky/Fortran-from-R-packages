! SPDX-License-Identifier: GPL-2.0-or-later
module fracdiff_linalg
   use fracdiff_kinds, only : dp
   implicit none
   private

   public :: solve_linear_system, invert_matrix, vector_norm2
   public :: symmetric_rank_k, matrix_is_finite

contains

   pure function vector_norm2(x) result(value)
      real(dp), intent(in) :: x(:)
      real(dp) :: value
      real(dp) :: scale, ssq, ax
      integer :: i

      scale = 0.0_dp
      ssq = 1.0_dp
      do i = 1, size(x)
         if (abs(x(i)) > 0.0_dp) then
            ax = abs(x(i))
            if (scale < ax) then
               ssq = 1.0_dp + ssq*(scale/ax)**2
               scale = ax
            else
               ssq = ssq + (ax/scale)**2
            end if
         end if
      end do
      if (scale <= 0.0_dp) then
         value = 0.0_dp
      else
         value = scale*sqrt(ssq)
      end if
   end function vector_norm2

   subroutine solve_linear_system(a, b, x, status)
      real(dp), intent(in) :: a(:,:)
      real(dp), intent(in) :: b(:)
      real(dp), intent(out) :: x(:)
      integer, intent(out) :: status

      real(dp), allocatable :: aa(:,:), bb(:)
      real(dp) :: factor, pivot_abs, temp
      integer :: n, i, j, k, pivot

      n = size(b)
      status = 0
      x = 0.0_dp
      if (size(a,1) /= n .or. size(a,2) /= n .or. size(x) /= n) then
         status = 1
         return
      end if
      if (n == 0) return

      allocate(aa(n,n), bb(n))
      aa = a
      bb = b

      do k = 1, n - 1
         pivot = k
         pivot_abs = abs(aa(k,k))
         do i = k + 1, n
            if (abs(aa(i,k)) > pivot_abs) then
               pivot = i
               pivot_abs = abs(aa(i,k))
            end if
         end do
         if (pivot_abs <= epsilon(1.0_dp)*max(1.0_dp, maxval(abs(aa)))) then
            status = 2
            return
         end if
         if (pivot /= k) then
            do j = k, n
               temp = aa(k,j)
               aa(k,j) = aa(pivot,j)
               aa(pivot,j) = temp
            end do
            temp = bb(k)
            bb(k) = bb(pivot)
            bb(pivot) = temp
         end if
         do i = k + 1, n
            factor = aa(i,k)/aa(k,k)
            aa(i,k) = 0.0_dp
            aa(i,k+1:n) = aa(i,k+1:n) - factor*aa(k,k+1:n)
            bb(i) = bb(i) - factor*bb(k)
         end do
      end do

      if (abs(aa(n,n)) <= epsilon(1.0_dp)*max(1.0_dp, maxval(abs(aa)))) then
         status = 2
         return
      end if

      x(n) = bb(n)/aa(n,n)
      do i = n - 1, 1, -1
         x(i) = (bb(i) - dot_product(aa(i,i+1:n), x(i+1:n)))/aa(i,i)
      end do
   end subroutine solve_linear_system

   subroutine invert_matrix(a, inverse, status)
      real(dp), intent(in) :: a(:,:)
      real(dp), intent(out) :: inverse(:,:)
      integer, intent(out) :: status

      real(dp), allocatable :: e(:), x(:)
      integer :: n, j, solve_status

      n = size(a,1)
      status = 0
      inverse = 0.0_dp
      if (size(a,2) /= n .or. size(inverse,1) /= n .or. size(inverse,2) /= n) then
         status = 1
         return
      end if
      allocate(e(n), x(n))
      do j = 1, n
         e = 0.0_dp
         e(j) = 1.0_dp
         call solve_linear_system(a, e, x, solve_status)
         if (solve_status /= 0) then
            status = solve_status
            inverse = 0.0_dp
            return
         end if
         inverse(:,j) = x
      end do
      inverse = 0.5_dp*(inverse + transpose(inverse))
   end subroutine invert_matrix

   pure subroutine symmetric_rank_k(jacobian, result)
      real(dp), intent(in) :: jacobian(:,:)
      real(dp), intent(out) :: result(:,:)
      integer :: i, j

      result = 0.0_dp
      do j = 1, size(jacobian,2)
         do i = 1, j
            result(i,j) = dot_product(jacobian(:,i), jacobian(:,j))
            result(j,i) = result(i,j)
         end do
      end do
   end subroutine symmetric_rank_k

   pure function matrix_is_finite(a) result(ok)
      use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
      real(dp), intent(in) :: a(:,:)
      logical :: ok
      integer :: i, j

      ok = .true.
      do j = 1, size(a,2)
         do i = 1, size(a,1)
            if (.not. ieee_is_finite(a(i,j))) then
               ok = .false.
               return
            end if
         end do
      end do
   end function matrix_is_finite

end module fracdiff_linalg
