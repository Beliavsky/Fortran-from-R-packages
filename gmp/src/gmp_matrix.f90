module gmp_matrix
   use, intrinsic :: iso_fortran_env, only : int64
   use gmp_bigz, only : bigz, bigz_from_int64, bigz_is_zero, bigz_add, bigz_mul
   use gmp_bigq, only : bigq, bigq_make, bigq_from_int64, bigq_is_zero, bigq_add, bigq_sub, bigq_mul, bigq_div
   implicit none
   private

   public :: bigz_matmul, bigq_matmul
   public :: bigz_crossprod, bigz_tcrossprod, bigq_crossprod, bigq_tcrossprod
   public :: bigz_transpose, bigq_transpose
   public :: bigq_solve, bigq_inverse, bigz_solve, bigz_inverse

contains

   pure subroutine bigz_matmul(a, b, c, info)
      type(bigz), intent(in) :: a(:, :) !! Left integer matrix with shape m by k.
      type(bigz), intent(in) :: b(:, :) !! Right integer matrix with shape k by n.
      type(bigz), allocatable, intent(out) :: c(:, :) !! Exact integer matrix product with shape m by n.
      integer, intent(out), optional :: info !! Zero on success and one when inner dimensions do not conform.
      integer :: i, j, k
      type(bigz) :: acc

      if (size(a, 2) /= size(b, 1)) then
         allocate(c(0, 0))
         if (present(info)) info = 1
         return
      end if
      allocate(c(size(a, 1), size(b, 2)))
      do j = 1, size(b, 2)
         do i = 1, size(a, 1)
            acc = bigz_from_int64(0_int64)
            do k = 1, size(a, 2)
               acc = bigz_add(acc, bigz_mul(a(i, k), b(k, j)))
            end do
            c(i, j) = acc
         end do
      end do
      if (present(info)) info = 0
   end subroutine bigz_matmul

   pure subroutine bigq_matmul(a, b, c, info)
      type(bigq), intent(in) :: a(:, :) !! Left rational matrix with shape m by k.
      type(bigq), intent(in) :: b(:, :) !! Right rational matrix with shape k by n.
      type(bigq), allocatable, intent(out) :: c(:, :) !! Exact rational matrix product with shape m by n.
      integer, intent(out), optional :: info !! Zero on success and one when inner dimensions do not conform.
      integer :: i, j, k
      type(bigq) :: acc

      if (size(a, 2) /= size(b, 1)) then
         allocate(c(0, 0))
         if (present(info)) info = 1
         return
      end if
      allocate(c(size(a, 1), size(b, 2)))
      do j = 1, size(b, 2)
         do i = 1, size(a, 1)
            acc = bigq_from_int64(0_int64)
            do k = 1, size(a, 2)
               acc = bigq_add(acc, bigq_mul(a(i, k), b(k, j)))
            end do
            c(i, j) = acc
         end do
      end do
      if (present(info)) info = 0
   end subroutine bigq_matmul

   pure function bigz_transpose(a) result(t)
      type(bigz), intent(in) :: a(:, :) !! Integer matrix to transpose.
      type(bigz), allocatable :: t(:, :)
      integer :: i, j

      allocate(t(size(a, 2), size(a, 1)))
      do j = 1, size(a, 2)
         do i = 1, size(a, 1)
            t(j, i) = a(i, j)
         end do
      end do
   end function bigz_transpose

   pure function bigq_transpose(a) result(t)
      type(bigq), intent(in) :: a(:, :) !! Rational matrix to transpose.
      type(bigq), allocatable :: t(:, :)
      integer :: i, j

      allocate(t(size(a, 2), size(a, 1)))
      do j = 1, size(a, 2)
         do i = 1, size(a, 1)
            t(j, i) = a(i, j)
         end do
      end do
   end function bigq_transpose

   pure subroutine bigz_crossprod(x, result, info)
      type(bigz), intent(in) :: x(:, :) !! Integer matrix x in t(x) times x.
      type(bigz), allocatable, intent(out) :: result(:, :) !! Exact cross-product matrix with shape ncol(x) by ncol(x).
      integer, intent(out), optional :: info !! Zero on success; retained for a uniform matrix API.
      type(bigz), allocatable :: tx(:, :)
      integer :: stat

      tx = bigz_transpose(x)
      call bigz_matmul(tx, x, result, stat)
      if (present(info)) info = stat
   end subroutine bigz_crossprod

   pure subroutine bigz_tcrossprod(x, result, info)
      type(bigz), intent(in) :: x(:, :) !! Integer matrix x in x times t(x).
      type(bigz), allocatable, intent(out) :: result(:, :) !! Exact tcross-product matrix with shape nrow(x) by nrow(x).
      integer, intent(out), optional :: info !! Zero on success; retained for a uniform matrix API.
      type(bigz), allocatable :: tx(:, :)
      integer :: stat

      tx = bigz_transpose(x)
      call bigz_matmul(x, tx, result, stat)
      if (present(info)) info = stat
   end subroutine bigz_tcrossprod

   pure subroutine bigq_crossprod(x, result, info)
      type(bigq), intent(in) :: x(:, :) !! Rational matrix x in t(x) times x.
      type(bigq), allocatable, intent(out) :: result(:, :) !! Exact cross-product matrix with shape ncol(x) by ncol(x).
      integer, intent(out), optional :: info !! Zero on success; retained for a uniform matrix API.
      type(bigq), allocatable :: tx(:, :)
      integer :: stat

      tx = bigq_transpose(x)
      call bigq_matmul(tx, x, result, stat)
      if (present(info)) info = stat
   end subroutine bigq_crossprod

   pure subroutine bigq_tcrossprod(x, result, info)
      type(bigq), intent(in) :: x(:, :) !! Rational matrix x in x times t(x).
      type(bigq), allocatable, intent(out) :: result(:, :) !! Exact tcross-product matrix with shape nrow(x) by nrow(x).
      integer, intent(out), optional :: info !! Zero on success; retained for a uniform matrix API.
      type(bigq), allocatable :: tx(:, :)
      integer :: stat

      tx = bigq_transpose(x)
      call bigq_matmul(x, tx, result, stat)
      if (present(info)) info = stat
   end subroutine bigq_tcrossprod

   pure subroutine bigq_solve(a, b, x, info)
      type(bigq), intent(in) :: a(:, :) !! Square rational coefficient matrix with shape n by n.
      type(bigq), intent(in) :: b(:, :) !! Rational right-hand side matrix with shape n by nrhs.
      type(bigq), allocatable, intent(out) :: x(:, :) !! Exact rational solution matrix with shape n by nrhs.
      integer, intent(out), optional :: info !! Zero on success; one for dimension mismatch; two for a singular matrix.
      type(bigq), allocatable :: aa(:, :)
      type(bigq) :: pivot, factor, tmp
      integer :: n, nrhs, i, j, k, p, stat

      n = size(a, 1)
      nrhs = size(b, 2)
      stat = 0
      if (size(a, 2) /= n .or. size(b, 1) /= n) then
         allocate(x(0, 0))
         stat = 1
         if (present(info)) info = stat
         return
      end if
      aa = a
      x = b
      do k = 1, n
         p = 0
         do i = k, n
            if (.not. bigq_is_zero(aa(i, k))) then
               p = i
               exit
            end if
         end do
         if (p == 0) then
            stat = 2
            if (present(info)) info = stat
            return
         end if
         if (p /= k) then
            do j = 1, n
               tmp = aa(k, j)
               aa(k, j) = aa(p, j)
               aa(p, j) = tmp
            end do
            do j = 1, nrhs
               tmp = x(k, j)
               x(k, j) = x(p, j)
               x(p, j) = tmp
            end do
         end if

         pivot = aa(k, k)
         do j = k, n
            aa(k, j) = bigq_div(aa(k, j), pivot)
         end do
         do j = 1, nrhs
            x(k, j) = bigq_div(x(k, j), pivot)
         end do

         do i = 1, n
            if (i == k) cycle
            factor = aa(i, k)
            if (bigq_is_zero(factor)) cycle
            do j = k, n
               aa(i, j) = bigq_sub(aa(i, j), bigq_mul(factor, aa(k, j)))
            end do
            do j = 1, nrhs
               x(i, j) = bigq_sub(x(i, j), bigq_mul(factor, x(k, j)))
            end do
         end do
      end do
      if (present(info)) info = stat
   end subroutine bigq_solve

   pure subroutine bigq_inverse(a, inverse, info)
      type(bigq), intent(in) :: a(:, :) !! Square rational matrix to invert.
      type(bigq), allocatable, intent(out) :: inverse(:, :) !! Exact rational inverse matrix.
      integer, intent(out), optional :: info !! Zero on success; one for nonsquare input; two for singular input.
      type(bigq), allocatable :: identity(:, :)
      integer :: i, n, stat

      n = size(a, 1)
      if (size(a, 2) /= n) then
         allocate(inverse(0, 0))
         if (present(info)) info = 1
         return
      end if
      allocate(identity(n, n))
      identity = bigq_from_int64(0_int64)
      do i = 1, n
         identity(i, i) = bigq_from_int64(1_int64)
      end do
      call bigq_solve(a, identity, inverse, stat)
      if (present(info)) info = stat
   end subroutine bigq_inverse

   pure subroutine bigz_solve(a, b, x, info)
      type(bigz), intent(in) :: a(:, :) !! Square integer coefficient matrix with shape n by n.
      type(bigz), intent(in) :: b(:, :) !! Integer right-hand side matrix with shape n by nrhs.
      type(bigq), allocatable, intent(out) :: x(:, :) !! Exact rational solution, matching non-modular solve.bigz behavior.
      integer, intent(out), optional :: info !! Zero on success; one for dimension mismatch; two for singular input.
      type(bigq), allocatable :: aq(:, :), bq(:, :)
      integer :: i, j, stat

      allocate(aq(size(a, 1), size(a, 2)))
      allocate(bq(size(b, 1), size(b, 2)))
      do j = 1, size(a, 2)
         do i = 1, size(a, 1)
            aq(i, j) = bigq_make(a(i, j), bigz_from_int64(1_int64))
         end do
      end do
      do j = 1, size(b, 2)
         do i = 1, size(b, 1)
            bq(i, j) = bigq_make(b(i, j), bigz_from_int64(1_int64))
         end do
      end do
      call bigq_solve(aq, bq, x, stat)
      if (present(info)) info = stat
   end subroutine bigz_solve

   pure subroutine bigz_inverse(a, inverse, info)
      type(bigz), intent(in) :: a(:, :) !! Square integer matrix to invert over the rationals.
      type(bigq), allocatable, intent(out) :: inverse(:, :) !! Exact rational inverse for non-modular input.
      integer, intent(out), optional :: info !! Zero on success; one for nonsquare input; two for singular input.
      type(bigq), allocatable :: aq(:, :)
      integer :: i, j, stat

      allocate(aq(size(a, 1), size(a, 2)))
      do j = 1, size(a, 2)
         do i = 1, size(a, 1)
            aq(i, j) = bigq_make(a(i, j), bigz_from_int64(1_int64))
         end do
      end do
      call bigq_inverse(aq, inverse, stat)
      if (present(info)) info = stat
   end subroutine bigz_inverse

end module gmp_matrix
