! Copyright (C) 1998 Douglas M. Bates and William N. Venables.
! Modern Fortran translation, 2026.
! SPDX-License-Identifier: GPL-2.0-or-later
module splines_linalg
   use splines_kinds, only : dp
   implicit none
   private
   public :: solve_linear, nullspace_transform

contains

   subroutine solve_linear(a, b, x, status)
      real(dp), intent(in) :: a(:, :), b(:)
      real(dp), allocatable, intent(out) :: x(:)
      integer, intent(out), optional :: status
      real(dp), allocatable :: lu(:, :), rhs(:), row(:)
      real(dp) :: pivot, factor, scale
      integer :: n, i, k, p

      if (present(status)) status = 0
      n = size(a, 1)
      allocate(x(max(0, n)))
      if (n == 0 .or. size(a, 2) /= n .or. size(b) /= n) then
         if (present(status)) status = 1
         if (n > 0) x = 0.0_dp
         return
      end if

      allocate(lu(n, n), rhs(n), row(n))
      lu = a
      rhs = b
      scale = max(1.0_dp, maxval(abs(lu)))

      do k = 1, n - 1
         p = k - 1 + maxloc(abs(lu(k:n, k)), dim=1)
         pivot = abs(lu(p, k))
         if (pivot <= epsilon(1.0_dp) * scale * real(n, dp)) then
            if (present(status)) status = 2
            x = 0.0_dp
            return
         end if
         if (p /= k) then
            row = lu(k, :)
            lu(k, :) = lu(p, :)
            lu(p, :) = row
            pivot = rhs(k)
            rhs(k) = rhs(p)
            rhs(p) = pivot
         end if
         do i = k + 1, n
            factor = lu(i, k) / lu(k, k)
            lu(i, k) = 0.0_dp
            lu(i, k + 1:n) = lu(i, k + 1:n) - factor * lu(k, k + 1:n)
            rhs(i) = rhs(i) - factor * rhs(k)
         end do
      end do

      if (abs(lu(n, n)) <= epsilon(1.0_dp) * scale * real(n, dp)) then
         if (present(status)) status = 2
         x = 0.0_dp
         return
      end if

      do i = n, 1, -1
         x(i) = rhs(i)
         if (i < n) x(i) = x(i) - dot_product(lu(i, i + 1:n), x(i + 1:n))
         x(i) = x(i) / lu(i, i)
      end do
   end subroutine solve_linear

   subroutine nullspace_transform(constraints, basis, transformed, status)
      ! Return coordinates of each row of BASIS in an orthonormal basis for
      ! null(transpose(CONSTRAINTS)).  CONSTRAINTS has m rows and p columns.
      real(dp), intent(in) :: constraints(:, :), basis(:, :)
      real(dp), allocatable, intent(out) :: transformed(:, :)
      integer, intent(out), optional :: status
      real(dp), allocatable :: q(:, :), v(:), candidate(:)
      real(dp) :: normv, tol
      integer :: p, m, i, j, k, rank, col

      if (present(status)) status = 0
      m = size(constraints, 1)
      p = size(constraints, 2)
      if (size(basis, 2) /= p .or. m > p) then
         allocate(transformed(0, 0))
         if (present(status)) status = 1
         return
      end if

      allocate(q(p, p), v(p), candidate(p))
      q = 0.0_dp
      rank = 0
      tol = 100.0_dp * epsilon(1.0_dp) * max(1.0_dp, maxval(abs(constraints)))

      ! Modified Gram-Schmidt for the constraint directions.
      do j = 1, m
         v = constraints(j, :)
         do k = 1, rank
            v = v - dot_product(q(:, k), v) * q(:, k)
         end do
         normv = sqrt(dot_product(v, v))
         if (normv > tol) then
            rank = rank + 1
            q(:, rank) = v / normv
         end if
      end do

      ! Deterministically complete the orthonormal basis with coordinate axes.
      do i = 1, p
         candidate = 0.0_dp
         candidate(i) = 1.0_dp
         do k = 1, rank
            candidate = candidate - dot_product(q(:, k), candidate) * q(:, k)
         end do
         ! Reorthogonalize.
         do k = 1, rank
            candidate = candidate - dot_product(q(:, k), candidate) * q(:, k)
         end do
         normv = sqrt(dot_product(candidate, candidate))
         if (normv > tol) then
            rank = rank + 1
            q(:, rank) = candidate / normv
            if (rank == p) exit
         end if
      end do

      if (rank /= p) then
         allocate(transformed(0, 0))
         if (present(status)) status = 2
         return
      end if

      allocate(transformed(size(basis, 1), p - m))
      do col = 1, p - m
         transformed(:, col) = matmul(basis, q(:, m + col))
      end do
   end subroutine nullspace_transform

end module splines_linalg
