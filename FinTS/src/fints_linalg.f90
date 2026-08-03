! SPDX-License-Identifier: GPL-2.0-or-later
module fints_linalg
   use fints_kinds, only : dp
   use fints_status, only : fints_ok, fints_invalid_input, fints_singular, &
      fints_iteration_limit
   implicit none
   private
   public :: solve_linear, invert_matrix, least_squares, symmetric_eigen

contains

   subroutine solve_linear(a, b, x, status)
      real(dp), intent(in) :: a(:,:), b(:)
      real(dp), allocatable, intent(out) :: x(:)
      integer, intent(out) :: status
      real(dp), allocatable :: aug(:,:)
      real(dp) :: pivot, factor, scale
      integer :: n, i, j, k, pivot_row

      n = size(b)
      status = fints_invalid_input
      allocate(x(0))
      if (size(a, 1) /= n .or. size(a, 2) /= n .or. n < 1) return

      deallocate(x)
      allocate(x(n), aug(n, n + 1))
      aug(:, 1:n) = a
      aug(:, n + 1) = b
      scale = max(1.0_dp, maxval(abs(a)))

      do k = 1, n
         pivot_row = k - 1 + maxloc(abs(aug(k:n, k)), dim=1)
         pivot = aug(pivot_row, k)
         if (abs(pivot) <= 100.0_dp * epsilon(1.0_dp) * scale) then
            x = 0.0_dp
            status = fints_singular
            return
         end if
         if (pivot_row /= k) call swap_rows(aug, k, pivot_row)

         do i = k + 1, n
            factor = aug(i, k) / aug(k, k)
            aug(i, k) = 0.0_dp
            do j = k + 1, n + 1
               aug(i, j) = aug(i, j) - factor * aug(k, j)
            end do
         end do
      end do

      do i = n, 1, -1
         x(i) = aug(i, n + 1)
         if (i < n) x(i) = x(i) - dot_product(aug(i, i + 1:n), x(i + 1:n))
         x(i) = x(i) / aug(i, i)
      end do
      status = fints_ok
   end subroutine solve_linear

   subroutine swap_rows(a, i, j)
      real(dp), intent(inout) :: a(:,:)
      integer, intent(in) :: i, j
      real(dp) :: temp(size(a, 2))

      temp = a(i, :)
      a(i, :) = a(j, :)
      a(j, :) = temp
   end subroutine swap_rows

   subroutine invert_matrix(a, ainv, status)
      real(dp), intent(in) :: a(:,:)
      real(dp), allocatable, intent(out) :: ainv(:,:)
      integer, intent(out) :: status
      real(dp), allocatable :: rhs(:), column(:)
      integer :: n, i, stat_local

      n = size(a, 1)
      status = fints_invalid_input
      allocate(ainv(0, 0))
      if (n < 1 .or. size(a, 2) /= n) return

      deallocate(ainv)
      allocate(ainv(n, n), rhs(n))
      ainv = 0.0_dp
      do i = 1, n
         rhs = 0.0_dp
         rhs(i) = 1.0_dp
         call solve_linear(a, rhs, column, stat_local)
         if (stat_local /= fints_ok) then
            status = stat_local
            return
         end if
         ainv(:, i) = column
      end do
      status = fints_ok
   end subroutine invert_matrix

   subroutine least_squares(design, y, beta, residuals, sse, status)
      real(dp), intent(in) :: design(:,:), y(:)
      real(dp), allocatable, intent(out) :: beta(:), residuals(:)
      real(dp), intent(out) :: sse
      integer, intent(out) :: status
      real(dp), allocatable :: xtx(:,:), xty(:)
      integer :: n, p

      n = size(y)
      p = size(design, 2)
      status = fints_invalid_input
      sse = huge(1.0_dp)
      allocate(beta(0), residuals(0))
      if (n < 1 .or. p < 1 .or. size(design, 1) /= n .or. n < p) return

      xtx = matmul(transpose(design), design)
      xty = matmul(transpose(design), y)
      call solve_linear(xtx, xty, beta, status)
      if (status /= fints_ok) return
      deallocate(residuals)
      allocate(residuals(n))
      residuals = y - matmul(design, beta)
      sse = dot_product(residuals, residuals)
   end subroutine least_squares

   subroutine symmetric_eigen(a, values, vectors, status, tolerance, max_sweeps)
      real(dp), intent(in) :: a(:,:)
      real(dp), allocatable, intent(out) :: values(:), vectors(:,:)
      integer, intent(out) :: status
      real(dp), intent(in), optional :: tolerance
      integer, intent(in), optional :: max_sweeps
      real(dp), allocatable :: work(:,:)
      real(dp) :: tol, app, aqq, apq, tau, t, c, s, off_max
      real(dp) :: wip, wiq, vip, viq
      integer :: n, p, q, i, sweep, sweeps

      n = size(a, 1)
      status = fints_invalid_input
      allocate(values(0), vectors(0, 0))
      if (n < 1 .or. size(a, 2) /= n) return

      tol = 100.0_dp * epsilon(1.0_dp)
      if (present(tolerance)) tol = max(tolerance, epsilon(1.0_dp))
      sweeps = max(50, 20 * n * n)
      if (present(max_sweeps)) sweeps = max_sweeps

      deallocate(values, vectors)
      allocate(values(n), vectors(n, n), work(n, n))
      work = 0.5_dp * (a + transpose(a))
      vectors = 0.0_dp
      do i = 1, n
         vectors(i, i) = 1.0_dp
      end do

      status = fints_iteration_limit
      do sweep = 1, sweeps
         off_max = 0.0_dp
         p = 1
         q = min(2, n)
         do i = 1, n - 1
            if (maxval(abs(work(i, i + 1:n))) > off_max) then
               q = i + maxloc(abs(work(i, i + 1:n)), dim=1)
               p = i
               off_max = abs(work(p, q))
            end if
         end do
         if (off_max <= tol * max(1.0_dp, maxval(abs(work)))) then
            status = fints_ok
            exit
         end if

         app = work(p, p)
         aqq = work(q, q)
         apq = work(p, q)
         tau = (aqq - app) / (2.0_dp * apq)
         if (tau >= 0.0_dp) then
            t = 1.0_dp / (tau + sqrt(1.0_dp + tau * tau))
         else
            t = -1.0_dp / (-tau + sqrt(1.0_dp + tau * tau))
         end if
         c = 1.0_dp / sqrt(1.0_dp + t * t)
         s = t * c

         do i = 1, n
            if (i /= p .and. i /= q) then
               wip = work(i, p)
               wiq = work(i, q)
               work(i, p) = c * wip - s * wiq
               work(p, i) = work(i, p)
               work(i, q) = s * wip + c * wiq
               work(q, i) = work(i, q)
            end if
         end do
         work(p, p) = c * c * app - 2.0_dp * s * c * apq + s * s * aqq
         work(q, q) = s * s * app + 2.0_dp * s * c * apq + c * c * aqq
         work(p, q) = 0.0_dp
         work(q, p) = 0.0_dp

         do i = 1, n
            vip = vectors(i, p)
            viq = vectors(i, q)
            vectors(i, p) = c * vip - s * viq
            vectors(i, q) = s * vip + c * viq
         end do
      end do

      do i = 1, n
         values(i) = work(i, i)
      end do
      call sort_descending(values, vectors)
   end subroutine symmetric_eigen

   subroutine sort_descending(values, vectors)
      real(dp), intent(inout) :: values(:), vectors(:,:)
      real(dp) :: value_temp
      real(dp) :: vector_temp(size(vectors, 1))
      integer :: i, j, k, n

      n = size(values)
      do i = 1, n - 1
         k = i
         do j = i + 1, n
            if (values(j) > values(k)) k = j
         end do
         if (k /= i) then
            value_temp = values(i)
            values(i) = values(k)
            values(k) = value_temp
            vector_temp = vectors(:, i)
            vectors(:, i) = vectors(:, k)
            vectors(:, k) = vector_temp
         end if
      end do
   end subroutine sort_descending

end module fints_linalg
