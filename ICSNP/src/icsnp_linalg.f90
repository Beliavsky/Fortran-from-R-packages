! SPDX-License-Identifier: GPL-2.0-or-later
module icsnp_linalg
   use icsnp_kinds, only : dp
   use icsnp_status, only : icsnp_ok, icsnp_invalid_input, icsnp_singular, &
      icsnp_iteration_limit
   implicit none
   private
   public :: solve_linear, invert_matrix, determinant, symmetric_eigen
   public :: matrix_sqrt, matrix_inv_sqrt, covariance_matrix, mahalanobis_squared
   public :: frobenius_norm, trace_matrix, identity_matrix, sample_mean

contains

   pure real(dp) function frobenius_norm(a) result(value)
      real(dp), intent(in) :: a(:,:)
      value = sqrt(sum(a * a))
   end function frobenius_norm

   pure real(dp) function trace_matrix(a) result(value)
      real(dp), intent(in) :: a(:,:)
      integer :: i, n
      n = min(size(a, 1), size(a, 2))
      value = 0.0_dp
      do i = 1, n
         value = value + a(i, i)
      end do
   end function trace_matrix

   pure function identity_matrix(n) result(a)
      integer, intent(in) :: n
      real(dp) :: a(n, n)
      integer :: i
      a = 0.0_dp
      do i = 1, n
         a(i, i) = 1.0_dp
      end do
   end function identity_matrix

   pure function sample_mean(x) result(mean_value)
      real(dp), intent(in) :: x(:,:)
      real(dp) :: mean_value(size(x, 2))
      if (size(x, 1) > 0) then
         mean_value = sum(x, dim=1) / real(size(x, 1), dp)
      else
         mean_value = 0.0_dp
      end if
   end function sample_mean

   subroutine solve_linear(a, b, x, status)
      real(dp), intent(in) :: a(:,:), b(:)
      real(dp), allocatable, intent(out) :: x(:)
      integer, intent(out) :: status
      real(dp), allocatable :: aug(:,:)
      real(dp) :: pivot, factor, scale
      integer :: n, i, j, k, pivot_row

      n = size(b)
      status = icsnp_invalid_input
      allocate(x(0))
      if (n < 1 .or. size(a, 1) /= n .or. size(a, 2) /= n) return

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
            status = icsnp_singular
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
      status = icsnp_ok
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
      integer :: n, j, local_status

      n = size(a, 1)
      allocate(ainv(0, 0))
      status = icsnp_invalid_input
      if (n < 1 .or. size(a, 2) /= n) return
      deallocate(ainv)
      allocate(ainv(n, n), rhs(n))
      do j = 1, n
         rhs = 0.0_dp
         rhs(j) = 1.0_dp
         call solve_linear(a, rhs, column, local_status)
         if (local_status /= icsnp_ok) then
            ainv = 0.0_dp
            status = local_status
            return
         end if
         ainv(:, j) = column
      end do
      status = icsnp_ok
   end subroutine invert_matrix

   subroutine determinant(a, det_value, status)
      real(dp), intent(in) :: a(:,:)
      real(dp), intent(out) :: det_value
      integer, intent(out) :: status
      real(dp), allocatable :: work(:,:)
      real(dp) :: factor, scale
      integer :: n, i, j, k, pivot_row, sign_value

      n = size(a, 1)
      det_value = 0.0_dp
      status = icsnp_invalid_input
      if (n < 1 .or. size(a, 2) /= n) return
      allocate(work(n, n))
      work = a
      scale = max(1.0_dp, maxval(abs(a)))
      sign_value = 1
      do k = 1, n - 1
         pivot_row = k - 1 + maxloc(abs(work(k:n, k)), dim=1)
         if (abs(work(pivot_row, k)) <= 100.0_dp * epsilon(1.0_dp) * scale) then
            status = icsnp_singular
            return
         end if
         if (pivot_row /= k) then
            call swap_rows(work, k, pivot_row)
            sign_value = -sign_value
         end if
         do i = k + 1, n
            factor = work(i, k) / work(k, k)
            do j = k + 1, n
               work(i, j) = work(i, j) - factor * work(k, j)
            end do
         end do
      end do
      if (abs(work(n, n)) <= 100.0_dp * epsilon(1.0_dp) * scale) then
         status = icsnp_singular
         return
      end if
      det_value = real(sign_value, dp)
      do i = 1, n
         det_value = det_value * work(i, i)
      end do
      status = icsnp_ok
   end subroutine determinant

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
      allocate(values(0), vectors(0, 0))
      status = icsnp_invalid_input
      if (n < 1 .or. size(a, 2) /= n) return
      tol = 100.0_dp * epsilon(1.0_dp)
      if (present(tolerance)) tol = max(tolerance, epsilon(1.0_dp))
      sweeps = max(50, 30 * n * n)
      if (present(max_sweeps)) sweeps = max_sweeps

      deallocate(values, vectors)
      allocate(values(n), vectors(n, n), work(n, n))
      work = 0.5_dp * (a + transpose(a))
      vectors = identity_matrix(n)
      status = icsnp_iteration_limit
      if (n == 1) then
         values(1) = work(1, 1)
         status = icsnp_ok
         return
      end if

      do sweep = 1, sweeps
         off_max = 0.0_dp
         p = 1
         q = 2
         do i = 1, n - 1
            if (maxval(abs(work(i, i + 1:n))) > off_max) then
               q = i + maxloc(abs(work(i, i + 1:n)), dim=1)
               p = i
               off_max = abs(work(p, q))
            end if
         end do
         if (off_max <= tol * max(1.0_dp, maxval(abs(work)))) then
            status = icsnp_ok
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
      real(dp) :: temp_value
      real(dp) :: temp_vector(size(vectors, 1))
      integer :: i, j, k
      do i = 1, size(values) - 1
         k = i
         do j = i + 1, size(values)
            if (values(j) > values(k)) k = j
         end do
         if (k /= i) then
            temp_value = values(i)
            values(i) = values(k)
            values(k) = temp_value
            temp_vector = vectors(:, i)
            vectors(:, i) = vectors(:, k)
            vectors(:, k) = temp_vector
         end if
      end do
   end subroutine sort_descending

   subroutine matrix_sqrt(a, root, status)
      real(dp), intent(in) :: a(:,:)
      real(dp), allocatable, intent(out) :: root(:,:)
      integer, intent(out) :: status
      real(dp), allocatable :: values(:), vectors(:,:), scaled(:,:)
      integer :: i, n
      call symmetric_eigen(a, values, vectors, status)
      if (status /= icsnp_ok) then
         allocate(root(0, 0))
         return
      end if
      n = size(values)
      if (minval(values) < -1000.0_dp * epsilon(1.0_dp) * max(1.0_dp, maxval(abs(values)))) then
         allocate(root(0, 0))
         status = icsnp_singular
         return
      end if
      allocate(scaled(n, n), root(n, n))
      scaled = vectors
      do i = 1, n
         scaled(:, i) = scaled(:, i) * sqrt(max(values(i), 0.0_dp))
      end do
      root = matmul(scaled, transpose(vectors))
      root = 0.5_dp * (root + transpose(root))
   end subroutine matrix_sqrt

   subroutine matrix_inv_sqrt(a, root, status)
      real(dp), intent(in) :: a(:,:)
      real(dp), allocatable, intent(out) :: root(:,:)
      integer, intent(out) :: status
      real(dp), allocatable :: values(:), vectors(:,:), scaled(:,:)
      real(dp) :: threshold
      integer :: i, n
      call symmetric_eigen(a, values, vectors, status)
      if (status /= icsnp_ok) then
         allocate(root(0, 0))
         return
      end if
      n = size(values)
      threshold = 1000.0_dp * epsilon(1.0_dp) * max(1.0_dp, maxval(abs(values)))
      if (minval(values) <= threshold) then
         allocate(root(0, 0))
         status = icsnp_singular
         return
      end if
      allocate(scaled(n, n), root(n, n))
      scaled = vectors
      do i = 1, n
         scaled(:, i) = scaled(:, i) / sqrt(values(i))
      end do
      root = matmul(scaled, transpose(vectors))
      root = 0.5_dp * (root + transpose(root))
   end subroutine matrix_inv_sqrt

   subroutine covariance_matrix(x, cov, status, center)
      real(dp), intent(in) :: x(:,:)
      real(dp), allocatable, intent(out) :: cov(:,:)
      integer, intent(out) :: status
      real(dp), intent(in), optional :: center(:)
      real(dp), allocatable :: centered(:,:)
      real(dp) :: mean_value(size(x, 2))
      integer :: n, p
      n = size(x, 1)
      p = size(x, 2)
      allocate(cov(0, 0))
      status = icsnp_invalid_input
      if (n < 2 .or. p < 1) return
      if (present(center)) then
         if (size(center) /= p) return
         mean_value = center
      else
         mean_value = sample_mean(x)
      end if
      allocate(centered(n, p))
      centered = x - spread(mean_value, dim=1, ncopies=n)
      deallocate(cov)
      allocate(cov(p, p))
      cov = matmul(transpose(centered), centered) / real(n - 1, dp)
      cov = 0.5_dp * (cov + transpose(cov))
      status = icsnp_ok
   end subroutine covariance_matrix

   subroutine mahalanobis_squared(x, center, covariance, distances, status)
      real(dp), intent(in) :: x(:,:), center(:), covariance(:,:)
      real(dp), allocatable, intent(out) :: distances(:)
      integer, intent(out) :: status
      real(dp), allocatable :: inv_cov(:,:), diff(:), transformed(:)
      integer :: i, n, p
      n = size(x, 1)
      p = size(x, 2)
      allocate(distances(0))
      status = icsnp_invalid_input
      if (n < 1 .or. p < 1 .or. size(center) /= p) return
      if (size(covariance, 1) /= p .or. size(covariance, 2) /= p) return
      call invert_matrix(covariance, inv_cov, status)
      if (status /= icsnp_ok) return
      deallocate(distances)
      allocate(distances(n), diff(p), transformed(p))
      do i = 1, n
         diff = x(i, :) - center
         transformed = matmul(inv_cov, diff)
         distances(i) = max(0.0_dp, dot_product(diff, transformed))
      end do
      status = icsnp_ok
   end subroutine mahalanobis_squared

end module icsnp_linalg
