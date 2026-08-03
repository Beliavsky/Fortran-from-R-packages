! SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only
module sandwich_linalg
   use sandwich_kinds, only : dp
   use sandwich_status, only : SANDWICH_SUCCESS, SANDWICH_INVALID_ARGUMENT, &
      SANDWICH_DIMENSION_MISMATCH, SANDWICH_SINGULAR_MATRIX, &
      SANDWICH_NUMERICAL_FAILURE
   implicit none
   private

   public :: identity_matrix, inverse_matrix, solve_linear
   public :: ols_coefficients, covariance_matrix
   public :: symmetric_eigen_jacobi, symmetric_matrix_power, project_psd
   public :: trace_matrix, outer_product

   interface solve_linear
      module procedure solve_linear_vector
      module procedure solve_linear_matrix
   end interface solve_linear

contains

   function identity_matrix(n) result(a)
      integer, intent(in) :: n
      real(dp), allocatable :: a(:, :)
      integer :: i

      allocate(a(max(n, 0), max(n, 0)))
      a = 0.0_dp
      do i = 1, n
         a(i, i) = 1.0_dp
      end do
   end function identity_matrix

   real(dp) function trace_matrix(a) result(value)
      real(dp), intent(in) :: a(:, :)
      integer :: i

      value = 0.0_dp
      do i = 1, min(size(a, 1), size(a, 2))
         value = value + a(i, i)
      end do
   end function trace_matrix

   function outer_product(x, y) result(a)
      real(dp), intent(in) :: x(:), y(:)
      real(dp), allocatable :: a(:, :)
      integer :: i, j

      allocate(a(size(x), size(y)))
      do j = 1, size(y)
         do i = 1, size(x)
            a(i, j) = x(i) * y(j)
         end do
      end do
   end function outer_product

   subroutine solve_linear_vector(a, b, x, status)
      real(dp), intent(in) :: a(:, :), b(:)
      real(dp), allocatable, intent(out) :: x(:)
      integer, intent(out), optional :: status
      real(dp), allocatable :: bm(:, :), xm(:, :)
      integer :: info

      if (size(a, 1) /= size(a, 2) .or. size(b) /= size(a, 1)) then
         allocate(x(0))
         if (present(status)) status = SANDWICH_DIMENSION_MISMATCH
         return
      end if

      allocate(bm(size(b), 1))
      bm(:, 1) = b
      call solve_linear_matrix(a, bm, xm, info)
      if (info == SANDWICH_SUCCESS) then
         allocate(x(size(b)))
         x = xm(:, 1)
      else
         allocate(x(0))
      end if
      if (present(status)) status = info
   end subroutine solve_linear_vector

   subroutine solve_linear_matrix(a, b, x, status)
      real(dp), intent(in) :: a(:, :), b(:, :)
      real(dp), allocatable, intent(out) :: x(:, :)
      integer, intent(out), optional :: status
      real(dp), allocatable :: aa(:, :), bb(:, :), scale(:), tmp_row(:)
      real(dp) :: factor, pivot_abs, candidate, tiny_pivot, tmp
      integer :: n, nrhs, i, j, k, pivot, info

      n = size(a, 1)
      nrhs = size(b, 2)
      info = SANDWICH_SUCCESS

      if (n <= 0 .or. size(a, 2) /= n .or. size(b, 1) /= n) then
         allocate(x(0, 0))
         if (present(status)) status = SANDWICH_DIMENSION_MISMATCH
         return
      end if

      allocate(aa(n, n), bb(n, nrhs), scale(n), tmp_row(max(n, nrhs)))
      aa = a
      bb = b
      do i = 1, n
         scale(i) = maxval(abs(aa(i, :)))
      end do
      if (any(scale <= tiny(1.0_dp))) then
         allocate(x(0, 0))
         if (present(status)) status = SANDWICH_SINGULAR_MATRIX
         return
      end if

      tiny_pivot = 100.0_dp * epsilon(1.0_dp) * max(1.0_dp, maxval(abs(aa)))

      do k = 1, n - 1
         pivot = k
         pivot_abs = abs(aa(k, k)) / scale(k)
         do i = k + 1, n
            candidate = abs(aa(i, k)) / scale(i)
            if (candidate > pivot_abs) then
               pivot = i
               pivot_abs = candidate
            end if
         end do

         if (abs(aa(pivot, k)) <= tiny_pivot) then
            info = SANDWICH_SINGULAR_MATRIX
            exit
         end if

         if (pivot /= k) then
            tmp_row(1:n) = aa(k, :)
            aa(k, :) = aa(pivot, :)
            aa(pivot, :) = tmp_row(1:n)
            if (nrhs > 0) then
               tmp_row(1:nrhs) = bb(k, :)
               bb(k, :) = bb(pivot, :)
               bb(pivot, :) = tmp_row(1:nrhs)
            end if
            tmp = scale(k)
            scale(k) = scale(pivot)
            scale(pivot) = tmp
         end if

         do i = k + 1, n
            factor = aa(i, k) / aa(k, k)
            aa(i, k) = 0.0_dp
            aa(i, k + 1:n) = aa(i, k + 1:n) - factor * aa(k, k + 1:n)
            bb(i, :) = bb(i, :) - factor * bb(k, :)
         end do
      end do

      if (info == SANDWICH_SUCCESS) then
         if (abs(aa(n, n)) <= tiny_pivot) info = SANDWICH_SINGULAR_MATRIX
      end if

      if (info /= SANDWICH_SUCCESS) then
         allocate(x(0, 0))
         if (present(status)) status = info
         return
      end if

      allocate(x(n, nrhs))
      x = 0.0_dp
      do j = 1, nrhs
         do i = n, 1, -1
            x(i, j) = bb(i, j)
            if (i < n) x(i, j) = x(i, j) - dot_product(aa(i, i + 1:n), x(i + 1:n, j))
            x(i, j) = x(i, j) / aa(i, i)
         end do
      end do

      if (present(status)) status = info
   end subroutine solve_linear_matrix

   subroutine inverse_matrix(a, ainv, status)
      real(dp), intent(in) :: a(:, :)
      real(dp), allocatable, intent(out) :: ainv(:, :)
      integer, intent(out), optional :: status
      real(dp), allocatable :: eye(:, :)
      integer :: info

      if (size(a, 1) /= size(a, 2) .or. size(a, 1) <= 0) then
         allocate(ainv(0, 0))
         if (present(status)) status = SANDWICH_DIMENSION_MISMATCH
         return
      end if

      eye = identity_matrix(size(a, 1))
      call solve_linear_matrix(a, eye, ainv, info)
      if (present(status)) status = info
   end subroutine inverse_matrix

   subroutine ols_coefficients(x, y, beta, status, weights)
      real(dp), intent(in) :: x(:, :), y(:)
      real(dp), allocatable, intent(out) :: beta(:)
      integer, intent(out), optional :: status
      real(dp), intent(in), optional :: weights(:)
      real(dp), allocatable :: xtx(:, :), xty(:), xw(:, :), yw(:)
      real(dp) :: sw
      integer :: i, n, k, info

      n = size(x, 1)
      k = size(x, 2)
      if (n <= 0 .or. k <= 0 .or. size(y) /= n) then
         allocate(beta(0))
         if (present(status)) status = SANDWICH_DIMENSION_MISMATCH
         return
      end if
      if (present(weights)) then
         if (size(weights) /= n .or. any(weights < 0.0_dp)) then
            allocate(beta(0))
            if (present(status)) status = SANDWICH_INVALID_ARGUMENT
            return
         end if
         allocate(xw(n, k), yw(n))
         do i = 1, n
            sw = sqrt(weights(i))
            xw(i, :) = sw * x(i, :)
            yw(i) = sw * y(i)
         end do
         xtx = matmul(transpose(xw), xw)
         xty = matmul(transpose(xw), yw)
      else
         xtx = matmul(transpose(x), x)
         xty = matmul(transpose(x), y)
      end if

      call solve_linear_vector(xtx, xty, beta, info)
      if (present(status)) status = info
   end subroutine ols_coefficients

   subroutine covariance_matrix(samples, covariance, status, center)
      real(dp), intent(in) :: samples(:, :)
      real(dp), allocatable, intent(out) :: covariance(:, :)
      integer, intent(out), optional :: status
      real(dp), intent(in), optional :: center(:)
      real(dp), allocatable :: centered(:, :), means(:)
      integer :: n, k, i

      n = size(samples, 1)
      k = size(samples, 2)
      if (n < 2 .or. k <= 0) then
         allocate(covariance(0, 0))
         if (present(status)) status = SANDWICH_INVALID_ARGUMENT
         return
      end if

      allocate(centered(n, k), means(k))
      if (present(center)) then
         if (size(center) /= k) then
            allocate(covariance(0, 0))
            if (present(status)) status = SANDWICH_DIMENSION_MISMATCH
            return
         end if
         means = center
      else
         means = sum(samples, dim = 1) / real(n, dp)
      end if
      do i = 1, n
         centered(i, :) = samples(i, :) - means
      end do
      covariance = matmul(transpose(centered), centered) / real(n - 1, dp)
      covariance = 0.5_dp * (covariance + transpose(covariance))
      if (present(status)) status = SANDWICH_SUCCESS
   end subroutine covariance_matrix

   subroutine symmetric_eigen_jacobi(a, values, vectors, status, tolerance, max_sweeps)
      real(dp), intent(in) :: a(:, :)
      real(dp), allocatable, intent(out) :: values(:), vectors(:, :)
      integer, intent(out), optional :: status
      real(dp), intent(in), optional :: tolerance
      integer, intent(in), optional :: max_sweeps
      real(dp), allocatable :: work(:, :)
      real(dp) :: tol, max_off, app, aqq, apq, tau, t, c, s
      real(dp) :: aip, aiq, vip, viq
      integer :: n, sweep_limit, sweep, p, q, i, j, imax
      logical :: converged

      n = size(a, 1)
      if (n <= 0 .or. size(a, 2) /= n) then
         allocate(values(0), vectors(0, 0))
         if (present(status)) status = SANDWICH_DIMENSION_MISMATCH
         return
      end if

      tol = 100.0_dp * epsilon(1.0_dp)
      if (present(tolerance)) tol = max(tolerance, epsilon(1.0_dp))
      sweep_limit = max(50, 20 * n * n)
      if (present(max_sweeps)) sweep_limit = max_sweeps

      allocate(work(n, n))
      work = 0.5_dp * (a + transpose(a))
      vectors = identity_matrix(n)
      converged = .false.

      do sweep = 1, sweep_limit
         max_off = 0.0_dp
         p = 1
         q = 1
         do j = 2, n
            do i = 1, j - 1
               if (abs(work(i, j)) > max_off) then
                  max_off = abs(work(i, j))
                  p = i
                  q = j
               end if
            end do
         end do

         if (max_off <= tol * max(1.0_dp, maxval(abs(work)))) then
            converged = .true.
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
               aip = work(i, p)
               aiq = work(i, q)
               work(i, p) = c * aip - s * aiq
               work(p, i) = work(i, p)
               work(i, q) = s * aip + c * aiq
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

      allocate(values(n))
      do i = 1, n
         values(i) = work(i, i)
      end do

      do i = 1, n - 1
         imax = i
         do j = i + 1, n
            if (values(j) > values(imax)) imax = j
         end do
         if (imax /= i) then
            app = values(i)
            values(i) = values(imax)
            values(imax) = app
            work(:, 1) = vectors(:, i)
            vectors(:, i) = vectors(:, imax)
            vectors(:, imax) = work(:, 1)
         end if
      end do

      if (present(status)) then
         if (converged) then
            status = SANDWICH_SUCCESS
         else
            status = SANDWICH_NUMERICAL_FAILURE
         end if
      end if
   end subroutine symmetric_eigen_jacobi

   subroutine symmetric_matrix_power(a, power, result_matrix, status, tolerance)
      real(dp), intent(in) :: a(:, :), power
      real(dp), allocatable, intent(out) :: result_matrix(:, :)
      integer, intent(out), optional :: status
      real(dp), intent(in), optional :: tolerance
      real(dp), allocatable :: values(:), vectors(:, :), scaled(:, :)
      real(dp) :: tol
      integer :: i, info

      tol = epsilon(1.0_dp)**(1.0_dp / 1.3_dp)
      if (present(tolerance)) tol = tolerance
      call symmetric_eigen_jacobi(a, values, vectors, info, tol)
      if (info /= SANDWICH_SUCCESS) then
         allocate(result_matrix(0, 0))
         if (present(status)) status = info
         return
      end if

      if (power < 0.0_dp .and. any(values <= tol)) then
         allocate(result_matrix(0, 0))
         if (present(status)) status = SANDWICH_SINGULAR_MATRIX
         return
      end if

      allocate(scaled(size(vectors, 1), size(vectors, 2)))
      scaled = vectors
      do i = 1, size(values)
         if (values(i) < tol) values(i) = 0.0_dp
         scaled(:, i) = scaled(:, i) * values(i)**power
      end do
      result_matrix = matmul(scaled, transpose(vectors))
      result_matrix = 0.5_dp * (result_matrix + transpose(result_matrix))
      if (present(status)) status = SANDWICH_SUCCESS
   end subroutine symmetric_matrix_power

   subroutine project_psd(a, result_matrix, status, tolerance)
      real(dp), intent(in) :: a(:, :)
      real(dp), allocatable, intent(out) :: result_matrix(:, :)
      integer, intent(out), optional :: status
      real(dp), intent(in), optional :: tolerance
      real(dp), allocatable :: values(:), vectors(:, :), scaled(:, :)
      real(dp) :: tol
      integer :: i, info

      tol = 100.0_dp * epsilon(1.0_dp)
      if (present(tolerance)) tol = tolerance
      call symmetric_eigen_jacobi(a, values, vectors, info, tol)
      if (info /= SANDWICH_SUCCESS) then
         allocate(result_matrix(0, 0))
         if (present(status)) status = info
         return
      end if

      scaled = vectors
      do i = 1, size(values)
         values(i) = max(values(i), 0.0_dp)
         scaled(:, i) = scaled(:, i) * values(i)
      end do
      result_matrix = matmul(scaled, transpose(vectors))
      result_matrix = 0.5_dp * (result_matrix + transpose(result_matrix))
      if (present(status)) status = SANDWICH_SUCCESS
   end subroutine project_psd

end module sandwich_linalg
