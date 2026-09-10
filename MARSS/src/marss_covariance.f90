! SPDX-License-Identifier: GPL-2.0-only
module marss_covariance
   use marss_kinds, only : dp
   implicit none
   private
   public :: psd_inverse
   public :: psd_inverse_logdet
   public :: psd_inverse_sqrt
   public :: psd_cholesky_standardize
   public :: psd_cholesky_inverse
   public :: is_positive_semidefinite
   public :: jacobi_symmetric_eigen

contains

   pure subroutine psd_inverse(a, inverse, rank, info)
      real(dp), intent(in) :: a(:, :) !! Symmetric positive-semidefinite matrix to pseudo-invert.
      real(dp), allocatable, intent(out) :: inverse(:, :) !! Moore-Penrose inverse formed from positive eigenvalues.
      integer, intent(out) :: rank !! Numerical rank under the package semidefinite tolerance.
      integer, intent(out) :: info !! Zero on success, nonzero for shape, convergence, or negative-eigenvalue failure.
      real(dp) :: logdet

      call psd_inverse_logdet(a, inverse, logdet, rank, info)
   end subroutine psd_inverse

   pure subroutine psd_inverse_logdet(a, inverse, logdet, rank, info)
      real(dp), intent(in) :: a(:, :) !! Symmetric positive-semidefinite matrix to pseudo-invert and factor spectrally.
      real(dp), allocatable, intent(out) :: inverse(:, :) !! Moore-Penrose inverse with zero weight on null directions.
      real(dp), intent(out) :: logdet !! Log pseudo-determinant, summing logs of positive eigenvalues only.
      integer, intent(out) :: rank !! Numerical rank under a scale-aware tolerance.
      integer, intent(out) :: info !! Zero on success, nonzero for nonsquare input, eigensolver failure, or indefiniteness.
      real(dp), allocatable :: values(:)
      real(dp), allocatable :: vectors(:, :)
      real(dp) :: scale
      real(dp) :: tol
      integer :: i
      integer :: n

      info = 0
      logdet = 0.0_dp
      rank = 0
      if (size(a, 1) /= size(a, 2)) then
         allocate(inverse(0, 0))
         info = 1
         return
      end if
      n = size(a, 1)
      allocate(inverse(n, n))
      inverse = 0.0_dp
      if (n == 0) return
      call jacobi_symmetric_eigen(a, values, vectors, info)
      if (info /= 0) return
      scale = max(1.0_dp, maxval(abs(values)))
      tol = 100.0_dp * epsilon(1.0_dp) * real(max(1, n), dp) * scale
      if (minval(values) < -tol) then
         info = 2
         return
      end if
      do i = 1, n
         if (values(i) > tol) then
            rank = rank + 1
            logdet = logdet + log(values(i))
            inverse = inverse + outer_product(vectors(:, i), vectors(:, i)) / values(i)
         end if
      end do
      inverse = 0.5_dp * (inverse + transpose(inverse))
   end subroutine psd_inverse_logdet


   pure subroutine psd_inverse_sqrt(a, inverse_sqrt, rank, info)
      real(dp), intent(in) :: a(:, :) !! Symmetric positive-semidefinite matrix whose spectral inverse square root is requested.
      real(dp), allocatable, intent(out) :: inverse_sqrt(:, :) !! Symmetric pseudoinverse square root with zero null-space weight.
      integer, intent(out) :: rank !! Numerical rank under the package scale-aware semidefinite tolerance.
      integer, intent(out) :: info !! Zero on success; nonzero for nonsquare input, eigensolver failure, or indefiniteness.
      real(dp), allocatable :: values(:)
      real(dp), allocatable :: vectors(:, :)
      real(dp) :: scale
      real(dp) :: tol
      integer :: i
      integer :: n

      info = 0
      rank = 0
      if (size(a, 1) /= size(a, 2)) then
         allocate(inverse_sqrt(0, 0))
         info = 1
         return
      end if
      n = size(a, 1)
      allocate(inverse_sqrt(n, n))
      inverse_sqrt = 0.0_dp
      if (n == 0) return
      call jacobi_symmetric_eigen(a, values, vectors, info)
      if (info /= 0) return
      scale = max(1.0_dp, maxval(abs(values)))
      tol = 100.0_dp * epsilon(1.0_dp) * real(max(1, n), dp) * scale
      if (minval(values) < -tol) then
         info = 2
         return
      end if
      do i = 1, n
         if (values(i) > tol) then
            rank = rank + 1
            inverse_sqrt = inverse_sqrt + outer_product(vectors(:, i), vectors(:, i)) / sqrt(values(i))
         end if
      end do
      inverse_sqrt = 0.5_dp * (inverse_sqrt + transpose(inverse_sqrt))
   end subroutine psd_inverse_sqrt


   pure subroutine psd_cholesky_standardize(a, x, y, info)
      real(dp), intent(in) :: a(:, :) !! Positive-semidefinite covariance matrix defining ordered Cholesky standardization.
      real(dp), intent(in) :: x(:) !! Residual vector ordered consistently with rows and columns of a.
      real(dp), allocatable, intent(out) :: y(:) !! Residual transformed by the inverse lower pseudo-Cholesky factor.
      integer, intent(out) :: info !! Zero on success; nonzero for dimensions or failure of the active principal factorization.
      real(dp), allocatable :: lower(:, :)
      real(dp), allocatable :: submatrix(:, :)
      real(dp), allocatable :: work(:, :)
      integer, allocatable :: active(:)
      real(dp) :: diagonal
      real(dp) :: scale
      real(dp) :: sum_value
      real(dp) :: tolerance
      integer :: i
      integer :: ia
      integer :: j
      integer :: ja
      integer :: k
      integer :: n
      integer :: nactive

      info = 0
      n = size(a, 1)
      if (size(a, 2) /= n .or. size(x) /= n) then
         allocate(y(0))
         info = 1
         return
      end if
      allocate(y(n))
      y = 0.0_dp
      if (n == 0) return
      allocate(work(n, n))
      work = 0.5_dp * (a + transpose(a))
      scale = max(1.0_dp, maxval(abs(work)))
      tolerance = sqrt(epsilon(1.0_dp)) * scale
      where (abs(work) < tolerance)
         work = 0.0_dp
      end where
      if (minval([(work(i, i), i=1,n)]) < -tolerance) then
         info = 2
         return
      end if
      nactive = count([(work(i, i) > tolerance, i=1,n)])
      allocate(active(nactive))
      k = 0
      do i = 1, n
         if (work(i, i) > tolerance) then
            k = k + 1
            active(k) = i
         end if
      end do
      if (nactive == 0) return
      allocate(submatrix(nactive, nactive), lower(nactive, nactive))
      do ja = 1, nactive
         j = active(ja)
         do ia = 1, nactive
            i = active(ia)
            submatrix(ia, ja) = work(i, j)
         end do
      end do
      lower = 0.0_dp
      do j = 1, nactive
         diagonal = submatrix(j, j)
         if (j > 1) diagonal = diagonal - dot_product(lower(j, 1:j - 1), lower(j, 1:j - 1))
         if (diagonal <= tolerance) then
            info = 3
            return
         end if
         lower(j, j) = sqrt(diagonal)
         do i = j + 1, nactive
            sum_value = submatrix(i, j)
            if (j > 1) sum_value = sum_value - dot_product(lower(i, 1:j - 1), lower(j, 1:j - 1))
            lower(i, j) = sum_value / lower(j, j)
         end do
      end do
      do ia = 1, nactive
         i = active(ia)
         sum_value = x(i)
         if (ia > 1) then
            do ja = 1, ia - 1
               sum_value = sum_value - lower(ia, ja) * y(active(ja))
            end do
         end if
         y(i) = sum_value / lower(ia, ia)
      end do
   end subroutine psd_cholesky_standardize


   pure subroutine psd_cholesky_inverse(a, inverse_factor, info)
      real(dp), intent(in) :: a(:, :) !! Positive-semidefinite covariance whose ordered lower pseudo-Cholesky inverse is requested.
      real(dp), allocatable, intent(out) :: inverse_factor(:, :) !! Ordered lower pseudo-Cholesky inverse transform.
      integer, intent(out) :: info !! Zero on success; nonzero when dimensions or active-factor construction fail.
      real(dp), allocatable :: basis(:)
      real(dp), allocatable :: column(:)
      integer :: j
      integer :: n

      info = 0
      n = size(a, 1)
      if (size(a, 2) /= n) then
         allocate(inverse_factor(0, 0))
         info = 1
         return
      end if
      allocate(inverse_factor(n, n), basis(n))
      inverse_factor = 0.0_dp
      basis = 0.0_dp
      do j = 1, n
         basis = 0.0_dp
         basis(j) = 1.0_dp
         call psd_cholesky_standardize(a, basis, column, info)
         if (info /= 0) return
         inverse_factor(:, j) = column
      end do
   end subroutine psd_cholesky_inverse

   pure logical function is_positive_semidefinite(a, tolerance) result(ok)
      real(dp), intent(in) :: a(:, :) !! Symmetric matrix tested for nonnegative eigenvalues within numerical tolerance.
      real(dp), intent(in), optional :: tolerance !! Absolute eigenvalue tolerance; default is scale-aware machine precision.
      real(dp), allocatable :: values(:)
      real(dp), allocatable :: vectors(:, :)
      real(dp) :: tol
      integer :: info

      ok = .false.
      if (size(a, 1) /= size(a, 2)) return
      if (size(a, 1) == 0) then
         ok = .true.
         return
      end if
      call jacobi_symmetric_eigen(a, values, vectors, info)
      if (info /= 0) return
      tol = 100.0_dp * epsilon(1.0_dp) * real(max(1, size(a, 1)), dp) * max(1.0_dp, maxval(abs(values)))
      if (present(tolerance)) tol = tolerance
      ok = minval(values) >= -tol
   end function is_positive_semidefinite

   pure subroutine jacobi_symmetric_eigen(a, values, vectors, info)
      real(dp), intent(in) :: a(:, :) !! Real symmetric matrix whose eigenpairs are computed by Jacobi rotations.
      real(dp), allocatable, intent(out) :: values(:) !! Eigenvalues in unsorted order.
      real(dp), allocatable, intent(out) :: vectors(:, :) !! Corresponding orthonormal eigenvectors stored by column.
      integer, intent(out) :: info !! Zero on convergence, nonzero for shape mismatch or iteration-limit failure.
      real(dp), allocatable :: work(:, :)
      real(dp) :: app
      real(dp) :: apq
      real(dp) :: aqq
      real(dp) :: c
      real(dp) :: maxoff
      real(dp) :: phi
      real(dp) :: s
      real(dp) :: scale
      real(dp) :: temp
      real(dp) :: tol
      integer :: i
      integer :: iter
      integer :: j
      integer :: max_iter
      integer :: n
      integer :: p
      integer :: q

      info = 0
      if (size(a, 1) /= size(a, 2)) then
         allocate(values(0), vectors(0, 0))
         info = 1
         return
      end if
      n = size(a, 1)
      allocate(values(n), vectors(n, n), work(n, n))
      if (n == 0) return
      work = 0.5_dp * (a + transpose(a))
      vectors = 0.0_dp
      do i = 1, n
         vectors(i, i) = 1.0_dp
      end do
      if (n == 1) then
         values(1) = work(1, 1)
         return
      end if
      scale = max(1.0_dp, maxval(abs(work)))
      tol = 50.0_dp * epsilon(1.0_dp) * scale
      max_iter = max(50, 100 * n * n)
      do iter = 1, max_iter
         maxoff = 0.0_dp
         p = 1
         q = 2
         do j = 2, n
            do i = 1, j - 1
               if (abs(work(i, j)) > maxoff) then
                  maxoff = abs(work(i, j))
                  p = i
                  q = j
               end if
            end do
         end do
         if (maxoff <= tol) exit
         app = work(p, p)
         aqq = work(q, q)
         apq = work(p, q)
         phi = 0.5_dp * atan2(2.0_dp * apq, aqq - app)
         c = cos(phi)
         s = sin(phi)
         do i = 1, n
            if (i /= p .and. i /= q) then
               temp = c * work(i, p) - s * work(i, q)
               work(i, q) = s * work(i, p) + c * work(i, q)
               work(q, i) = work(i, q)
               work(i, p) = temp
               work(p, i) = temp
            end if
         end do
         work(p, p) = c*c*app - 2.0_dp*s*c*apq + s*s*aqq
         work(q, q) = s*s*app + 2.0_dp*s*c*apq + c*c*aqq
         work(p, q) = 0.0_dp
         work(q, p) = 0.0_dp
         do i = 1, n
            temp = c * vectors(i, p) - s * vectors(i, q)
            vectors(i, q) = s * vectors(i, p) + c * vectors(i, q)
            vectors(i, p) = temp
         end do
      end do
      if (iter > max_iter) then
         info = 2
         return
      end if
      do i = 1, n
         values(i) = work(i, i)
      end do
   end subroutine jacobi_symmetric_eigen

   pure function outer_product(x, y) result(a)
      real(dp), intent(in) :: x(:) !! Left vector in the outer product.
      real(dp), intent(in) :: y(:) !! Right vector in the outer product.
      real(dp) :: a(size(x), size(y))

      a = spread(x, 2, size(y)) * spread(y, 1, size(x))
   end function outer_product

end module marss_covariance
