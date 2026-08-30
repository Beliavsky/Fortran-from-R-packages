module learnbayes_linalg
   use learnbayes_kinds, only: dp
   implicit none
   private

   public :: cholesky_lower
   public :: determinant_spd
   public :: inverse_matrix
   public :: least_squares
   public :: solve_linear
   public :: quadratic_form

contains

   subroutine cholesky_lower(a, l, info)
      real(dp), intent(in) :: a(:, :) !! Symmetric positive-definite matrix to factor.
      real(dp), intent(out) :: l(:, :) !! Lower-triangular Cholesky factor satisfying A = L*transpose(L).
      integer, intent(out) :: info !! Zero on success; otherwise the first non-positive pivot index.
      integer :: i
      integer :: j
      integer :: k
      real(dp) :: s
      integer :: n

      n = size(a, 1)
      l = 0.0_dp
      info = 0
      if (size(a, 2) /= n .or. size(l, 1) /= n .or. size(l, 2) /= n) then
         info = -1
         return
      end if
      do i = 1, n
         do j = 1, i
            s = a(i, j)
            do k = 1, j - 1
               s = s - l(i, k)*l(j, k)
            end do
            if (i == j) then
               if (s <= 0.0_dp) then
                  info = i
                  return
               end if
               l(i, j) = sqrt(s)
            else
               l(i, j) = s/l(j, j)
            end if
         end do
      end do
   end subroutine cholesky_lower

   subroutine solve_linear(a, b, x, info)
      real(dp), intent(in) :: a(:, :) !! Square coefficient matrix.
      real(dp), intent(in) :: b(:) !! Right-hand-side vector with length equal to the matrix order.
      real(dp), intent(out) :: x(:) !! Solution vector to A*x = b.
      integer, intent(out) :: info !! Zero on success; nonzero if dimensions fail or a pivot is singular.
      real(dp), allocatable :: aa(:, :)
      real(dp), allocatable :: bb(:)
      real(dp) :: factor
      real(dp) :: pivot
      real(dp) :: tmp
      integer :: i
      integer :: j
      integer :: k
      integer :: imax
      integer :: n

      n = size(a, 1)
      info = 0
      if (size(a, 2) /= n .or. size(b) /= n .or. size(x) /= n) then
         info = -1
         return
      end if
      allocate(aa(n, n), bb(n))
      aa = a
      bb = b
      do k = 1, n - 1
         imax = k
         pivot = abs(aa(k, k))
         do i = k + 1, n
            if (abs(aa(i, k)) > pivot) then
               pivot = abs(aa(i, k))
               imax = i
            end if
         end do
         if (pivot <= tiny(1.0_dp)) then
            info = k
            return
         end if
         if (imax /= k) then
            do j = k, n
               tmp = aa(k, j)
               aa(k, j) = aa(imax, j)
               aa(imax, j) = tmp
            end do
            tmp = bb(k)
            bb(k) = bb(imax)
            bb(imax) = tmp
         end if
         do i = k + 1, n
            factor = aa(i, k)/aa(k, k)
            aa(i, k) = 0.0_dp
            do j = k + 1, n
               aa(i, j) = aa(i, j) - factor*aa(k, j)
            end do
            bb(i) = bb(i) - factor*bb(k)
         end do
      end do
      if (abs(aa(n, n)) <= tiny(1.0_dp)) then
         info = n
         return
      end if
      x(n) = bb(n)/aa(n, n)
      do i = n - 1, 1, -1
         x(i) = (bb(i) - dot_product(aa(i, i + 1:n), x(i + 1:n)))/aa(i, i)
      end do
   end subroutine solve_linear

   subroutine inverse_matrix(a, ainv, info)
      real(dp), intent(in) :: a(:, :) !! Square matrix to invert.
      real(dp), intent(out) :: ainv(:, :) !! Matrix inverse on successful return.
      integer, intent(out) :: info !! Zero on success; nonzero if dimensions fail or the matrix is singular.
      real(dp), allocatable :: e(:)
      real(dp), allocatable :: col(:)
      integer :: j
      integer :: n

      n = size(a, 1)
      info = 0
      if (size(a, 2) /= n .or. size(ainv, 1) /= n .or. size(ainv, 2) /= n) then
         info = -1
         return
      end if
      allocate(e(n), col(n))
      do j = 1, n
         e = 0.0_dp
         e(j) = 1.0_dp
         call solve_linear(a, e, col, info)
         if (info /= 0) return
         ainv(:, j) = col
      end do
   end subroutine inverse_matrix

   function determinant_spd(a, info) result(value)
      real(dp), intent(in) :: a(:, :) !! Symmetric positive-definite matrix whose determinant is requested.
      integer, intent(out) :: info !! Zero on success; nonzero when Cholesky factorization fails.
      real(dp) :: value
      real(dp), allocatable :: l(:, :)
      integer :: i
      integer :: n

      n = size(a, 1)
      allocate(l(n, n))
      call cholesky_lower(a, l, info)
      if (info /= 0) then
         value = 0.0_dp
         return
      end if
      value = 1.0_dp
      do i = 1, n
         value = value*l(i, i)*l(i, i)
      end do
   end function determinant_spd

   function quadratic_form(x, a) result(value)
      real(dp), intent(in) :: x(:) !! Vector used on both sides of the quadratic form.
      real(dp), intent(in) :: a(:, :) !! Square matrix with order equal to size(x).
      real(dp) :: value

      value = dot_product(x, matmul(a, x))
   end function quadratic_form

   subroutine least_squares(x, y, beta, residuals, sigma2, xtx_inv, rank, info)
      real(dp), intent(in) :: x(:, :) !! Design matrix with observations in rows and predictors in columns.
      real(dp), intent(in) :: y(:) !! Response vector with one value for each design-matrix row.
      real(dp), intent(out) :: beta(:) !! Ordinary least-squares coefficient estimate.
      real(dp), intent(out) :: residuals(:) !! Residual vector y - X*beta.
      real(dp), intent(out) :: sigma2 !! Residual mean square using n-p degrees of freedom.
      real(dp), intent(out) :: xtx_inv(:, :) !! Inverse of transpose(X)*X for full-rank designs.
      integer, intent(out) :: rank !! Effective rank; this implementation returns p for full-rank designs or zero on failure.
      integer, intent(out) :: info !! Zero on success; nonzero for dimension or singularity failures.
      real(dp), allocatable :: xtx(:, :)
      real(dp), allocatable :: xty(:)
      integer :: n
      integer :: p

      n = size(x, 1)
      p = size(x, 2)
      info = 0
      rank = 0
      if (size(y) /= n .or. size(beta) /= p .or. size(residuals) /= n) then
         info = -1
         return
      end if
      if (size(xtx_inv, 1) /= p .or. size(xtx_inv, 2) /= p) then
         info = -2
         return
      end if
      allocate(xtx(p, p), xty(p))
      xtx = matmul(transpose(x), x)
      xty = matmul(transpose(x), y)
      call solve_linear(xtx, xty, beta, info)
      if (info /= 0) return
      call inverse_matrix(xtx, xtx_inv, info)
      if (info /= 0) return
      residuals = y - matmul(x, beta)
      if (n > p) then
         sigma2 = dot_product(residuals, residuals)/real(n - p, dp)
      else
         sigma2 = 0.0_dp
      end if
      rank = p
   end subroutine least_squares

end module learnbayes_linalg
