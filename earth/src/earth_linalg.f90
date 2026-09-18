module earth_linalg
   use earth_kinds, only : dp
   implicit none
   private

   public :: least_squares
   public :: weighted_mean
   public :: weighted_centered_sos

contains

   subroutine least_squares(a, y, coefficients, fitted, residuals, rss, rank, leverages, ok)
      real(dp), intent(in) :: a(:, :) !! Design matrix with observations in rows and regression columns in columns.
      real(dp), intent(in) :: y(:, :) !! Response matrix with the same observation count as a and one or more responses.
      real(dp), allocatable, intent(out) :: coefficients(:, :) !! Least-squares coefficient matrix, one column per response.
      real(dp), allocatable, intent(out) :: fitted(:, :) !! Fitted response values at the rows of the supplied design matrix.
      real(dp), allocatable, intent(out) :: residuals(:, :) !! Response residuals y minus fitted values for every response.
      real(dp), intent(out) :: rss !! Residual sum of squares accumulated across all response columns.
      integer, intent(out) :: rank !! Numerical rank detected by modified Gram-Schmidt without column pivoting.
      real(dp), allocatable, intent(out) :: leverages(:) !! Diagonal of the projection matrix for a full-rank design.
      logical, intent(out) :: ok !! True when the design is full rank and the least-squares solution was formed.

      integer :: n
      integer :: p
      integer :: nr
      integer :: i
      integer :: j
      integer :: k
      integer :: ir
      real(dp) :: col_norm
      real(dp) :: v_norm
      real(dp) :: proj
      real(dp) :: tol
      real(dp), allocatable :: q(:, :)
      real(dp), allocatable :: r(:, :)
      real(dp), allocatable :: v(:)
      real(dp), allocatable :: qty(:, :)

      n = size(a, 1)
      p = size(a, 2)
      nr = size(y, 2)
      tol = 100.0_dp * epsilon(1.0_dp)
      ok = .false.
      rank = 0
      rss = huge(1.0_dp)

      if (n < 1 .or. p < 1 .or. nr < 1) return
      if (size(y, 1) /= n) return
      if (p > n) return

      allocate(q(n, p))
      allocate(r(p, p))
      allocate(v(n))
      q = 0.0_dp
      r = 0.0_dp

      do j = 1, p
         v = a(:, j)
         do k = 1, j - 1
            proj = dot_product(q(:, k), v)
            r(k, j) = r(k, j) + proj
            v = v - proj * q(:, k)
         end do
         do k = 1, j - 1
            proj = dot_product(q(:, k), v)
            r(k, j) = r(k, j) + proj
            v = v - proj * q(:, k)
         end do
         col_norm = sqrt(sum(a(:, j) * a(:, j)))
         v_norm = sqrt(sum(v * v))
         if (v_norm <= tol * max(1.0_dp, col_norm)) then
            rank = j - 1
            return
         end if
         r(j, j) = v_norm
         q(:, j) = v / v_norm
         rank = j
      end do

      allocate(qty(p, nr))
      qty = matmul(transpose(q), y)
      allocate(coefficients(p, nr))
      coefficients = qty

      do ir = 1, nr
         do i = p, 1, -1
            if (i < p) then
               coefficients(i, ir) = coefficients(i, ir) - &
                  dot_product(r(i, i + 1:p), coefficients(i + 1:p, ir))
            end if
            coefficients(i, ir) = coefficients(i, ir) / r(i, i)
         end do
      end do

      allocate(fitted(n, nr))
      allocate(residuals(n, nr))
      fitted = matmul(a, coefficients)
      residuals = y - fitted
      rss = sum(residuals * residuals)
      allocate(leverages(n))
      leverages = 0.0_dp
      do j = 1, p
         leverages = leverages + q(:, j) * q(:, j)
      end do
      ok = .true.
   end subroutine least_squares

   pure real(dp) function weighted_mean(x, w) result(mean_value)
      real(dp), intent(in) :: x(:) !! Numeric observations whose weighted arithmetic mean is required.
      real(dp), intent(in) :: w(:) !! Nonnegative case weights with the same length as x.

      real(dp) :: sw

      sw = sum(w)
      if (sw > 0.0_dp) then
         mean_value = dot_product(w, x) / sw
      else
         mean_value = 0.0_dp
      end if
   end function weighted_mean

   pure real(dp) function weighted_centered_sos(x, w) result(sos)
      real(dp), intent(in) :: x(:) !! Numeric observations whose weighted centered sum of squares is required.
      real(dp), intent(in) :: w(:) !! Nonnegative case weights with the same length as x.

      real(dp) :: center

      center = weighted_mean(x, w)
      sos = sum(w * (x - center) * (x - center))
   end function weighted_centered_sos

end module earth_linalg
