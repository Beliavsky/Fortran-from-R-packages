module matrixcalc_core
   use r_kinds, only : dp
   implicit none
   private

   public :: direct_prod, direct_sum, hadamard_prod, frobenius_prod
   public :: set_submatrix, shift_down, shift_up, shift_left, shift_right
   public :: vec, vech, lower_triangle, upper_triangle, matrix_trace
   public :: entrywise_norm, frobenius_norm, hilbert_schmidt_norm
   public :: maximum_norm, one_norm, inf_norm

   interface direct_prod
      module procedure direct_prod_mm
      module procedure direct_prod_vm
      module procedure direct_prod_mv
      module procedure direct_prod_vv
   end interface direct_prod

   interface direct_sum
      module procedure direct_sum_mm
      module procedure direct_sum_vm
      module procedure direct_sum_mv
      module procedure direct_sum_vv
   end interface direct_sum

   interface hadamard_prod
      module procedure hadamard_prod_mm
      module procedure hadamard_prod_vv
   end interface hadamard_prod

   interface frobenius_prod
      module procedure frobenius_prod_mm
      module procedure frobenius_prod_vv
   end interface frobenius_prod

   interface entrywise_norm
      module procedure entrywise_norm_matrix
      module procedure entrywise_norm_vector
   end interface entrywise_norm

   interface frobenius_norm
      module procedure frobenius_norm_matrix
      module procedure frobenius_norm_vector
   end interface frobenius_norm

   interface hilbert_schmidt_norm
      module procedure hilbert_schmidt_norm_matrix
      module procedure hilbert_schmidt_norm_vector
   end interface hilbert_schmidt_norm

   interface maximum_norm
      module procedure maximum_norm_matrix
      module procedure maximum_norm_vector
   end interface maximum_norm

   interface one_norm
      module procedure one_norm_matrix
      module procedure one_norm_vector
   end interface one_norm

   interface inf_norm
      module procedure inf_norm_matrix
      module procedure inf_norm_vector
   end interface inf_norm

contains

   pure function direct_prod_mm(x, y) result(z)
      real(dp), intent(in) :: x(:, :) !! Left matrix in the Kronecker product.
      real(dp), intent(in) :: y(:, :) !! Right matrix in the Kronecker product.
      real(dp), allocatable :: z(:, :)
      integer :: i
      integer :: j
      integer :: r1
      integer :: r2
      integer :: c1
      integer :: c2

      r1 = size(x, 1)
      c1 = size(x, 2)
      r2 = size(y, 1)
      c2 = size(y, 2)
      allocate(z(r1 * r2, c1 * c2))
      do j = 1, c1
         do i = 1, r1
            z((i - 1) * r2 + 1:i * r2, (j - 1) * c2 + 1:j * c2) = x(i, j) * y
         end do
      end do
   end function direct_prod_mm

   pure function direct_prod_vm(x, y) result(z)
      real(dp), intent(in) :: x(:) !! Left vector, interpreted as a column matrix as in R matrix(x).
      real(dp), intent(in) :: y(:, :) !! Right matrix in the Kronecker product.
      real(dp), allocatable :: z(:, :)
      real(dp), allocatable :: xm(:, :)

      allocate(xm(size(x), 1))
      xm(:, 1) = x
      z = direct_prod_mm(xm, y)
   end function direct_prod_vm

   pure function direct_prod_mv(x, y) result(z)
      real(dp), intent(in) :: x(:, :) !! Left matrix in the Kronecker product.
      real(dp), intent(in) :: y(:) !! Right vector, interpreted as a column matrix as in R matrix(y).
      real(dp), allocatable :: z(:, :)
      real(dp), allocatable :: ym(:, :)

      allocate(ym(size(y), 1))
      ym(:, 1) = y
      z = direct_prod_mm(x, ym)
   end function direct_prod_mv

   pure function direct_prod_vv(x, y) result(z)
      real(dp), intent(in) :: x(:) !! Left vector, interpreted as a column matrix.
      real(dp), intent(in) :: y(:) !! Right vector, interpreted as a column matrix.
      real(dp), allocatable :: z(:, :)
      real(dp), allocatable :: xm(:, :)
      real(dp), allocatable :: ym(:, :)

      allocate(xm(size(x), 1), ym(size(y), 1))
      xm(:, 1) = x
      ym(:, 1) = y
      z = direct_prod_mm(xm, ym)
   end function direct_prod_vv

   pure function direct_sum_mm(x, y) result(z)
      real(dp), intent(in) :: x(:, :) !! Upper-left matrix block.
      real(dp), intent(in) :: y(:, :) !! Lower-right matrix block.
      real(dp), allocatable :: z(:, :)
      integer :: rx
      integer :: cx
      integer :: ry
      integer :: cy

      rx = size(x, 1)
      cx = size(x, 2)
      ry = size(y, 1)
      cy = size(y, 2)
      allocate(z(rx + ry, cx + cy))
      z = 0.0_dp
      z(1:rx, 1:cx) = x
      z(rx + 1:rx + ry, cx + 1:cx + cy) = y
   end function direct_sum_mm

   pure function direct_sum_vm(x, y) result(z)
      real(dp), intent(in) :: x(:) !! Upper-left vector, interpreted as a column matrix.
      real(dp), intent(in) :: y(:, :) !! Lower-right matrix block.
      real(dp), allocatable :: z(:, :)
      real(dp), allocatable :: xm(:, :)

      allocate(xm(size(x), 1))
      xm(:, 1) = x
      z = direct_sum_mm(xm, y)
   end function direct_sum_vm

   pure function direct_sum_mv(x, y) result(z)
      real(dp), intent(in) :: x(:, :) !! Upper-left matrix block.
      real(dp), intent(in) :: y(:) !! Lower-right vector, interpreted as a column matrix.
      real(dp), allocatable :: z(:, :)
      real(dp), allocatable :: ym(:, :)

      allocate(ym(size(y), 1))
      ym(:, 1) = y
      z = direct_sum_mm(x, ym)
   end function direct_sum_mv

   pure function direct_sum_vv(x, y) result(z)
      real(dp), intent(in) :: x(:) !! Upper-left vector, interpreted as a column matrix.
      real(dp), intent(in) :: y(:) !! Lower-right vector, interpreted as a column matrix.
      real(dp), allocatable :: z(:, :)
      real(dp), allocatable :: xm(:, :)
      real(dp), allocatable :: ym(:, :)

      allocate(xm(size(x), 1), ym(size(y), 1))
      xm(:, 1) = x
      ym(:, 1) = y
      z = direct_sum_mm(xm, ym)
   end function direct_sum_vv

   pure function hadamard_prod_mm(x, y) result(z)
      real(dp), intent(in) :: x(:, :) !! First matrix; shape must match y.
      real(dp), intent(in) :: y(:, :) !! Second matrix; shape must match x.
      real(dp), allocatable :: z(:, :)

      allocate(z(size(x, 1), size(x, 2)))
      if (any(shape(x) /= shape(y))) then
         z = 0.0_dp
         return
      end if
      z = x * y
   end function hadamard_prod_mm

   pure function hadamard_prod_vv(x, y) result(z)
      real(dp), intent(in) :: x(:) !! First vector.
      real(dp), intent(in) :: y(:) !! Second vector; length must match x.
      real(dp), allocatable :: z(:, :)

      allocate(z(size(x), 1))
      if (size(x) /= size(y)) then
         z = 0.0_dp
         return
      end if
      z(:, 1) = x * y
   end function hadamard_prod_vv

   pure real(dp) function frobenius_prod_mm(x, y) result(value)
      real(dp), intent(in) :: x(:, :) !! First matrix; shape must match y.
      real(dp), intent(in) :: y(:, :) !! Second matrix; shape must match x.

      if (any(shape(x) /= shape(y))) then
         value = 0.0_dp
      else
         value = sum(x * y)
      end if
   end function frobenius_prod_mm

   pure real(dp) function frobenius_prod_vv(x, y) result(value)
      real(dp), intent(in) :: x(:) !! First vector.
      real(dp), intent(in) :: y(:) !! Second vector; length must match x.

      if (size(x) /= size(y)) then
         value = 0.0_dp
      else
         value = sum(x * y)
      end if
   end function frobenius_prod_vv

   pure function set_submatrix(x, y, row, col) result(z)
      real(dp), intent(in) :: x(:, :) !! Matrix to copy and modify.
      real(dp), intent(in) :: y(:, :) !! Submatrix inserted into x.
      integer, intent(in) :: row !! One-based first destination row.
      integer, intent(in) :: col !! One-based first destination column.
      real(dp), allocatable :: z(:, :)
      integer :: row2
      integer :: col2

      allocate(z(size(x, 1), size(x, 2)))
      z = x
      row2 = row + size(y, 1) - 1
      col2 = col + size(y, 2) - 1
      if (row < 1 .or. col < 1 .or. row2 > size(x, 1) .or. col2 > size(x, 2)) return
      z(row:row2, col:col2) = y
   end function set_submatrix

   pure function shift_down(a, rows, fill) result(b)
      real(dp), intent(in) :: a(:, :) !! Matrix to shift downward.
      integer, intent(in), optional :: rows !! Number of rows shifted; default is one.
      real(dp), intent(in), optional :: fill !! Fill value for new rows; default is zero.
      real(dp), allocatable :: b(:, :)
      integer :: nshift
      real(dp) :: fill_value
      integer :: nr

      nshift = 1
      if (present(rows)) nshift = max(0, rows)
      fill_value = 0.0_dp
      if (present(fill)) fill_value = fill
      nr = size(a, 1)
      allocate(b(nr, size(a, 2)))
      b = fill_value
      if (nshift < nr) b(nshift + 1:nr, :) = a(1:nr - nshift, :)
   end function shift_down

   pure function shift_up(a, rows, fill) result(b)
      real(dp), intent(in) :: a(:, :) !! Matrix to shift upward.
      integer, intent(in), optional :: rows !! Number of rows shifted; default is one.
      real(dp), intent(in), optional :: fill !! Fill value for new rows; default is zero.
      real(dp), allocatable :: b(:, :)
      integer :: nshift
      real(dp) :: fill_value
      integer :: nr

      nshift = 1
      if (present(rows)) nshift = max(0, rows)
      fill_value = 0.0_dp
      if (present(fill)) fill_value = fill
      nr = size(a, 1)
      allocate(b(nr, size(a, 2)))
      b = fill_value
      if (nshift < nr) b(1:nr - nshift, :) = a(nshift + 1:nr, :)
   end function shift_up

   pure function shift_left(a, cols, fill) result(b)
      real(dp), intent(in) :: a(:, :) !! Matrix to shift left.
      integer, intent(in), optional :: cols !! Number of columns shifted; default is one.
      real(dp), intent(in), optional :: fill !! Fill value for new columns; default is zero.
      real(dp), allocatable :: b(:, :)
      integer :: nshift
      real(dp) :: fill_value
      integer :: nc

      nshift = 1
      if (present(cols)) nshift = max(0, cols)
      fill_value = 0.0_dp
      if (present(fill)) fill_value = fill
      nc = size(a, 2)
      allocate(b(size(a, 1), nc))
      b = fill_value
      if (nshift < nc) b(:, 1:nc - nshift) = a(:, nshift + 1:nc)
   end function shift_left

   pure function shift_right(a, cols, fill) result(b)
      real(dp), intent(in) :: a(:, :) !! Matrix to shift right.
      integer, intent(in), optional :: cols !! Number of columns shifted; default is one.
      real(dp), intent(in), optional :: fill !! Fill value for new columns; default is zero.
      real(dp), allocatable :: b(:, :)
      integer :: nshift
      real(dp) :: fill_value
      integer :: nc

      nshift = 1
      if (present(cols)) nshift = max(0, cols)
      fill_value = 0.0_dp
      if (present(fill)) fill_value = fill
      nc = size(a, 2)
      allocate(b(size(a, 1), nc))
      b = fill_value
      if (nshift < nc) b(:, nshift + 1:nc) = a(:, 1:nc - nshift)
   end function shift_right

   pure function vec(x) result(v)
      real(dp), intent(in) :: x(:, :) !! Matrix whose columns are stacked in column-major order.
      real(dp), allocatable :: v(:)

      allocate(v(size(x)))
      v = reshape(x, [size(x)])
   end function vec

   pure function vech(x) result(v)
      real(dp), intent(in) :: x(:, :) !! Square matrix whose lower triangle is stacked by columns.
      real(dp), allocatable :: v(:)
      integer :: n
      integer :: i
      integer :: j
      integer :: k

      n = min(size(x, 1), size(x, 2))
      allocate(v(n * (n + 1) / 2))
      k = 0
      do j = 1, n
         do i = j, n
            k = k + 1
            v(k) = x(i, j)
         end do
      end do
   end function vech

   pure function lower_triangle(x) result(y)
      real(dp), intent(in) :: x(:, :) !! Square matrix whose upper triangle is zeroed.
      real(dp), allocatable :: y(:, :)
      integer :: i
      integer :: j

      allocate(y(size(x, 1), size(x, 2)))
      y = x
      do j = 1, size(y, 2)
         do i = 1, min(j - 1, size(y, 1))
            y(i, j) = 0.0_dp
         end do
      end do
   end function lower_triangle

   pure function upper_triangle(x) result(y)
      real(dp), intent(in) :: x(:, :) !! Square matrix whose lower triangle is zeroed.
      real(dp), allocatable :: y(:, :)
      integer :: i
      integer :: j

      allocate(y(size(x, 1), size(x, 2)))
      y = x
      do j = 1, size(y, 2)
         do i = j + 1, size(y, 1)
            y(i, j) = 0.0_dp
         end do
      end do
   end function upper_triangle

   pure real(dp) function matrix_trace(x) result(value)
      real(dp), intent(in) :: x(:, :) !! Square matrix whose diagonal is summed.
      integer :: i

      value = 0.0_dp
      do i = 1, min(size(x, 1), size(x, 2))
         value = value + x(i, i)
      end do
   end function matrix_trace

   pure real(dp) function entrywise_norm_matrix(x, p) result(value)
      real(dp), intent(in) :: x(:, :) !! Matrix whose entrywise p-norm is computed.
      real(dp), intent(in) :: p !! Nonzero exponent; positive infinity requests the maximum norm.

      if (.not. (p < huge(p))) then
         value = maxval(abs(x))
      else
         value = sum(abs(x)**p)**(1.0_dp / p)
      end if
   end function entrywise_norm_matrix

   pure real(dp) function entrywise_norm_vector(x, p) result(value)
      real(dp), intent(in) :: x(:) !! Vector whose entrywise p-norm is computed.
      real(dp), intent(in) :: p !! Nonzero exponent; positive infinity requests the maximum norm.

      if (.not. (p < huge(p))) then
         value = maxval(abs(x))
      else
         value = sum(abs(x)**p)**(1.0_dp / p)
      end if
   end function entrywise_norm_vector

   pure real(dp) function frobenius_norm_matrix(x) result(value)
      real(dp), intent(in) :: x(:, :) !! Matrix whose Frobenius norm is computed.

      value = sqrt(sum(x * x))
   end function frobenius_norm_matrix

   pure real(dp) function frobenius_norm_vector(x) result(value)
      real(dp), intent(in) :: x(:) !! Vector whose Euclidean/Frobenius norm is computed.

      value = sqrt(sum(x * x))
   end function frobenius_norm_vector

   pure real(dp) function hilbert_schmidt_norm_matrix(x) result(value)
      real(dp), intent(in) :: x(:, :) !! Matrix whose Hilbert-Schmidt norm is computed.

      value = sqrt(sum(x * x))
   end function hilbert_schmidt_norm_matrix

   pure real(dp) function hilbert_schmidt_norm_vector(x) result(value)
      real(dp), intent(in) :: x(:) !! Vector whose Hilbert-Schmidt norm is computed.

      value = sqrt(sum(x * x))
   end function hilbert_schmidt_norm_vector

   pure real(dp) function maximum_norm_matrix(x) result(value)
      real(dp), intent(in) :: x(:, :) !! Matrix whose largest absolute entry is returned.

      value = maxval(abs(x))
   end function maximum_norm_matrix

   pure real(dp) function maximum_norm_vector(x) result(value)
      real(dp), intent(in) :: x(:) !! Vector whose largest absolute entry is returned.

      value = maxval(abs(x))
   end function maximum_norm_vector

   pure real(dp) function one_norm_matrix(x) result(value)
      real(dp), intent(in) :: x(:, :) !! Matrix whose induced one-norm is computed.
      integer :: j

      value = 0.0_dp
      do j = 1, size(x, 2)
         value = max(value, sum(abs(x(:, j))))
      end do
   end function one_norm_matrix

   pure real(dp) function one_norm_vector(x) result(value)
      real(dp), intent(in) :: x(:) !! Vector interpreted as one column, matching upstream behavior.

      value = sum(abs(x))
   end function one_norm_vector

   pure real(dp) function inf_norm_matrix(x) result(value)
      real(dp), intent(in) :: x(:, :) !! Matrix whose induced infinity-norm is computed.
      integer :: i

      value = 0.0_dp
      do i = 1, size(x, 1)
         value = max(value, sum(abs(x(i, :))))
      end do
   end function inf_norm_matrix

   pure real(dp) function inf_norm_vector(x) result(value)
      real(dp), intent(in) :: x(:) !! Vector interpreted as one column, matching upstream behavior.

      value = maxval(abs(x))
   end function inf_norm_vector

end module matrixcalc_core
