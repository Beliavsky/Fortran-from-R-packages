module psych_matrix
   use r_kinds, only : dp
   use r_linalg, only : symmetric_eigen
   implicit none
   private
   public :: vech_values, matrix_addition_operator, matrix_square_root, lower_upper_merge
   public :: super_matrix_block, super_matrix_pair, missing_matmul, sort_matrix_by_loadings
contains
   pure subroutine vech_values(matrix, values, lower, diagonal)
      real(dp), intent(in) :: matrix(:, :) !! Square matrix whose selected triangle is vectorized in Fortran column order.
      real(dp), allocatable, intent(out) :: values(:) !! Selected strict or inclusive triangular entries.
      logical, intent(in), optional :: lower !! Select lower triangle when true; upper triangle otherwise.
      logical, intent(in), optional :: diagonal !! Include diagonal entries when true.
      logical :: use_lower, include_diagonal
      integer :: i, j, k, n, count_values
      use_lower = .true.
      include_diagonal = .false.
      if (present(lower)) use_lower = lower
      if (present(diagonal)) include_diagonal = diagonal
      n = size(matrix, 1)
      count_values = 0
      do j = 1, n
         do i = 1, n
            if (selected(i, j, use_lower, include_diagonal)) count_values = count_values + 1
         end do
      end do
      allocate(values(count_values))
      k = 0
      do j = 1, n
         do i = 1, n
            if (selected(i, j, use_lower, include_diagonal)) then
               k = k + 1
               values(k) = matrix(i, j)
            end if
         end do
      end do
   end subroutine vech_values

   pure subroutine matrix_addition_operator(x, y, result)
      real(dp), intent(in) :: x(:, :) !! Left matrix; its row sums form the row contribution.
      real(dp), intent(in) :: y(:, :) !! Right matrix; its column sums form the column contribution.
      real(dp), allocatable, intent(out) :: result(:, :) !! Outer addition of row sums of x and column sums of y.
      real(dp), allocatable :: row_sum(:), column_sum(:)
      integer :: i, j
      allocate(row_sum(size(x, 1)), column_sum(size(y, 2)))
      row_sum = sum(x, dim=2)
      column_sum = sum(y, dim=1)
      allocate(result(size(x, 1), size(y, 2)))
      do j = 1, size(y, 2)
         do i = 1, size(x, 1)
            result(i, j) = row_sum(i) + column_sum(j)
         end do
      end do
   end subroutine matrix_addition_operator

   subroutine matrix_square_root(matrix, root, info)
      real(dp), intent(in) :: matrix(:, :) !! Symmetric positive-semidefinite matrix.
      real(dp), allocatable, intent(out) :: root(:, :) !! Symmetric principal square root of matrix.
      integer, intent(out) :: info !! Zero on success or the eigensolver status.
      real(dp), allocatable :: values(:), vectors(:, :), scaled(:, :)
      integer :: j, n
      call symmetric_eigen(matrix, values, vectors, info, descending=.true.)
      if (info /= 0) then
         allocate(root(0, 0))
         return
      end if
      n = size(values)
      allocate(scaled(n, n))
      scaled = vectors
      do j = 1, n
         scaled(:, j) = scaled(:, j) * sqrt(max(values(j), 0.0_dp))
      end do
      root = matmul(scaled, transpose(vectors))
   end subroutine matrix_square_root

   pure subroutine lower_upper_merge(lower_matrix, upper_matrix, result, difference)
      real(dp), intent(in) :: lower_matrix(:, :) !! Matrix supplying lower-triangular values and the diagonal.
      real(dp), intent(in), optional :: upper_matrix(:, :) !! Optional matrix supplying upper-triangular values.
      real(dp), allocatable, intent(out) :: result(:, :) !! Matrix combining the requested triangular pieces.
      logical, intent(in), optional :: difference !! Use lower minus upper for the upper contribution when true.
      integer :: i, j, n
      logical :: use_difference
      n = size(lower_matrix, 1)
      allocate(result(n, n))
      use_difference = .false.
      if (present(difference)) use_difference = difference
      if (.not. present(upper_matrix)) then
         result = lower_matrix
         do j = 1, n
            do i = 1, j - 1
               result(i, j) = lower_matrix(j, i)
            end do
         end do
      else
         result = lower_matrix
         do j = 1, n
            do i = 1, j - 1
               if (use_difference) then
                  result(i, j) = lower_matrix(j, i) - upper_matrix(j, i)
               else
                  result(i, j) = upper_matrix(j, i)
               end if
            end do
         end do
      end if
   end subroutine lower_upper_merge

   pure subroutine super_matrix_block(blocks, result)
      real(dp), intent(in) :: blocks(:, :, :, :) !! Grid of equal-size matrix blocks: row, column, block-row, block-column.
      real(dp), allocatable, intent(out) :: result(:, :) !! Large matrix obtained by tiling the supplied blocks.
      integer :: br, bc, nr, nc, i0, j0
      nr = size(blocks, 1)
      nc = size(blocks, 2)
      allocate(result(nr * size(blocks, 3), nc * size(blocks, 4)))
      do bc = 1, size(blocks, 4)
         do br = 1, size(blocks, 3)
            i0 = (br - 1) * nr
            j0 = (bc - 1) * nc
            result(i0 + 1:i0 + nr, j0 + 1:j0 + nc) = blocks(:, :, br, bc)
         end do
      end do
   end subroutine super_matrix_block


   pure subroutine super_matrix_pair(x, y, result)
      real(dp), intent(in) :: x(:, :) !! First matrix placed in the upper-left block.
      real(dp), intent(in) :: y(:, :) !! Second matrix placed in the lower-right block.
      real(dp), allocatable, intent(out) :: result(:, :) !! Block-diagonal supermatrix with zero off-diagonal blocks.
      allocate(result(size(x, 1) + size(y, 1), size(x, 2) + size(y, 2)))
      result = 0.0_dp
      result(1:size(x, 1), 1:size(x, 2)) = x
      result(size(x, 1) + 1:, size(x, 2) + 1:) = y
   end subroutine super_matrix_pair

   pure subroutine missing_matmul(x, y, result)
      use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
      real(dp), intent(in) :: x(:, :) !! Left matrix; IEEE NaN products are omitted from each dot product.
      real(dp), intent(in) :: y(:, :) !! Right matrix; IEEE NaN products are omitted from each dot product.
      real(dp), allocatable, intent(out) :: result(:, :) !! Missing-tolerant matrix product using finite element products only.
      real(dp) :: product
      integer :: i, j, k
      allocate(result(size(x, 1), size(y, 2)))
      result = 0.0_dp
      do j = 1, size(y, 2)
         do i = 1, size(x, 1)
            do k = 1, size(x, 2)
               if (ieee_is_nan(x(i, k)) .or. ieee_is_nan(y(k, j))) cycle
               product = x(i, k) * y(k, j)
               result(i, j) = result(i, j) + product
            end do
         end do
      end do
   end subroutine missing_matmul

   pure subroutine sort_matrix_by_loadings(matrix, loadings, sorted_matrix, order)
      real(dp), intent(in) :: matrix(:, :) !! Square item matrix to reorder.
      real(dp), intent(in) :: loadings(:, :) !! Item-by-factor loading matrix used to define clusters and within-cluster order.
      real(dp), allocatable, intent(out) :: sorted_matrix(:, :) !! Matrix reordered by dominant factor and loading magnitude.
      integer, allocatable, intent(out) :: order(:) !! One-based item order applied to rows and columns.
      integer, allocatable :: cluster(:)
      integer :: i, j, k, n, best_factor, tmp_i
      real(dp) :: best_value
      n = size(loadings, 1)
      allocate(cluster(n), order(n))
      do i = 1, n
         best_factor = 1
         best_value = abs(loadings(i, 1))
         do j = 2, size(loadings, 2)
            if (abs(loadings(i, j)) > best_value) then
               best_value = abs(loadings(i, j))
               best_factor = j
            end if
         end do
         cluster(i) = best_factor
         order(i) = i
      end do
      do i = 2, n
         tmp_i = order(i)
         j = i - 1
         do while (j >= 1)
            if (cluster(order(j)) < cluster(tmp_i)) exit
            if (cluster(order(j)) == cluster(tmp_i)) then
               k = cluster(tmp_i)
               if (abs(loadings(order(j), k)) >= abs(loadings(tmp_i, k))) exit
            end if
            order(j + 1) = order(j)
            j = j - 1
         end do
         order(j + 1) = tmp_i
      end do
      allocate(sorted_matrix(n, n))
      do j = 1, n
         do i = 1, n
            sorted_matrix(i, j) = matrix(order(i), order(j))
         end do
      end do
   end subroutine sort_matrix_by_loadings

   pure logical function selected(i, j, lower, diagonal) result(use_entry)
      integer, intent(in) :: i !! Row index under consideration.
      integer, intent(in) :: j !! Column index under consideration.
      logical, intent(in) :: lower !! Whether lower-triangular entries are requested.
      logical, intent(in) :: diagonal !! Whether diagonal entries are included.
      if (lower) then
         use_entry = i > j .or. (diagonal .and. i == j)
      else
         use_entry = i < j .or. (diagonal .and. i == j)
      end if
   end function selected
end module psych_matrix
