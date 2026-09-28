! SPDX-License-Identifier: GPL-3.0-or-later
! Sparse direct factorization helpers for sanic.
module sanic_sparse_direct
   use sanic_kinds, only: dp
   use sanic_sparse, only: sparse_matrix, canonical_sparse_entries
   implicit none
   private

   type :: sparse_row
      integer, allocatable :: col(:)
      real(dp), allocatable :: value(:)
   end type sparse_row

   type :: givens_rotation
      integer :: first_row = 0
      integer :: second_row = 0
      real(dp) :: cosine = 1.0_dp
      real(dp) :: sine = 0.0_dp
   end type givens_rotation

   type :: index_set
      integer, allocatable :: item(:)
   end type index_set

   public :: sparse_chol_solve_core
   public :: sparse_lu_solve_core
   public :: sparse_qr_solve_core
   public :: sparse_column_permutation

contains

   pure subroutine sparse_lu_solve_core(a, b, x, info, ordering)
      !! Solve a sparse square system by fill-permitting sparse LU with partial row pivoting.
      type(sparse_matrix), intent(in) :: a !! Sparse square coefficient matrix in coordinate storage.
      real(dp), intent(in) :: b(:, :) !! Dense RHS matrix with a%nrow rows.
      real(dp), allocatable, intent(out) :: x(:, :) !! Solution matrix in the original column ordering.
      integer, intent(out) :: info !! Status: zero success, one invalid input, two numerical singularity.
      integer, intent(in) :: ordering !! 0 symmetric-degree, 1 column-degree, 2 natural ordering selector.
      type(sparse_row), allocatable :: rows(:)
      type(sparse_row) :: temporary_row
      real(dp), allocatable :: rhs(:, :), solution_perm(:, :), temporary_rhs(:)
      integer, allocatable :: inverse_perm(:), perm(:)
      real(dp) :: factor_value, pivot_value, scale, threshold, candidate
      integer :: i, j, k, p, pivot_row

      info = 0
      if (a%nrow /= a%ncol .or. size(b, 1) /= a%nrow .or. a%nrow < 1) then
         allocate(x(0, 0))
         info = 1
         return
      end if
      call sparse_column_permutation(a, ordering, perm, inverse_perm)
      call build_permuted_rows(a, inverse_perm, .false., rows)
      rhs = b
      scale = sparse_value_scale(a)
      threshold = epsilon(1.0_dp)*scale*real(max(1, a%ncol), dp)

      do k = 1, a%ncol
         pivot_row = k
         pivot_value = abs(row_get(rows(k), k))
         do i = k + 1, a%nrow
            candidate = abs(row_get(rows(i), k))
            if (candidate > pivot_value) then
               pivot_value = candidate
               pivot_row = i
            end if
         end do
         if (pivot_value <= threshold) then
            allocate(x(0, 0))
            info = 2
            return
         end if
         if (pivot_row /= k) then
            temporary_row = rows(k)
            rows(k) = rows(pivot_row)
            rows(pivot_row) = temporary_row
            temporary_rhs = rhs(k, :)
            rhs(k, :) = rhs(pivot_row, :)
            rhs(pivot_row, :) = temporary_rhs
         end if

         pivot_value = row_get(rows(k), k)
         do i = k + 1, a%nrow
            candidate = row_get(rows(i), k)
            if (abs(candidate) <= tiny(1.0_dp)) cycle
            factor_value = candidate/pivot_value
            call row_set(rows(i), k, factor_value)
            if (allocated(rows(k)%col)) then
               do p = 1, size(rows(k)%col)
                  j = rows(k)%col(p)
                  if (j <= k) cycle
                  candidate = row_get(rows(i), j) - factor_value*rows(k)%value(p)
                  call row_set(rows(i), j, candidate)
               end do
            end if
         end do
      end do

      allocate(solution_perm(a%ncol, size(b, 2)))
      solution_perm = rhs
      do i = 1, a%ncol
         if (allocated(rows(i)%col)) then
            do p = 1, size(rows(i)%col)
               j = rows(i)%col(p)
               if (j >= i) exit
               solution_perm(i, :) = solution_perm(i, :) - rows(i)%value(p)*solution_perm(j, :)
            end do
         end if
      end do
      do i = a%ncol, 1, -1
         pivot_value = row_get(rows(i), i)
         if (abs(pivot_value) <= threshold) then
            deallocate(solution_perm)
            allocate(x(0, 0))
            info = 2
            return
         end if
         if (allocated(rows(i)%col)) then
            do p = 1, size(rows(i)%col)
               j = rows(i)%col(p)
               if (j > i) solution_perm(i, :) = solution_perm(i, :) - rows(i)%value(p)*solution_perm(j, :)
            end do
         end if
         solution_perm(i, :) = solution_perm(i, :)/pivot_value
      end do

      allocate(x(a%ncol, size(b, 2)))
      x = 0.0_dp
      do i = 1, a%ncol
         x(perm(i), :) = solution_perm(i, :)
      end do
   end subroutine sparse_lu_solve_core

   pure subroutine sparse_chol_solve_core(a, b, x, info, llt, ordering)
      !! Solve a sparse symmetric system by fill-permitting LLT or LDLT factorization.
      type(sparse_matrix), intent(in) :: a !! Sparse symmetric square coefficient matrix in coordinate storage.
      real(dp), intent(in) :: b(:, :) !! Dense RHS matrix with a%nrow rows.
      real(dp), allocatable, intent(out) :: x(:, :) !! Solution matrix in the original ordering.
      integer, intent(out) :: info !! Status: zero success, one invalid input, two failed numerical factorization.
      logical, intent(in) :: llt !! True selects positive-definite LLT; false selects LDLT.
      integer, intent(in) :: ordering !! 0/1 degree-based symmetric ordering, 2 natural ordering.
      type(sparse_row), allocatable :: a_rows(:), l_rows(:)
      real(dp), allocatable :: d(:), rhs(:, :), work(:, :)
      integer, allocatable :: inverse_perm(:), perm(:)
      real(dp) :: diagonal, lij, scale, sum_value, threshold
      integer :: i, j, k, p

      info = 0
      if (a%nrow /= a%ncol .or. size(b, 1) /= a%nrow .or. a%nrow < 1) then
         allocate(x(0, 0))
         info = 1
         return
      end if
      call sparse_symmetric_permutation(a, ordering, perm, inverse_perm)
      call build_permuted_rows(a, inverse_perm, .true., a_rows)
      allocate(l_rows(a%ncol), d(a%ncol))
      d = 0.0_dp
      scale = sparse_value_scale(a)
      threshold = epsilon(1.0_dp)*scale*real(max(1, a%ncol), dp)

      do i = 1, a%ncol
         do k = 1, i - 1
            sum_value = row_get(a_rows(i), k)
            if (allocated(l_rows(i)%col)) then
               do p = 1, size(l_rows(i)%col)
                  j = l_rows(i)%col(p)
                  if (j >= k) exit
                  if (llt) then
                     sum_value = sum_value - l_rows(i)%value(p)*row_get(l_rows(k), j)
                  else
                     sum_value = sum_value - l_rows(i)%value(p)*row_get(l_rows(k), j)*d(j)
                  end if
               end do
            end if
            if (llt) then
               diagonal = row_get(l_rows(k), k)
               if (diagonal <= threshold) then
                  allocate(x(0, 0))
                  info = 2
                  return
               end if
               lij = sum_value/diagonal
            else
               if (abs(d(k)) <= threshold) then
                  allocate(x(0, 0))
                  info = 2
                  return
               end if
               lij = sum_value/d(k)
            end if
            if (abs(lij) > tiny(1.0_dp)) call row_set(l_rows(i), k, lij)
         end do

         diagonal = row_get(a_rows(i), i)
         if (allocated(l_rows(i)%col)) then
            do p = 1, size(l_rows(i)%col)
               j = l_rows(i)%col(p)
               if (j >= i) exit
               if (llt) then
                  diagonal = diagonal - l_rows(i)%value(p)*l_rows(i)%value(p)
               else
                  diagonal = diagonal - l_rows(i)%value(p)*l_rows(i)%value(p)*d(j)
               end if
            end do
         end if
         if (llt) then
            if (diagonal <= threshold) then
               allocate(x(0, 0))
               info = 2
               return
            end if
            call row_set(l_rows(i), i, sqrt(diagonal))
         else
            d(i) = diagonal
            if (abs(d(i)) <= threshold) then
               allocate(x(0, 0))
               info = 2
               return
            end if
            call row_set(l_rows(i), i, 1.0_dp)
         end if
      end do

      allocate(rhs(a%nrow, size(b, 2)))
      do i = 1, a%nrow
         rhs(i, :) = b(perm(i), :)
      end do
      allocate(work(a%ncol, size(b, 2)))
      work = rhs
      do i = 1, a%ncol
         if (allocated(l_rows(i)%col)) then
            do p = 1, size(l_rows(i)%col)
               j = l_rows(i)%col(p)
               if (j >= i) exit
               work(i, :) = work(i, :) - l_rows(i)%value(p)*work(j, :)
            end do
         end if
         if (llt) work(i, :) = work(i, :)/row_get(l_rows(i), i)
      end do
      if (.not. llt) then
         do i = 1, a%ncol
            work(i, :) = work(i, :)/d(i)
         end do
      end if
      do i = a%ncol, 1, -1
         do j = i + 1, a%ncol
            lij = row_get(l_rows(j), i)
            if (abs(lij) > tiny(1.0_dp)) work(i, :) = work(i, :) - lij*work(j, :)
         end do
         if (llt) work(i, :) = work(i, :)/row_get(l_rows(i), i)
      end do

      allocate(x(a%ncol, size(b, 2)))
      x = 0.0_dp
      do i = 1, a%ncol
         x(perm(i), :) = work(i, :)
      end do
   end subroutine sparse_chol_solve_core

   pure subroutine sparse_qr_solve_core(a, b, x, info, ordering)
      !! Solve a sparse least-squares or underdetermined problem without densifying the coefficient matrix.
      type(sparse_matrix), intent(in) :: a !! Sparse m-by-n coefficient matrix in coordinate storage.
      real(dp), intent(in) :: b(:, :) !! Dense RHS matrix with a%nrow rows.
      real(dp), allocatable, intent(out) :: x(:, :) !! Eigen-style basic solution in original column order.
      integer, intent(out) :: info !! Status: zero success, one invalid input, two numerical rank failure.
      integer, intent(in) :: ordering !! 0 dynamic symmetric degree, 1 dynamic column-intersection degree, 2 natural order.
      type(sparse_row), allocatable :: rows(:)
      type(givens_rotation), allocatable :: rotations(:)
      real(dp), allocatable :: rhs(:, :), solution_perm(:, :)
      integer, allocatable :: inverse_perm(:), perm(:)
      real(dp) :: diagonal, pivot_threshold, threshold
      integer :: i, j, nrhs, rank_value

      info = 0
      if (a%nrow < 1 .or. a%ncol < 1 .or. size(b, 1) /= a%nrow) then
         allocate(x(0, 0))
         info = 1
         return
      end if
      nrhs = size(b, 2)
      call sparse_column_permutation(a, ordering, perm, inverse_perm)
      call build_permuted_rows(a, inverse_perm, .false., rows)
      rhs = b
      threshold = epsilon(1.0_dp)*sparse_value_scale(a)*real(max(a%nrow, a%ncol), dp)
      pivot_threshold = sparse_qr_default_threshold(rows, a%nrow, a%ncol)
      call sparse_givens_factor(rows, rhs, perm, a%ncol, .true., threshold, pivot_threshold, rotations, rank_value)

      allocate(solution_perm(a%ncol, nrhs))
      solution_perm = 0.0_dp
      if (rank_value > 0) then
         do i = rank_value, 1, -1
            solution_perm(i, :) = rhs(i, :)
            do j = i + 1, rank_value
               solution_perm(i, :) = solution_perm(i, :) - row_get(rows(i), j)*solution_perm(j, :)
            end do
            diagonal = row_get(rows(i), i)
            if (abs(diagonal) <= threshold) then
               deallocate(solution_perm)
               allocate(x(0, 0))
               info = 2
               return
            end if
            solution_perm(i, :) = solution_perm(i, :)/diagonal
         end do
      end if

      allocate(x(a%ncol, nrhs))
      x = 0.0_dp
      do i = 1, a%ncol
         x(perm(i), :) = solution_perm(i, :)
      end do
   end subroutine sparse_qr_solve_core

   pure subroutine sparse_column_permutation(a, ordering, perm, inverse_perm)
      !! Build natural or sparse dynamic minimum-degree structural column ordering.
      type(sparse_matrix), intent(in) :: a !! Sparse matrix whose nonzero graph defines the ordering.
      integer, intent(in) :: ordering !! 0 symmetrized graph, 1 column-intersection graph, 2 natural order.
      integer, allocatable, intent(out) :: perm(:) !! New-index to old-index column permutation.
      integer, allocatable, intent(out) :: inverse_perm(:) !! Old-index to new-index inverse permutation.
      type(index_set), allocatable :: graph(:)
      integer :: i, j, k, mode

      allocate(perm(a%ncol), inverse_perm(a%ncol))
      perm = [(i, i=1, a%ncol)]
      mode = ordering
      if (mode < 0 .or. mode > 2) mode = 1
      if (mode == 2 .or. a%ncol <= 1) then
         inverse_perm = perm
         return
      end if

      allocate(graph(a%ncol))
      if (mode == 0 .and. a%nrow == a%ncol) then
         if (allocated(a%value)) then
            do k = 1, size(a%value)
               i = a%row(k)
               j = a%col(k)
               if (i < 1 .or. i > a%ncol .or. j < 1 .or. j > a%ncol .or. i == j) cycle
               call graph_add_edge(graph, i, j)
            end do
         end if
      else
         call build_column_intersection_graph(a, graph)
      end if
      call exact_minimum_degree_order(graph, perm)
      do i = 1, a%ncol
         inverse_perm(perm(i)) = i
      end do
   end subroutine sparse_column_permutation

   pure subroutine sparse_symmetric_permutation(a, ordering, perm, inverse_perm)
      !! Build the symmetric ordering used by sparse Cholesky; ordering=1 falls back to the default degree ordering.
      type(sparse_matrix), intent(in) :: a !! Sparse symmetric matrix whose graph defines the ordering.
      integer, intent(in) :: ordering !! 0 default degree ordering, 2 natural; ordering=1 follows upstream default fallback.
      integer, allocatable, intent(out) :: perm(:) !! New-index to old-index symmetric permutation.
      integer, allocatable, intent(out) :: inverse_perm(:) !! Old-index to new-index inverse permutation.
      integer :: mode

      mode = ordering
      if (mode == 2) then
         call sparse_column_permutation(a, 2, perm, inverse_perm)
      else
         call sparse_column_permutation(a, 0, perm, inverse_perm)
      end if
   end subroutine sparse_symmetric_permutation

   pure subroutine build_permuted_rows(a, inverse_perm, symmetric_rows, rows)
      !! Convert coordinate storage into sorted dynamic sparse rows after column or symmetric permutation.
      type(sparse_matrix), intent(in) :: a !! Coordinate sparse matrix to permute.
      integer, intent(in) :: inverse_perm(:) !! Old-index to new-index permutation mapping.
      logical, intent(in) :: symmetric_rows !! True applies the same permutation to rows and columns.
      type(sparse_row), allocatable, intent(out) :: rows(:) !! Sorted row representation with duplicates accumulated.
      integer :: i, j, k, new_i, new_j
      real(dp) :: accumulated

      allocate(rows(a%nrow))
      if (.not. allocated(a%value)) return
      do k = 1, size(a%value)
         i = a%row(k)
         j = a%col(k)
         if (i < 1 .or. i > a%nrow .or. j < 1 .or. j > a%ncol) cycle
         if (symmetric_rows) then
            new_i = inverse_perm(i)
         else
            new_i = i
         end if
         new_j = inverse_perm(j)
         accumulated = row_get(rows(new_i), new_j) + a%value(k)
         call row_set(rows(new_i), new_j, accumulated)
      end do
   end subroutine build_permuted_rows

   pure subroutine build_column_intersection_graph(a, graph)
      !! Build sparse column intersections from canonical row-sorted entries without repeated triplet scans.
      type(sparse_matrix), intent(in) :: a !! Sparse matrix whose row structure defines column intersections.
      type(index_set), intent(inout) :: graph(:) !! Sparse undirected adjacency lists of length a%ncol.
      integer, allocatable :: col_index(:), row_index(:)
      real(dp), allocatable :: values(:)
      integer :: first, i, j, last, nnz

      call canonical_sparse_entries(a, row_index, col_index, values)
      nnz = size(values)
      first = 1
      do while (first <= nnz)
         last = first
         do while (last < nnz)
            if (row_index(last + 1) /= row_index(first)) exit
            last = last + 1
         end do
         do i = first, last - 1
            do j = i + 1, last
               if (col_index(i) /= col_index(j)) call graph_add_edge(graph, col_index(i), col_index(j))
            end do
         end do
         first = last + 1
      end do
   end subroutine build_column_intersection_graph

   pure subroutine exact_minimum_degree_order(graph, perm)
      !! Greedily eliminate the current minimum-degree vertex while adding sparse fill-clique edges.
      type(index_set), intent(inout) :: graph(:) !! Sparse undirected adjacency lists; fill edges are added during elimination.
      integer, intent(out) :: perm(:) !! New-index to old-index elimination ordering.
      logical, allocatable :: active(:)
      integer, allocatable :: neighbors(:)
      integer :: best, best_degree, degree_value, i, j, k, nneighbor, n, p

      n = size(graph)
      allocate(active(n), neighbors(n))
      active = .true.
      do k = 1, n
         best = 0
         best_degree = huge(1)
         do i = 1, n
            if (.not. active(i)) cycle
            degree_value = 0
            if (allocated(graph(i)%item)) then
               do p = 1, size(graph(i)%item)
                  if (active(graph(i)%item(p))) degree_value = degree_value + 1
               end do
            end if
            if (degree_value < best_degree .or. (degree_value == best_degree .and. (best == 0 .or. i < best))) then
               best = i
               best_degree = degree_value
            end if
         end do
         perm(k) = best
         nneighbor = 0
         if (allocated(graph(best)%item)) then
            do p = 1, size(graph(best)%item)
               i = graph(best)%item(p)
               if (.not. active(i)) cycle
               nneighbor = nneighbor + 1
               neighbors(nneighbor) = i
            end do
         end if
         do i = 1, nneighbor
            do j = i + 1, nneighbor
               call graph_add_edge(graph, neighbors(i), neighbors(j))
            end do
         end do
         active(best) = .false.
      end do
   end subroutine exact_minimum_degree_order

   pure subroutine graph_add_edge(graph, first, second)
      !! Add one undirected edge to sparse sorted adjacency lists when it is not already present.
      type(index_set), intent(inout) :: graph(:) !! Sparse adjacency lists to modify.
      integer, intent(in) :: first !! First one-based vertex index.
      integer, intent(in) :: second !! Second one-based vertex index.

      if (first == second) return
      call index_set_add(graph(first), second)
      call index_set_add(graph(second), first)
   end subroutine graph_add_edge

   pure subroutine index_set_add(set_value, new_item)
      !! Insert one integer into a sorted unique sparse index set.
      type(index_set), intent(inout) :: set_value !! Sorted unique adjacency set to modify.
      integer, intent(in) :: new_item !! Vertex index to add when absent.
      integer, allocatable :: extended(:)
      integer :: insert_at, n

      if (.not. allocated(set_value%item)) then
         allocate(set_value%item(1))
         set_value%item(1) = new_item
         return
      end if
      n = size(set_value%item)
      insert_at = 1
      do while (insert_at <= n)
         if (set_value%item(insert_at) >= new_item) exit
         insert_at = insert_at + 1
      end do
      if (insert_at <= n) then
         if (set_value%item(insert_at) == new_item) return
      end if
      allocate(extended(n + 1))
      if (insert_at > 1) extended(1:insert_at - 1) = set_value%item(1:insert_at - 1)
      extended(insert_at) = new_item
      if (insert_at <= n) extended(insert_at + 1:n + 1) = set_value%item(insert_at:n)
      call move_alloc(extended, set_value%item)
   end subroutine index_set_add



   pure subroutine sparse_givens_factor(rows, rhs, perm, ncol_active, pivot_columns, threshold, &
      pivot_threshold, rotations, rank_value)
      !! Reduce sparse rows to upper trapezoidal form with Givens rotations and optional numerical column pivoting.
      type(sparse_row), intent(inout) :: rows(:) !! Sparse matrix rows overwritten by the R factor.
      real(dp), intent(inout) :: rhs(:, :) !! Dense RHS transformed by the same left orthogonal rotations.
      integer, intent(inout) :: perm(:) !! Current new-index to original-column mapping, updated by column pivots.
      integer, intent(in) :: ncol_active !! Number of active matrix columns represented in rows.
      logical, intent(in) :: pivot_columns !! True enables rank-revealing numerical column pivoting.
      real(dp), intent(in) :: threshold !! Scale-aware pruning threshold for sparse arithmetic.
      real(dp), intent(in) :: pivot_threshold !! Eigen-style column threshold below which a direction is treated as rank deficient.
      type(givens_rotation), allocatable, intent(out) :: rotations(:) !! Ordered Givens rotations defining Q transpose.
      integer, intent(out) :: rank_value !! Number of numerically independent leading columns.
      real(dp), allocatable :: tmp_rhs(:)
      real(dp) :: a_value, b_value, cosine, current_norm, norm_value, replacement_norm, sine
      integer :: best, i, k, limit_value, old_perm

      allocate(rotations(0))
      rank_value = 0
      limit_value = min(size(rows), ncol_active)
      do k = 1, limit_value
         current_norm = sqrt(max(0.0_dp, trailing_column_norm2(rows, k, k)))
         if (pivot_columns .and. current_norm <= pivot_threshold) then
            best = trailing_best_column(rows, k, ncol_active)
            replacement_norm = sqrt(max(0.0_dp, trailing_column_norm2(rows, k, best)))
            if (best /= k .and. replacement_norm > pivot_threshold) then
               call swap_sparse_columns(rows, k, best)
               old_perm = perm(k)
               perm(k) = perm(best)
               perm(best) = old_perm
               current_norm = replacement_norm
            end if
         end if
         if (current_norm <= pivot_threshold) exit

         do i = k + 1, size(rows)
            b_value = row_get(rows(i), k)
            if (abs(b_value) <= threshold) cycle
            a_value = row_get(rows(k), k)
            norm_value = hypot(a_value, b_value)
            if (norm_value <= threshold) cycle
            cosine = a_value/norm_value
            sine = b_value/norm_value
            call rotate_sparse_rows(rows(k), rows(i), cosine, sine, threshold)
            if (size(rhs, 2) > 0) then
               tmp_rhs = cosine*rhs(k, :) + sine*rhs(i, :)
               rhs(i, :) = -sine*rhs(k, :) + cosine*rhs(i, :)
               rhs(k, :) = tmp_rhs
            end if
            call append_rotation(rotations, k, i, cosine, sine)
         end do
         if (abs(row_get(rows(k), k)) > pivot_threshold) rank_value = rank_value + 1
      end do
   end subroutine sparse_givens_factor


   pure real(dp) function sparse_qr_default_threshold(rows, nrow, ncol) result(value)
      !! Return Eigen SparseQR's default rank threshold 20*(m+n)*epsilon*max column 2-norm.
      type(sparse_row), intent(in) :: rows(:) !! Sparse rows before numerical QR elimination.
      integer, intent(in) :: nrow !! Number of rows in the represented matrix.
      integer, intent(in) :: ncol !! Number of columns in the represented matrix.
      real(dp) :: max_norm
      integer :: j

      max_norm = 0.0_dp
      do j = 1, ncol
         max_norm = max(max_norm, sqrt(max(0.0_dp, trailing_column_norm2(rows, 1, j))))
      end do
      value = 20.0_dp*real(nrow + ncol, dp)*epsilon(1.0_dp)*max_norm
      if (value <= 0.0_dp) value = 20.0_dp*real(nrow + ncol, dp)*epsilon(1.0_dp)
   end function sparse_qr_default_threshold

   pure integer function trailing_best_column(rows, start_index, ncol) result(best)
      !! Select the remaining column with greatest trailing Euclidean norm for numerical rank revealing.
      type(sparse_row), intent(in) :: rows(:) !! Current sparse matrix rows.
      integer, intent(in) :: start_index !! First active row and column index.
      integer, intent(in) :: ncol !! Number of matrix columns represented by the rows.
      real(dp) :: best_norm, current_norm
      integer :: j

      best = start_index
      best_norm = -1.0_dp
      do j = start_index, ncol
         current_norm = trailing_column_norm2(rows, start_index, j)
         if (current_norm > best_norm) then
            best_norm = current_norm
            best = j
         end if
      end do
   end function trailing_best_column

   pure real(dp) function trailing_column_norm2(rows, start_row, target_col) result(value)
      !! Return squared 2-norm of one sparse column below and including start_row.
      type(sparse_row), intent(in) :: rows(:) !! Sparse matrix rows to inspect.
      integer, intent(in) :: start_row !! First row included in the norm.
      integer, intent(in) :: target_col !! Column whose trailing norm is requested.
      real(dp) :: entry
      integer :: i

      value = 0.0_dp
      do i = start_row, size(rows)
         entry = row_get(rows(i), target_col)
         value = value + entry*entry
      end do
   end function trailing_column_norm2

   pure subroutine swap_sparse_columns(rows, first_col, second_col)
      !! Exchange two columns in a dynamic sparse-row matrix without forming dense storage.
      type(sparse_row), intent(inout) :: rows(:) !! Sparse rows whose columns are exchanged.
      integer, intent(in) :: first_col !! First one-based column index.
      integer, intent(in) :: second_col !! Second one-based column index.
      real(dp) :: first_value, second_value
      integer :: i

      if (first_col == second_col) return
      do i = 1, size(rows)
         first_value = row_get(rows(i), first_col)
         second_value = row_get(rows(i), second_col)
         call row_assign(rows(i), first_col, second_value)
         call row_assign(rows(i), second_col, first_value)
      end do
   end subroutine swap_sparse_columns

   pure subroutine rotate_sparse_rows(first, second, cosine, sine, threshold)
      !! Apply one real Givens rotation to two sparse rows using a sorted union of their column indices.
      type(sparse_row), intent(inout) :: first !! First sparse row, overwritten by c*first+s*second.
      type(sparse_row), intent(inout) :: second !! Second sparse row, overwritten by -s*first+c*second.
      real(dp), intent(in) :: cosine !! Givens cosine coefficient.
      real(dp), intent(in) :: sine !! Givens sine coefficient.
      real(dp), intent(in) :: threshold !! Absolute pruning threshold for generated roundoff entries.
      integer, allocatable :: col_first(:), col_second(:), union_col(:)
      real(dp), allocatable :: value_first(:), value_second(:), new_first(:), new_second(:)
      real(dp) :: left_value, right_value
      integer :: i, j, nfirst, nsecond, nunion

      call copy_row_arrays(first, col_first, value_first)
      call copy_row_arrays(second, col_second, value_second)
      nfirst = size(col_first)
      nsecond = size(col_second)
      allocate(union_col(nfirst + nsecond), new_first(nfirst + nsecond), new_second(nfirst + nsecond))
      i = 1
      j = 1
      nunion = 0
      do while (i <= nfirst .or. j <= nsecond)
         nunion = nunion + 1
         if (i > nfirst) then
            union_col(nunion) = col_second(j)
            left_value = 0.0_dp
            right_value = value_second(j)
            j = j + 1
         else if (j > nsecond) then
            union_col(nunion) = col_first(i)
            left_value = value_first(i)
            right_value = 0.0_dp
            i = i + 1
         else if (col_first(i) < col_second(j)) then
            union_col(nunion) = col_first(i)
            left_value = value_first(i)
            right_value = 0.0_dp
            i = i + 1
         else if (col_second(j) < col_first(i)) then
            union_col(nunion) = col_second(j)
            left_value = 0.0_dp
            right_value = value_second(j)
            j = j + 1
         else
            union_col(nunion) = col_first(i)
            left_value = value_first(i)
            right_value = value_second(j)
            i = i + 1
            j = j + 1
         end if
         new_first(nunion) = cosine*left_value + sine*right_value
         new_second(nunion) = -sine*left_value + cosine*right_value
      end do
      call replace_row_from_union(first, union_col(1:nunion), new_first(1:nunion), threshold)
      call replace_row_from_union(second, union_col(1:nunion), new_second(1:nunion), threshold)
   end subroutine rotate_sparse_rows

   pure subroutine copy_row_arrays(row, columns, values)
      !! Copy a sparse row into always-allocated arrays, using zero length for an empty row.
      type(sparse_row), intent(in) :: row !! Sparse row to copy.
      integer, allocatable, intent(out) :: columns(:) !! Sorted copied column indices.
      real(dp), allocatable, intent(out) :: values(:) !! Copied numerical values.

      if (allocated(row%col)) then
         columns = row%col
         values = row%value
      else
         allocate(columns(0), values(0))
      end if
   end subroutine copy_row_arrays

   pure subroutine replace_row_from_union(row, columns, values, threshold)
      !! Replace a sparse row from a sorted candidate union while pruning tiny entries.
      type(sparse_row), intent(inout) :: row !! Sparse row to replace.
      integer, intent(in) :: columns(:) !! Sorted candidate columns.
      real(dp), intent(in) :: values(:) !! Candidate values aligned with columns.
      real(dp), intent(in) :: threshold !! Entries with absolute value at or below threshold are omitted.
      integer, allocatable :: kept_col(:)
      real(dp), allocatable :: kept_value(:)
      integer :: count_value, i

      count_value = count(abs(values) > threshold)
      if (allocated(row%col)) deallocate(row%col)
      if (allocated(row%value)) deallocate(row%value)
      allocate(kept_col(count_value), kept_value(count_value))
      count_value = 0
      do i = 1, size(columns)
         if (abs(values(i)) <= threshold) cycle
         count_value = count_value + 1
         kept_col(count_value) = columns(i)
         kept_value(count_value) = values(i)
      end do
      call move_alloc(kept_col, row%col)
      call move_alloc(kept_value, row%value)
   end subroutine replace_row_from_union

   pure subroutine row_assign(row, target_col, new_value)
      !! Assign one row entry and remove the structural entry when the new value is exactly zero.
      type(sparse_row), intent(inout) :: row !! Sparse row to modify.
      integer, intent(in) :: target_col !! One-based column index to assign.
      real(dp), intent(in) :: new_value !! Replacement value; exact zero removes the entry.

      if (abs(new_value) <= tiny(1.0_dp)) then
         call row_remove(row, target_col)
      else
         call row_set(row, target_col, new_value)
      end if
   end subroutine row_assign

   pure subroutine row_remove(row, target_col)
      !! Remove one column entry from a sparse row when present.
      type(sparse_row), intent(inout) :: row !! Sparse row to modify.
      integer, intent(in) :: target_col !! One-based column index to remove.
      integer, allocatable :: new_col(:)
      real(dp), allocatable :: new_value(:)
      integer :: i, index_value, n

      if (.not. allocated(row%col)) return
      n = size(row%col)
      index_value = 0
      do i = 1, n
         if (row%col(i) == target_col) then
            index_value = i
            exit
         end if
      end do
      if (index_value == 0) return
      allocate(new_col(n - 1), new_value(n - 1))
      if (index_value > 1) then
         new_col(1:index_value - 1) = row%col(1:index_value - 1)
         new_value(1:index_value - 1) = row%value(1:index_value - 1)
      end if
      if (index_value < n) then
         new_col(index_value:n - 1) = row%col(index_value + 1:n)
         new_value(index_value:n - 1) = row%value(index_value + 1:n)
      end if
      call move_alloc(new_col, row%col)
      call move_alloc(new_value, row%value)
   end subroutine row_remove

   pure subroutine append_rotation(rotations, first_row, second_row, cosine, sine)
      !! Append one Givens rotation to a dynamically sized rotation sequence.
      type(givens_rotation), allocatable, intent(inout) :: rotations(:) !! Rotation sequence accumulated in factor order.
      integer, intent(in) :: first_row !! First rotated row index.
      integer, intent(in) :: second_row !! Second rotated row index.
      real(dp), intent(in) :: cosine !! Givens cosine coefficient.
      real(dp), intent(in) :: sine !! Givens sine coefficient.
      type(givens_rotation), allocatable :: extended(:)
      integer :: n

      n = size(rotations)
      allocate(extended(n + 1))
      if (n > 0) extended(1:n) = rotations
      extended(n + 1)%first_row = first_row
      extended(n + 1)%second_row = second_row
      extended(n + 1)%cosine = cosine
      extended(n + 1)%sine = sine
      call move_alloc(extended, rotations)
   end subroutine append_rotation



   pure real(dp) function sparse_value_scale(a) result(scale)
      !! Return max(1,max(abs(nonzero values))) for numerical pivot thresholds.
      type(sparse_matrix), intent(in) :: a !! Sparse matrix whose stored values set the numerical scale.

      scale = 1.0_dp
      if (allocated(a%value)) then
         if (size(a%value) > 0) scale = max(1.0_dp, maxval(abs(a%value)))
      end if
   end function sparse_value_scale

   pure real(dp) function row_get(row, target_col) result(value)
      !! Read one entry from a sorted sparse row, returning exact zero for an absent structural entry.
      type(sparse_row), intent(in) :: row !! Sorted sparse row to query.
      integer, intent(in) :: target_col !! One-based column index to query.
      integer :: left, middle, right

      value = 0.0_dp
      if (.not. allocated(row%col)) return
      left = 1
      right = size(row%col)
      do while (left <= right)
         middle = (left + right)/2
         if (row%col(middle) == target_col) then
            value = row%value(middle)
            return
         else if (row%col(middle) < target_col) then
            left = middle + 1
         else
            right = middle - 1
         end if
      end do
   end function row_get

   pure subroutine row_set(row, target_col, new_value)
      !! Insert or replace one entry while preserving sorted sparse-row column order.
      type(sparse_row), intent(inout) :: row !! Sparse row to modify.
      integer, intent(in) :: target_col !! One-based column index to insert or replace.
      real(dp), intent(in) :: new_value !! Numerical value to store, including explicit zeros when generated.
      integer, allocatable :: new_col(:)
      real(dp), allocatable :: new_values(:)
      integer :: insert_at, n

      if (.not. allocated(row%col)) then
         allocate(row%col(1), row%value(1))
         row%col(1) = target_col
         row%value(1) = new_value
         return
      end if
      n = size(row%col)
      insert_at = 1
      do while (insert_at <= n)
         if (row%col(insert_at) >= target_col) exit
         insert_at = insert_at + 1
      end do
      if (insert_at <= n) then
         if (row%col(insert_at) == target_col) then
            row%value(insert_at) = new_value
            return
         end if
      end if

      allocate(new_col(n + 1), new_values(n + 1))
      if (insert_at > 1) then
         new_col(1:insert_at - 1) = row%col(1:insert_at - 1)
         new_values(1:insert_at - 1) = row%value(1:insert_at - 1)
      end if
      new_col(insert_at) = target_col
      new_values(insert_at) = new_value
      if (insert_at <= n) then
         new_col(insert_at + 1:n + 1) = row%col(insert_at:n)
         new_values(insert_at + 1:n + 1) = row%value(insert_at:n)
      end if
      call move_alloc(new_col, row%col)
      call move_alloc(new_values, row%value)
   end subroutine row_set

end module sanic_sparse_direct
