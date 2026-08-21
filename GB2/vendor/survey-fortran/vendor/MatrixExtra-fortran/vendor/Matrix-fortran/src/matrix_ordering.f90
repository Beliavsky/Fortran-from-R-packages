! SPDX-License-Identifier: GPL-3.0-only
module matrix_ordering
   use matrix_sparse, only : csr_matrix, csr_transpose
   use matrix_status, only : matrix_success, matrix_err_shape
   implicit none
   private
   public :: reverse_cuthill_mckee, minimum_degree_ordering, column_degree_ordering

contains

   subroutine reverse_cuthill_mckee(a, perm, info)
      type(csr_matrix), intent(in) :: a
      integer, allocatable, intent(out) :: perm(:)
      integer, intent(out) :: info
      integer, allocatable :: degree(:), queue(:), nbr(:), order(:)
      logical, allocatable :: visited(:)
      integer :: n, i, k, qhead, qtail, count_order, start, nnbr, v, u
      if (a%nrow /= a%ncol) then
         allocate(perm(0))
         info = matrix_err_shape
         return
      end if
      n = a%nrow
      allocate(degree(n), queue(n), nbr(n), order(n))
      allocate(visited(n), source=.false.)
      do i = 1, n
         degree(i) = a%row_ptr(i + 1) - a%row_ptr(i)
      end do
      count_order = 0
      do while (count_order < n)
         start = 0
         do i = 1, n
            if (.not. visited(i)) then
               if (start == 0) then
                  start = i
               else if (degree(i) < degree(start)) then
                  start = i
               end if
            end if
         end do
         qhead = 1
         qtail = 1
         queue(1) = start
         visited(start) = .true.
         do while (qhead <= qtail)
            v = queue(qhead)
            qhead = qhead + 1
            count_order = count_order + 1
            order(count_order) = v
            nnbr = 0
            do k = a%row_ptr(v), a%row_ptr(v + 1) - 1
               u = a%col_ind(k)
               if (u /= v .and. .not. visited(u)) then
                  visited(u) = .true.
                  nnbr = nnbr + 1
                  nbr(nnbr) = u
               end if
            end do
            call sort_by_degree(nbr(:nnbr), degree)
            do i = 1, nnbr
               qtail = qtail + 1
               queue(qtail) = nbr(i)
            end do
         end do
      end do
      allocate(perm(n))
      do i = 1, n
         perm(i) = order(n - i + 1)
      end do
      info = matrix_success
   end subroutine reverse_cuthill_mckee

   subroutine sort_by_degree(nodes, degree)
      integer, intent(inout) :: nodes(:)
      integer, intent(in) :: degree(:)
      integer :: i, j, key
      do i = 2, size(nodes)
         key = nodes(i)
         j = i - 1
         do while (j >= 1)
            if (degree(nodes(j)) < degree(key)) exit
            if (degree(nodes(j)) == degree(key) .and. nodes(j) <= key) exit
            nodes(j + 1) = nodes(j)
            j = j - 1
         end do
         nodes(j + 1) = key
      end do
   end subroutine sort_by_degree

   subroutine minimum_degree_ordering(a, perm, info)
      type(csr_matrix), intent(in) :: a
      integer, allocatable, intent(out) :: perm(:)
      integer, intent(out) :: info
      logical, allocatable :: active(:), adjacency(:,:)
      integer :: n, i, j, k, step, v, u, w, best_degree, degree
      if (a%nrow /= a%ncol) then
         allocate(perm(0))
         info = matrix_err_shape
         return
      end if
      n = a%nrow
      allocate(adjacency(n, n), source=.false.)
      do i = 1, n
         do k = a%row_ptr(i), a%row_ptr(i + 1) - 1
            j = a%col_ind(k)
            if (i /= j) then
               adjacency(i, j) = .true.
               adjacency(j, i) = .true.
            end if
         end do
      end do
      allocate(active(n), source=.true.)
      allocate(perm(n))
      do step = 1, n
         v = 0
         best_degree = huge(1)
         do i = 1, n
            if (active(i)) then
               degree = count(adjacency(i, :) .and. active)
               if (degree < best_degree) then
                  best_degree = degree
                  v = i
               end if
            end if
         end do
         perm(step) = v
         do u = 1, n
            if (active(u) .and. adjacency(v, u)) then
               do w = u + 1, n
                  if (active(w) .and. adjacency(v, w)) then
                     adjacency(u, w) = .true.
                     adjacency(w, u) = .true.
                  end if
               end do
            end if
         end do
         active(v) = .false.
      end do
      info = matrix_success
   end subroutine minimum_degree_ordering

   subroutine column_degree_ordering(a, perm)
      type(csr_matrix), intent(in) :: a
      integer, allocatable, intent(out) :: perm(:)
      type(csr_matrix) :: at
      integer, allocatable :: degree(:)
      integer :: i
      call csr_transpose(a, at)
      allocate(degree(a%ncol), perm(a%ncol))
      do i = 1, a%ncol
         degree(i) = at%row_ptr(i + 1) - at%row_ptr(i)
         perm(i) = i
      end do
      call sort_by_degree(perm, degree)
   end subroutine column_degree_ordering

end module matrix_ordering
