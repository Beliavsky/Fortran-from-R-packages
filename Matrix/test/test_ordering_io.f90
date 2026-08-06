! SPDX-License-Identifier: GPL-3.0-only
program test_ordering_io
   use matrix, only : dp, csr_matrix, csr_from_dense, csr_to_dense, reverse_cuthill_mckee, &
      minimum_degree_ordering, column_degree_ordering, write_matrix_market, read_matrix_market, &
      write_dense_text, read_dense_text, matrix_success
   implicit none
   type(csr_matrix) :: a, b
   real(dp), allocatable :: dense(:,:), dense2(:,:)
   integer, allocatable :: p(:), q(:), cperm(:)
   integer :: info
   character(len=*), parameter :: mmfile = 'test_matrix_market.mtx'
   character(len=*), parameter :: txtfile = 'test_dense_matrix.txt'

   dense = reshape([1.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
                    1.0_dp, 1.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, &
                    0.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 0.0_dp, &
                    0.0_dp, 0.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, &
                    0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, 1.0_dp], [5, 5])
   call csr_from_dense(dense, a)
   call reverse_cuthill_mckee(a, p, info)
   call check(info == matrix_success .and. valid_perm(p), 'rcm')
   call minimum_degree_ordering(a, q, info)
   call check(info == matrix_success .and. valid_perm(q), 'minimum degree')
   call column_degree_ordering(a, cperm)
   call check(valid_perm(cperm), 'column ordering')

   call write_matrix_market(mmfile, a, info, symmetric=.true.)
   call check(info == matrix_success, 'write matrix market')
   call read_matrix_market(mmfile, b, info)
   call check(info == matrix_success, 'read matrix market')
   call check_close(csr_to_dense(b), dense, 1.0e-14_dp, 'matrix market roundtrip')

   call write_dense_text(txtfile, dense, info)
   call read_dense_text(txtfile, dense2, info)
   call check_close(dense2, dense, 1.0e-14_dp, 'dense text roundtrip')
   call delete_file(mmfile)
   call delete_file(txtfile)

   print '(a)', 'test_ordering_io: PASS'
contains
   logical function valid_perm(x) result(ok)
      integer, intent(in) :: x(:)
      logical, allocatable :: seen(:)
      integer :: i, n
      n = size(x)
      if (any(x < 1) .or. any(x > n)) then
         ok = .false.
         return
      end if
      allocate(seen(n), source=.false.)
      do i = 1, n
         if (seen(x(i))) then
            ok = .false.
            return
         end if
         seen(x(i)) = .true.
      end do
      ok = .true.
   end function valid_perm


   subroutine delete_file(filename)
      character(len=*), intent(in) :: filename
      integer :: unit, ios
      open(newunit=unit, file=filename, status='old', iostat=ios)
      if (ios == 0) close(unit, status='delete')
   end subroutine delete_file

   subroutine check(condition, name)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: name
      if (.not. condition) then
         print '(a)', 'FAIL: ' // name
         error stop 1
      end if
   end subroutine check

   subroutine check_close(actual, expected, tol, name)
      real(dp), intent(in) :: actual(:,:), expected(:,:)
      real(dp), intent(in) :: tol
      character(len=*), intent(in) :: name
      call check(all(shape(actual) == shape(expected)) .and. maxval(abs(actual - expected)) <= tol, name)
   end subroutine check_close
end program test_ordering_io
