! SPDX-License-Identifier: GPL-3.0-only
program matrix_market_roundtrip
   use matrix, only : dp, csr_matrix, csr_from_dense, csr_to_dense, &
      write_matrix_market, read_matrix_market, matrix_success
   implicit none
   type(csr_matrix) :: a, b
   real(dp), allocatable :: dense(:,:)
   integer :: info, unit
   character(len=*), parameter :: filename = 'example_matrix.mtx'

   dense = reshape([2.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 3.0_dp, 0.0_dp, &
      1.0_dp, 0.0_dp, 4.0_dp], [3, 3])
   call csr_from_dense(dense, a)
   call write_matrix_market(filename, a, info, symmetric=.true.)
   if (info /= matrix_success) error stop 'write failed'
   call read_matrix_market(filename, b, info)
   if (info /= matrix_success) error stop 'read failed'
   print '(a,1x,es12.4)', 'maximum roundtrip error:', maxval(abs(csr_to_dense(b) - dense))
   open(newunit=unit, file=filename, status='old')
   close(unit, status='delete')
end program matrix_market_roundtrip
