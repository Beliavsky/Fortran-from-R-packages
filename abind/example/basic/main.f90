program abind_example
   use abind, only : dp, array_value, init_array, abind_arrays
   implicit none

   type(array_value) :: pieces(2)
   type(array_value) :: joined
   integer :: status

   call init_array(pieces(1), real([1, 2, 3, 4], dp), [2, 2], status)
   if (status /= 0) error stop 'failed to initialize first matrix'
   call init_array(pieces(2), real([5, 6, 7, 8], dp), [2, 2], status)
   if (status /= 0) error stop 'failed to initialize second matrix'
   call abind_arrays(pieces, joined, status, along=3.0_dp)
   if (status /= 0) error stop 'abind_arrays failed'

   print '(a,*(1x,i0))', 'shape:', joined%dims
   print '(a,*(1x,g0))', 'data :', joined%data
end program abind_example
