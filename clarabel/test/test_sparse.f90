program test_sparse
   use, intrinsic :: iso_c_binding, only : c_size_t
   use clarabel
   implicit none
   real(dp) :: x(3,2), s(3,3)
   real(dp), allocatable :: y(:,:)
   type(csc_matrix) :: a, p
   logical :: ok
   character(len=:), allocatable :: message
   x = reshape([1.0_dp,0.0_dp,2.0_dp, 0.0_dp,3.0_dp,0.0_dp], shape(x))
   a = csc_from_dense(x)
   call a%validate(ok,message); if(.not.ok) error stop message
   if(a%nnz()/=3) error stop "bad nnz"
   if(any(a%colptr /= [0_c_size_t,2_c_size_t,3_c_size_t])) error stop "bad colptr"
   y=a%to_dense(); if(maxval(abs(y-x))>1e-14_dp) error stop "roundtrip"
   s=reshape([2.0_dp,1.0_dp,0.0_dp,1.0_dp,3.0_dp,4.0_dp,0.0_dp,4.0_dp,5.0_dp],shape(s))
   p=csc_from_symmetric_upper(s); y=p%to_dense(.true.)
   if(maxval(abs(y-s))>1e-14_dp) error stop "symmetric roundtrip"
   print *, "test_sparse: PASS"
end program test_sparse
