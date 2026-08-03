program test_sparse
   use osqp
   implicit none
   real(dp) :: a(3,3), expected(3,3)
   integer(osqp_int) :: row(5), col(5), status
   real(dp) :: value(5)
   type(osqp_sparse_matrix) :: s, t

   a = reshape([1.0_dp,0.0_dp,4.0_dp, 2.0_dp,3.0_dp,0.0_dp, 0.0_dp,5.0_dp,6.0_dp],[3,3])
   s = osqp_csc_from_dense(a)
   call check(s%valid(), "dense CSC valid")
   call check(maxval(abs(s%to_dense()-a)) < 1.0e-14_dp, "dense round trip")

   row = [1,1,2,3,3]
   col = [1,1,2,1,3]
   value = [0.25_dp,0.75_dp,3.0_dp,4.0_dp,6.0_dp]
   t = osqp_csc_from_triplet(3,3,row,col,value,status=status)
   expected = 0.0_dp
   expected(1,1)=1.0_dp; expected(2,2)=3.0_dp
   expected(3,1)=4.0_dp; expected(3,3)=6.0_dp
   call check(status == 0 .and. t%valid(), "triplet construction")
   call check(maxval(abs(t%to_dense()-expected)) < 1.0e-14_dp, "duplicate aggregation")

   t = osqp_csc_from_dense(a, upper_only=.true.)
   expected = a
   expected(2,1)=0.0_dp; expected(3,1)=0.0_dp; expected(3,2)=0.0_dp
   call check(maxval(abs(t%to_dense()-expected)) < 1.0e-14_dp, "upper triangular extraction")
   print *, "PASS test_sparse"
contains
   subroutine check(ok, message)
      logical, intent(in) :: ok
      character(len=*), intent(in) :: message
      if (.not. ok) error stop message
   end subroutine check
end program test_sparse
