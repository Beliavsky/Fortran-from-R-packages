program test_map_filter
   use matrixextra
   implicit none
   real(dp) :: x(3,3), want(3,3)
   type(csr_matrix) :: a,b
   x=reshape([1.0_dp,-2.0_dp,0.0_dp,0.0_dp,3.0_dp,-4.0_dp,5.0_dp,0.0_dp,6.0_dp],[3,3])
   call csr_from_dense(x,a)
   call csr_filter(a,keep_positive,b)
   want=0.0_dp; where(x>0.0_dp) want=x
   if (maxval(abs(csr_to_dense(b)-want))>1e-12_dp) error stop 1
   call csr_map(a,square_value,b)
   want=x*x
   if (maxval(abs(csr_to_dense(b)-want))>1e-12_dp) error stop 2
   print *, 'PASS test_map_filter'
contains
   logical function keep_positive(i,j,v)
      integer, intent(in) :: i,j
      real(dp), intent(in) :: v
      keep_positive=v>0.0_dp .and. i>=1 .and. j>=1
   end function keep_positive
   real(dp) function square_value(i,j,v)
      integer, intent(in) :: i,j
      real(dp), intent(in) :: v
      square_value=v*v+0.0_dp*real(i+j,dp)
   end function square_value
end program
