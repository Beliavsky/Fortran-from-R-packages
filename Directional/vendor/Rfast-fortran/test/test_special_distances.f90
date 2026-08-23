program test_special_distances
   use rfast
   implicit none
   real(dp) :: x(3,2), y(3,2), d(3,3), c

   call assert_close(normal_cdf(0.0_dp),0.5_dp,1e-15_dp,'normal cdf')
   call assert_close(normal_quantile(0.975_dp),1.95996398454005_dp,2e-9_dp,'normal quantile')
   call assert_close(digamma_r(1.0_dp),-0.577215664901533_dp,2e-12_dp,'digamma')
   call assert_close(trigamma_r(1.0_dp),1.644934066848226_dp,2e-12_dp,'trigamma')
   call assert_close(chisq_cdf(3.84145882069412_dp,1.0_dp),0.95_dp,2e-10_dp,'chi square')

   x = reshape([0.0_dp,1.0_dp,2.0_dp,0.0_dp,0.0_dp,0.0_dp],[3,2])
   d = dist_matrix(x)
   call assert_close(d(1,2),1.0_dp,1e-12_dp,'distance 1')
   call assert_close(d(1,3),2.0_dp,1e-12_dp,'distance 2')
   y = 2.0_dp*x
   c = distance_correlation(x,y)
   call assert_close(c,1.0_dp,1e-10_dp,'distance correlation')
   call assert_close(energy_distance(x,x),0.0_dp,1e-12_dp,'energy identical')

   print *, 'test_special_distances: PASS'
contains
   subroutine assert_close(got,want,tol,msg)
      real(dp),intent(in)::got,want,tol
      character(*),intent(in)::msg
      if(abs(got-want)>tol)then
         print *, 'FAIL ',trim(msg),got,want
         error stop 1
      end if
   end subroutine
end program test_special_distances
