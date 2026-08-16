program test_distribution
   use ccd, only : dp, i8, dcc, pcc, qcc
   implicit none
   integer :: fail
   integer(i8) :: q(5)
   real(dp), parameter :: tol = 5.0e-13_dp
   fail = 0

   call check_close(dcc(0.0_dp, 0.0_dp, 1.5_dp), 0.21217234361392046_dp, tol, fail)
   call check_close(dcc(2.0_dp, 0.0_dp, 1.5_dp), 0.07638204370101137_dp, tol, fail)
   call check_close(pcc(-3_i8, 0.0_dp, 1.5_dp), 0.17064323891316036_dp, tol, fail)
   call check_close(pcc(0_i8, 0.0_dp, 1.5_dp), 0.6060861718069602_dp, tol, fail)
   call check_close(pcc(2_i8, 0.0_dp, 1.5_dp), 0.8293567610868396_dp, tol, fail)

   q = qcc([0.10_dp, 0.25_dp, 0.50_dp, 0.75_dp, 0.90_dp], 0.0_dp, 1.5_dp)
   if (any(q /= [-5_i8, -1_i8, 0_i8, 1_i8, 5_i8])) fail = fail + 1

   if (fail == 0) then
      print *, 'test_distribution: PASS'
   else
      print *, 'test_distribution: FAIL', fail
      error stop 1
   end if
contains
   subroutine check_close(x, ref, eps, fail)
      real(dp), intent(in) :: x, ref, eps
      integer, intent(inout) :: fail
      if (abs(x-ref) > eps) then
         print *, 'mismatch:', x, ref
         fail = fail + 1
      end if
   end subroutine check_close
end program test_distribution
