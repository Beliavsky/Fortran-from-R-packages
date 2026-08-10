program test_reference
   use nnls, only : dp, nnls_result, nnls_fit, NNLS_SUCCESS
   implicit none
   integer, parameter :: m = 6, n = 4
   real(dp) :: a(m,n), b(m), expected(n)
   type(nnls_result) :: result

   a = reshape([ &
      1.0_dp, 2.0_dp, 0.0_dp, 1.0_dp, 3.0_dp, -1.0_dp, &
      0.0_dp, 1.0_dp, 2.0_dp, 1.0_dp, -1.0_dp, 2.0_dp, &
      2.0_dp, -1.0_dp, 1.0_dp, 0.0_dp, 1.0_dp, 1.0_dp, &
      1.0_dp, 0.0_dp, 1.0_dp, 2.0_dp, 1.0_dp, -2.0_dp ], [m,n])
   b = [1.0_dp, 2.0_dp, -1.0_dp, 0.5_dp, 3.0_dp, -0.2_dp]
   expected = [0.91875_dp, 0.0_dp, 0.0_dp, 0.0_dp]

   call nnls_fit(a, b, result)
   call require(result%mode == NNLS_SUCCESS, 'mode')
   call require(maxval(abs(result%x - expected)) < 1.0e-12_dp, 'coefficients')
   call require(abs(result%rnorm - 1.3358050007392546_dp) < 1.0e-12_dp, 'rnorm')
   call require(result%nsetp == 1, 'nsetp')
   call require(size(result%passive) == 1 .and. result%passive(1) == 1, 'passive set')
   print *, 'PASS test_reference'
contains
   subroutine require(ok, what)
      logical, intent(in) :: ok
      character(*), intent(in) :: what
      if (.not. ok) then
         print *, 'FAIL: ', what
         error stop 1
      end if
   end subroutine require
end program test_reference
