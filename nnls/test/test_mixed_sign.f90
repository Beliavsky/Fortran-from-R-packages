program test_mixed_sign
   use nnls, only : dp, nnls_result, nnnpls_fit, NNLS_SUCCESS
   implicit none
   integer, parameter :: m = 6, n = 4
   real(dp) :: a(m,n), b(m), con(n), expected(n)
   type(nnls_result) :: result

   a = reshape([ &
      1.0_dp, 2.0_dp, 0.0_dp, 1.0_dp, 3.0_dp, -1.0_dp, &
      0.0_dp, 1.0_dp, 2.0_dp, 1.0_dp, -1.0_dp, 2.0_dp, &
      2.0_dp, -1.0_dp, 1.0_dp, 0.0_dp, 1.0_dp, 1.0_dp, &
      1.0_dp, 0.0_dp, 1.0_dp, 2.0_dp, 1.0_dp, -2.0_dp ], [m,n])
   b = [1.0_dp, 2.0_dp, -1.0_dp, 0.5_dp, 3.0_dp, -0.2_dp]
   con = [1.0_dp, -1.0_dp, 1.0_dp, -1.0_dp]
   expected = [1.1103351955307263_dp, -0.11662011173184360_dp, &
      0.08496275605214158_dp, -0.43356610800744877_dp]

   call nnnpls_fit(a, b, con, result)
   call require(result%mode == NNLS_SUCCESS, 'mode')
   call require(maxval(abs(result%x - expected)) < 2.0e-12_dp, 'coefficients')
   call require(abs(result%rnorm - 0.62018724270569914_dp) < 2.0e-12_dp, 'rnorm')
   call require(result%nsetp == 4, 'nsetp')
   call require(all(result%passive == [1,4,2,3]), 'passive order')
   call require(result%x(2) <= 0.0_dp .and. result%x(4) <= 0.0_dp, 'negative constraints')
   print *, 'PASS test_mixed_sign'
contains
   subroutine require(ok, what)
      logical, intent(in) :: ok
      character(*), intent(in) :: what
      if (.not. ok) then
         print *, 'FAIL: ', what
         error stop 1
      end if
   end subroutine require
end program test_mixed_sign
