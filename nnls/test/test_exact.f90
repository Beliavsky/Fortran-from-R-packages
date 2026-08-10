program test_exact
   use nnls, only : dp, nnls_result, nnls_fit, nnnpls_fit, NNLS_SUCCESS
   implicit none
   real(dp) :: a(5,3), xtrue(3), b(5), con(3), xmix(3)
   type(nnls_result) :: result

   a = reshape([ &
      1.0_dp, 0.0_dp, 1.0_dp, 2.0_dp, -1.0_dp, &
      0.0_dp, 1.0_dp, 1.0_dp, -1.0_dp, 2.0_dp, &
      1.0_dp, 1.0_dp, 0.0_dp, 1.0_dp, 1.0_dp ], [5,3])
   xtrue = [0.5_dp, 1.25_dp, 2.0_dp]
   b = matmul(a, xtrue)
   call nnls_fit(a, b, result)
   call require(result%mode == NNLS_SUCCESS, 'nnls mode')
   call require(maxval(abs(result%x - xtrue)) < 2.0e-12_dp, 'exact NNLS')
   call require(result%deviance < 1.0e-24_dp, 'exact deviance')

   xmix = [0.5_dp, -1.25_dp, 2.0_dp]
   con = [1.0_dp, -1.0_dp, 1.0_dp]
   b = matmul(a, xmix)
   call nnnpls_fit(a, b, con, result)
   call require(maxval(abs(result%x - xmix)) < 2.0e-12_dp, 'exact NNNPLS')
   print *, 'PASS test_exact'
contains
   subroutine require(ok, what)
      logical, intent(in) :: ok
      character(*), intent(in) :: what
      if (.not. ok) then
         print *, 'FAIL: ', what
         error stop 1
      end if
   end subroutine require
end program test_exact
