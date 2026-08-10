program test_wide_and_rank
   use nnls, only : dp, nnls_result, nnls_fit, NNLS_SUCCESS
   implicit none
   real(dp) :: a(4,6), b(4), ar(4,3), br(4)
   type(nnls_result) :: result

   a = reshape([ &
      1.0_dp,0.0_dp,1.0_dp,2.0_dp, &
      0.0_dp,1.0_dp,1.0_dp,0.0_dp, &
      1.0_dp,1.0_dp,0.0_dp,1.0_dp, &
      2.0_dp,0.0_dp,1.0_dp,1.0_dp, &
      0.0_dp,2.0_dp,1.0_dp,1.0_dp, &
      1.0_dp,0.0_dp,2.0_dp,0.0_dp ], [4,6])
   b = [1.0_dp, 2.0_dp, 1.0_dp, 2.0_dp]
   call nnls_fit(a, b, result)
   call require(result%mode == NNLS_SUCCESS, 'wide mode')
   call require(all(result%x >= -1.0e-14_dp), 'wide nonnegative')
   call require(result%nsetp <= 4, 'at most m passive variables')

   ar(:,1) = [1.0_dp,2.0_dp,3.0_dp,4.0_dp]
   ar(:,2) = 2.0_dp * ar(:,1)
   ar(:,3) = [0.0_dp,1.0_dp,0.0_dp,1.0_dp]
   br = ar(:,1)
   call nnls_fit(ar, br, result)
   call require(result%mode == NNLS_SUCCESS, 'rank-deficient mode')
   call require(result%deviance < 1.0e-22_dp, 'rank-deficient exact fit')
   print *, 'PASS test_wide_and_rank'
contains
   subroutine require(ok, what)
      logical, intent(in) :: ok
      character(*), intent(in) :: what
      if (.not. ok) then
         print *, 'FAIL: ', what
         error stop 1
      end if
   end subroutine require
end program test_wide_and_rank
