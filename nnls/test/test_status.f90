program test_status
   use nnls, only : dp, nnls_result, nnls_fit, nnnpls_fit, &
      NNLS_BAD_DIMENSIONS, NNLS_ITERATION_LIMIT
   implicit none
   real(dp), allocatable :: a0(:,:), b0(:)
   real(dp) :: a(3,2), b(3), conbad(1)
   type(nnls_result) :: result

   allocate(a0(0,2), b0(0))
   call nnls_fit(a0, b0, result)
   call require(result%mode == NNLS_BAD_DIMENSIONS, 'bad dimensions')

   a = reshape([1.0_dp,0.0_dp,1.0_dp, 0.0_dp,1.0_dp,1.0_dp], [3,2])
   b = [1.0_dp,1.0_dp,1.0_dp]
   call nnls_fit(a, b, result, max_iter=0)
   call require(result%mode == NNLS_ITERATION_LIMIT, 'iteration limit')

   conbad = [1.0_dp]
   call nnnpls_fit(a, b, conbad, result)
   call require(result%mode == NNLS_BAD_DIMENSIONS, 'bad con length')
   print *, 'PASS test_status'
contains
   subroutine require(ok, what)
      logical, intent(in) :: ok
      character(*), intent(in) :: what
      if (.not. ok) then
         print *, 'FAIL: ', what
         error stop 1
      end if
   end subroutine require
end program test_status
