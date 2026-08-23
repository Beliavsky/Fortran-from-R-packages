program test_arrays_linalg
   use rfast
   implicit none
   real(dp) :: x(5), a(4,3), s(3,3), l(3,3), inv(3,3), eye(3,3)
   real(dp), allocatable :: r(:)
   real(dp) :: cm(3)
   integer, allocatable :: ord(:)
   integer :: info, i

   x = [3.0_dp, 1.0_dp, 4.0_dp, 1.0_dp, 5.0_dp]
   call assert_close(mean_r(x), 2.8_dp, 1e-12_dp, 'mean')
   call assert_close(median_r(x), 3.0_dp, 1e-12_dp, 'median')
   call assert_close(variance_r(x), 3.2_dp, 1e-12_dp, 'variance')
   ord = order_real(x)
   call assert_true(all(ord == [2,4,1,3,5]), 'order')
   r = rank_average(x)
   call assert_close(r(2), 1.5_dp, 1e-12_dp, 'rank tie')
   call assert_close(r(5), 5.0_dp, 1e-12_dp, 'rank max')

   a = reshape([1.0_dp,2.0_dp,3.0_dp,4.0_dp, &
                2.0_dp,4.0_dp,6.0_dp,8.0_dp, &
                1.0_dp,1.0_dp,2.0_dp,2.0_dp], [4,3])
   cm = colmeans(a)
   call assert_close(cm(1), 2.5_dp, 1e-12_dp, 'colmeans')

   s = reshape([4.0_dp,1.0_dp,1.0_dp, &
                1.0_dp,3.0_dp,0.5_dp, &
                1.0_dp,0.5_dp,2.0_dp], [3,3])
   call cholesky_lower(s,l,info)
   call assert_true(info == 0, 'cholesky info')
   call assert_true(maxval(abs(matmul(l,transpose(l))-s)) < 1e-11_dp, 'cholesky')
   call inverse_matrix(s,inv,info)
   call assert_true(info == 0, 'inverse info')
   eye = matmul(s,inv)
   do i=1,3
      eye(i,i) = eye(i,i) - 1.0_dp
   end do
   call assert_true(maxval(abs(eye)) < 1e-11_dp, 'inverse')
   call assert_true(matrix_rank(a) == 2, 'rank matrix')

   print *, 'test_arrays_linalg: PASS'
contains
   subroutine assert_close(got,want,tol,msg)
      real(dp),intent(in)::got,want,tol
      character(*),intent(in)::msg
      if(abs(got-want)>tol)then
         print *, 'FAIL ',trim(msg),got,want
         error stop 1
      end if
   end subroutine
   subroutine assert_true(ok,msg)
      logical,intent(in)::ok
      character(*),intent(in)::msg
      if(.not.ok)then
         print *, 'FAIL ',trim(msg)
         error stop 1
      end if
   end subroutine
end program test_arrays_linalg
