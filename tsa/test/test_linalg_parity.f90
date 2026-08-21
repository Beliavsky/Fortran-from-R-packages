program test_linalg_parity
  use tsa, only : dp
  use tseries_linalg, only : right_singular_vectors, invert_matrix_lu
  implicit none
  real(dp) :: a(7,3), v(3,3), s(3), gramv(3,3), b(7,3)
  real(dp) :: m(3,3), inv(3,3), ident(3,3)
  integer :: i, status

  do i = 1, 7
    a(i,1) = 1.0_dp + 0.1_dp*real(i,dp)
    a(i,2) = 0.7_dp*a(i,1) + 0.03_dp*sin(real(i,dp))
    a(i,3) = -0.2_dp*a(i,1) + 0.4_dp*cos(0.7_dp*real(i,dp))
  end do
  call right_singular_vectors(a,v,s,status)
  call check(status == 0,'Jacobi SVD status')
  gramv = matmul(transpose(v),v)
  call matrix_close(gramv,eye3(),2.0e-10_dp,'V orthogonality')
  b = matmul(a,v)
  call close(dot_product(b(:,1),b(:,2)),0.0_dp,2.0e-10_dp,'SVD col12')
  call close(dot_product(b(:,1),b(:,3)),0.0_dp,2.0e-10_dp,'SVD col13')
  call close(dot_product(b(:,2),b(:,3)),0.0_dp,2.0e-10_dp,'SVD col23')
  call check(s(1) >= s(2) .and. s(2) >= s(3),'singular-value order')

  m = reshape([4.0_dp,1.0_dp,2.0_dp, &
               1.0_dp,3.0_dp,0.5_dp, &
               2.0_dp,0.5_dp,5.0_dp],[3,3])
  call invert_matrix_lu(m,inv,status)
  call check(status == 0,'pivoted LU inverse status')
  ident = matmul(m,inv)
  call matrix_close(ident,eye3(),2.0e-12_dp,'pivoted LU inverse')

  print '(a)', 'test_linalg_parity: PASS'
contains
  function eye3() result(e)
    real(dp) :: e(3,3)
    e = 0.0_dp
    e(1,1)=1.0_dp; e(2,2)=1.0_dp; e(3,3)=1.0_dp
  end function eye3

  subroutine matrix_close(x,y,tol,msg)
    real(dp), intent(in) :: x(:,:), y(:,:), tol
    character(len=*), intent(in) :: msg
    if (maxval(abs(x-y)) > tol) then
      print '(a,es24.14)', trim(msg)//' FAIL maxerr: ',maxval(abs(x-y))
      error stop 1
    end if
  end subroutine matrix_close

  subroutine close(x,y,tol,msg)
    real(dp), intent(in) :: x,y,tol
    character(len=*), intent(in) :: msg
    if (abs(x-y) > tol*max(1.0_dp,abs(y))) then
      print '(a,2es24.14)', trim(msg)//' FAIL: ',x,y
      error stop 1
    end if
  end subroutine close

  subroutine check(ok,msg)
    logical, intent(in) :: ok
    character(len=*), intent(in) :: msg
    if (.not. ok) then
      print '(a)', trim(msg)//' FAIL'
      error stop 1
    end if
  end subroutine check
end program test_linalg_parity
