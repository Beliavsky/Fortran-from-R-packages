program test_sequences
use qrng, only: dp, korobov, ghalton, sobol
implicit none
real(dp), allocatable :: u(:,:)
integer :: shifts(3,32)
real(dp), parameter :: tol = 2.0e-14_dp

u = korobov(7, 3, 3)
call assert_close(u(:,1), [0.0_dp,1.0_dp/7,2.0_dp/7,3.0_dp/7,4.0_dp/7,5.0_dp/7,6.0_dp/7], tol, "korobov col1")
call assert_close(u(:,2), [0.0_dp,3.0_dp/7,6.0_dp/7,2.0_dp/7,5.0_dp/7,1.0_dp/7,4.0_dp/7], tol, "korobov col2")
call assert_close(u(:,3), [0.0_dp,2.0_dp/7,4.0_dp/7,6.0_dp/7,1.0_dp/7,3.0_dp/7,5.0_dp/7], tol, "korobov col3")

shifts = 0
u = ghalton(6, 3, method="halton", shift_coeff=shifts)
call assert_close(u(:,1), [0.0_dp,0.5_dp,0.25_dp,0.75_dp,0.125_dp,0.625_dp], tol, "halton dim1")
call assert_close(u(:,2), [0.0_dp,1.0_dp/3,2.0_dp/3,1.0_dp/9,4.0_dp/9,7.0_dp/9], tol, "halton dim2")
call assert_close(u(:,3), [0.0_dp,0.2_dp,0.4_dp,0.6_dp,0.8_dp,0.04_dp], tol, "halton dim3")

u = ghalton(6, 3, method="generalized", shift_coeff=shifts)
call assert_close(u(:,1), [0.0_dp,0.5_dp,0.25_dp,0.75_dp,0.125_dp,0.625_dp], tol, "ghalton dim1")
call assert_close(u(:,2), [0.0_dp,1.0_dp/3,2.0_dp/3,1.0_dp/9,4.0_dp/9,7.0_dp/9], tol, "ghalton dim2")
call assert_close(u(:,3), [0.0_dp,0.6_dp,0.2_dp,0.8_dp,0.4_dp,0.12_dp], tol, "ghalton dim3")

u = sobol(8,4)
call assert_close(u(:,1), [0.0_dp,0.5_dp,0.75_dp,0.25_dp,0.375_dp,0.875_dp,0.625_dp,0.125_dp], tol, "sobol dim1")
call assert_close(u(:,2), [0.0_dp,0.5_dp,0.25_dp,0.75_dp,0.375_dp,0.875_dp,0.125_dp,0.625_dp], tol, "sobol dim2")
call assert_close(u(:,3), [0.0_dp,0.5_dp,0.25_dp,0.75_dp,0.625_dp,0.125_dp,0.875_dp,0.375_dp], tol, "sobol dim3")
call assert_close(u(:,4), [0.0_dp,0.5_dp,0.75_dp,0.25_dp,0.875_dp,0.375_dp,0.125_dp,0.625_dp], tol, "sobol dim4")

u = sobol(5,3,skip=5)
call assert_close(u(:,1), [0.875_dp,0.625_dp,0.125_dp,0.1875_dp,0.6875_dp], tol, "sobol skip dim1")
call assert_close(u(:,2), [0.875_dp,0.125_dp,0.625_dp,0.3125_dp,0.8125_dp], tol, "sobol skip dim2")
call assert_close(u(:,3), [0.125_dp,0.875_dp,0.375_dp,0.9375_dp,0.4375_dp], tol, "sobol skip dim3")

u = sobol(8,16510)
call assert_close(u(:,16510), [0.0_dp,0.5_dp,0.25_dp,0.75_dp,0.625_dp,0.125_dp,0.875_dp,0.375_dp], &
   tol, "sobol max dimension")

print '(a)', 'test_sequences: PASS'

contains
subroutine assert_close(x,y,tol,label)
real(dp), intent(in) :: x(:), y(:), tol
character(len=*), intent(in) :: label
if (size(x) /= size(y) .or. any(abs(x-y) > tol)) then
   print '(a)', 'FAILED: '//label
   if (size(x)==size(y)) then
      print *, x
      print *, y
   end if
   error stop 1
end if
end subroutine assert_close
end program test_sequences
