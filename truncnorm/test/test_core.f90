program test_core
use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_positive_inf, ieee_negative_inf
use truncnorm, only: dp, dtruncnorm, ptruncnorm, etruncnorm, vtruncnorm
implicit none
real(dp) :: pinf, ninf
pinf = ieee_value(0.0_dp, ieee_positive_inf)
ninf = ieee_value(0.0_dp, ieee_negative_inf)

call check_close(dtruncnorm(0.2_dp, -1.0_dp, 2.0_dp, 0.3_dp, 1.4_dp), &
   0.39969448016080844_dp, 2.0e-12_dp, 'density general')
call check_close(ptruncnorm(0.2_dp, -1.0_dp, 2.0_dp, 0.3_dp, 1.4_dp), &
   0.41479724658819545_dp, 2.0e-12_dp, 'cdf general')
call check_close(etruncnorm(-1.0_dp, 2.0_dp, 0.3_dp, 1.4_dp), &
   0.4345779731683568_dp, 2.0e-12_dp, 'mean general')
call check_close(vtruncnorm(-1.0_dp, 2.0_dp, 0.3_dp, 1.4_dp), &
   0.6396569422409290_dp, 2.0e-11_dp, 'variance general')

call check_close(dtruncnorm(-4.5_dp, -5.0_dp, -4.0_dp, 0.0_dp, 1.0_dp), &
   0.5092862762858839_dp, 5.0e-10_dp, 'left-tail density')
call check_close(etruncnorm(-5.0_dp, -4.0_dp, 0.0_dp, 1.0_dp), &
   -4.216830780601033_dp, 5.0e-10_dp, 'left-tail mean')
call check_close(vtruncnorm(-5.0_dp, -4.0_dp, 0.0_dp, 1.0_dp), &
   0.03829028732908668_dp, 2.0e-9_dp, 'left-tail variance')

call check_close(etruncnorm(0.0_dp, 1.0_dp, 5.0_dp, 0.1_dp), 0.5_dp, 0.0_dp, 'upstream extreme mean fallback')
call check_close(vtruncnorm(0.0_dp, 1.0_dp, 5.0_dp, 0.1_dp), 1.0_dp/12.0_dp, 0.0_dp, &
   'upstream extreme variance fallback')
call check_close(dtruncnorm(0.25_dp, 0.0_dp, 1.0_dp, 5.0_dp, 0.1_dp), 1.0_dp, 0.0_dp, &
   'upstream extreme density fallback')

call check_close(etruncnorm(ninf, pinf, 1.25_dp, 2.5_dp), 1.25_dp, 0.0_dp, 'untruncated mean')
call check_close(vtruncnorm(ninf, pinf, 1.25_dp, 2.5_dp), 6.25_dp, 0.0_dp, 'untruncated variance')
call check_close(ptruncnorm(-2.0_dp, -1.0_dp, 1.0_dp, 0.0_dp, 1.0_dp), 0.0_dp, 0.0_dp, 'cdf below')
call check_close(ptruncnorm(2.0_dp, -1.0_dp, 1.0_dp, 0.0_dp, 1.0_dp), 1.0_dp, 0.0_dp, 'cdf above')
call check_close(dtruncnorm(2.0_dp, -1.0_dp, 1.0_dp, 0.0_dp, 1.0_dp), 0.0_dp, 0.0_dp, 'density outside')

print *, 'test_core: PASS'
contains
subroutine check_close(x, ref, tol, label)
real(dp), intent(in) :: x, ref, tol
character(*), intent(in) :: label
if (abs(x-ref) > tol) then
   write(*,'(a,1x,a,2(1x,es24.16))') 'FAIL', trim(label), x, ref
   error stop 1
end if
end subroutine check_close
end program test_core
