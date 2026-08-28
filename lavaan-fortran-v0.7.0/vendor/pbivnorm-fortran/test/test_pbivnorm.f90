! SPDX-License-Identifier: GPL-2.0-or-later
program test_pbivnorm
   use pbivnorm_mod, only : dp, pbivnorm, pbivnorm_recycle
   use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
   implicit none
   real(dp), parameter :: tol = 2.0e-12_dp
   real(dp) :: p, x(4), y(4), r(4), out(4)
   integer :: status

   call check_close(pbivnorm(0.0_dp, 0.0_dp, 0.0_dp), 0.25_dp, tol, 'independent origin')
   call check_close(pbivnorm(0.0_dp, 0.0_dp, 0.5_dp), 1.0_dp/3.0_dp, tol, 'rho=.5 origin')
   call check_close(pbivnorm(0.0_dp, 0.0_dp, -0.5_dp), 1.0_dp/6.0_dp, tol, 'rho=-.5 origin')
   call check_close(pbivnorm(1.0_dp, -0.5_dp, 0.0_dp), &
                    0.8413447460685429_dp*0.3085375387259869_dp, tol, 'independence')

   call check_close(pbivnorm(0.4_dp, -0.7_dp, 1.0_dp), &
                    0.24196365222307303_dp, 5.0e-15_dp, 'rho=1')
   call check_close(pbivnorm(0.4_dp, -0.7_dp, -1.0_dp), 0.0_dp, 5.0e-15_dp, 'rho=-1 empty')
   call check_close(pbivnorm(1.2_dp, 0.4_dp, -1.0_dp), &
                    0.5403520713886159_dp, 2.0e-15_dp, 'rho=-1 interval')

   p = pbivnorm(huge(1.0_dp), 0.3_dp, 0.8_dp)
   call check_close(p, 0.6179114221889526_dp, 2.0e-15_dp, 'large finite x')

   if (.not. ieee_is_nan(pbivnorm(0.0_dp, 0.0_dp, 1.1_dp))) error stop 'invalid rho should return NaN'

   x = [-1.0_dp, 0.0_dp, 0.5_dp, 1.0_dp]
   y = [ 0.5_dp, 0.0_dp, 1.0_dp, 2.0_dp]
   r = [ 0.0_dp, 0.5_dp, -0.2_dp, 0.9_dp]
   out = pbivnorm(x, y, r)
   call check_close(out(1), pbivnorm(x(1), y(1), r(1)), tol, 'elemental 1')
   call check_close(out(4), pbivnorm(x(4), y(4), r(4)), tol, 'elemental 4')

   call pbivnorm_recycle(x(1:2), y, [0.25_dp], out, status)
   if (status /= 0) error stop 'recycle status'
   call check_close(out(3), pbivnorm(x(1), y(3), 0.25_dp), tol, 'recycle')

   print '(a)', 'test_pbivnorm: PASS'

contains

   subroutine check_close(actual, expected, atol, label)
      real(dp), intent(in) :: actual, expected, atol
      character(*), intent(in) :: label
      if (abs(actual - expected) > atol) then
         print '(a,2(1x,es24.15))', trim(label)//' actual expected:', actual, expected
         error stop 1
      end if
   end subroutine check_close

end program test_pbivnorm
