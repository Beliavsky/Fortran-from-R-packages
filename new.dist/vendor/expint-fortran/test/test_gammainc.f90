program test_gammainc
   use expint_mod
   use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
   implicit none

   real(dp), parameter :: tol = 2.0e-11_dp
   real(dp), parameter :: a(8) = [ &
      -3.2_dp, -2.0_dp, -1.2_dp, -0.499_dp, -0.25_dp, 0.2_dp, 1.2_dp, 5.5_dp ]
   real(dp), parameter :: ref02(8) = [ &
      40.616540196120045952_dp, 8.798632802871763995_dp, 3.533666787642764806_dp, &
      1.791447907955587255_dp, 1.462097540200405147_dp, 1.081443946818191376_dp, &
      0.809688189222084353_dp, 52.34275580804438042_dp ]
   real(dp), parameter :: ref25(8) = [ &
      0.000712096836339791193_dp, 0.002607259100267012327_dp, &
      0.006325292202403112735_dp, 0.01399234045691350902_dp, &
      0.01863567728445995988_dp, 0.03149504331501303358_dp, &
      0.1048933061734996029_dp, 48.73984697226371190_dp ]
   integer :: i

   do i = 1, size(a)
      call check_close(gamma_inc(a(i), 0.2_dp), ref02(i), tol, 'Gamma(a,0.2)')
      call check_close(gamma_inc(a(i), 2.5_dp), ref25(i), tol, 'Gamma(a,2.5)')
   end do

   call check_close(gammainc(0.0_dp, 2.5_dp), expint_e1(2.5_dp), tol, 'Gamma(0,x)=E1')
   call check_close(gammainc(1.0_dp, 2.5_dp), exp(-2.5_dp), tol, 'Gamma(1,x)')
   call check_close(gammainc(2.0_dp, 2.5_dp), 3.5_dp * exp(-2.5_dp), tol, 'Gamma(2,x)')

   call check_close(gammainc(-0.499_dp, 1.0e-5_dp), 622.92999782140934_dp, 5.0e-10_dp, &
      'small-x recursion near -0.5')

   if (.not. ieee_is_nan(gammainc(1.2_dp, -1.0_dp))) error stop 'negative x should be NaN'

   print '(a)', 'test_gammainc: PASS'

contains

   subroutine check_close(value, reference, reltol, label)
      real(dp), intent(in) :: value, reference, reltol
      character(*), intent(in) :: label
      real(dp) :: scale

      scale = max(1.0_dp, abs(reference))
      if (abs(value - reference) > reltol * scale) then
         print '(a)', trim(label)
         print '(a,es24.16)', '  value:     ', value
         print '(a,es24.16)', '  reference: ', reference
         error stop 'numeric comparison failed'
      end if
   end subroutine check_close

end program test_gammainc
