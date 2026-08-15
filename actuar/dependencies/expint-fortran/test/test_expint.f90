program test_expint
   use expint_mod
   use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
   implicit none

   real(dp), parameter :: tol = 5.0e-12_dp
   real(dp) :: x(5), expected(5), got(5), scaled(5)
   real(dp) :: xv(6)
   integer :: ord(6)

   x = [0.2_dp, 1.0_dp, 2.5_dp, 5.0_dp, 10.0_dp]
   expected = [ &
      1.2226505441838930429_dp, &
      0.21938393439552027368_dp, &
      0.024914917870269735496_dp, &
      0.0011482955912753257973_dp, &
      0.0000041569689296853242774_dp ]
   got = expint_e1(x)
   call check_vec(got, expected, tol, 'E1 values')

   scaled = expint_e1(x, .true.)
   call check_vec(scaled, expected * exp(x), tol, 'scaled E1 identity')

   call check_close(expint_ei(0.5_dp), 0.45421990486317357992_dp, tol, 'Ei(0.5)')
   call check_close(expint_e1(-2.0_dp), -4.9542343560018901634_dp, tol, 'E1(-2)')

   call check_close(expint_en(1.275_dp, 3), 0.076030309585570385476_dp, tol, 'E3(1.275)')
   call check_close(expint_en(1.275_dp, 10), 0.026846986405719063706_dp, tol, 'E10(1.275)')
   call check_close(expint_en(10.0_dp, 10), 2.3253026570282108178e-6_dp, tol, 'E10(10)')
   call check_close(expint_en(0.0_dp, 4), 1.0_dp / 3.0_dp, tol, 'En(0)')
   call check_close(expint_en(2.0_dp, 0), exp(-2.0_dp) / 2.0_dp, tol, 'E0 identity')

   xv = [0.2_dp, 0.4_dp, 0.7_dp, 1.0_dp, 1.5_dp, 2.0_dp]
   ord = [1, 2, 3, 4, 5, 6]
   got = 0.0_dp
   call check_close(expint(xv(1), ord(1)), expint_e1(xv(1)), tol, 'generic E1')
   call check_close(expint(xv(2), ord(2)), expint_e2(xv(2)), tol, 'generic E2')
   call check_close(expint(xv(3), ord(3)), expint_en(xv(3), 3), tol, 'generic E3')
   call check_vec(expint(xv, ord), [ &
      expint_en(xv(1), 1), expint_en(xv(2), 2), expint_en(xv(3), 3), &
      expint_en(xv(4), 4), expint_en(xv(5), 5), expint_en(xv(6), 6)], tol, &
      'elemental vectorization')

   if (.not. ieee_is_nan(expint_e1(0.0_dp))) error stop 'E1(0) should be NaN'
   if (.not. ieee_is_nan(expint_en(-1.0_dp, 3))) error stop 'En negative x, n>2 should be NaN'
   if (.not. ieee_is_nan(expint_en(1.0_dp, -1))) error stop 'negative order should be NaN'

   call check_vec(expint_recycle([0.2_dp, 0.4_dp, 0.6_dp, 0.8_dp], [1, 2]), [ &
      expint_e1(0.2_dp), expint_e2(0.4_dp), expint_e1(0.6_dp), expint_e2(0.8_dp)], &
      tol, 'R-style argument recycling')

   print '(a)', 'test_expint: PASS'

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

   subroutine check_vec(value, reference, reltol, label)
      real(dp), intent(in) :: value(:), reference(:), reltol
      character(*), intent(in) :: label
      integer :: i

      if (size(value) /= size(reference)) error stop 'vector size mismatch'
      do i = 1, size(value)
         call check_close(value(i), reference(i), reltol, label)
      end do
   end subroutine check_vec

end program test_expint
