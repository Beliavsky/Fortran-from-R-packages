program test_core
   use countdm
   implicit none
   integer :: failures, q, i
   integer, allocatable :: x(:), y(:), seed(:)
   real(dp) :: s, p, ref

   failures = 0
   call check_close('TP(0,3)', touchard_polynomial(0, 3.0_dp), 1.0_dp, 1.0e-14_dp)
   call check_close('TP(2,3)', touchard_polynomial(2, 3.0_dp), 12.0_dp, 1.0e-14_dp)
   call check_close('TP(4,3)', touchard_polynomial(4, 3.0_dp), 309.0_dp, 1.0e-12_dp)
   call check_close('Bell number B5', bell_number(5), 52.0_dp, 1.0e-12_dp)

   ref = 3.386239172414399e-5_dp
   call check_close('Bell-Touchard PMF reference', dbellt(2, 2.0_dp, 2.0_dp), ref, 2.0e-16_dp)

   s = 0.0_dp
   do i = 0, 100
      s = s + dbellt(i, 0.5_dp, 1.2_dp)
   end do
   call check_close('Bell-Touchard normalization', s, 1.0_dp, 5.0e-13_dp)

   do i = 0, 12
      call check_close('theta=1 Bell reduction', dbellt(i, 0.7_dp, 1.0_dp), dbell(i, 0.7_dp), 2.0e-13_dp)
   end do

   p = pbellt(3, 0.8_dp, 1.3_dp)
   q = qbellt(p, 0.8_dp, 1.3_dp)
   call check_true('Bell-Touchard quantile inversion', q <= 3)
   q = qbellt(0.75_dp, 0.8_dp, 1.3_dp)
   call check_true('Bell-Touchard quantile bracket lower', pbellt(q, 0.8_dp, 1.3_dp) >= 0.75_dp)
   if (q > 0) call check_true('Bell-Touchard quantile bracket upper', pbellt(q-1, 0.8_dp, 1.3_dp) < 0.75_dp)

   do i = 0, 8
      call check_close('ZIBT mixture PMF', dzibellt(i, 0.8_dp, 1.3_dp, 0.2_dp), &
         merge(0.2_dp + 0.8_dp * dbellt(0, 0.8_dp, 1.3_dp), &
         0.8_dp * dbellt(i, 0.8_dp, 1.3_dp), i == 0), 2.0e-13_dp)
   end do
   call check_close('ZIBT CDF identity', pzibellt(4, 0.8_dp, 1.3_dp, 0.2_dp), &
      0.2_dp + 0.8_dp * pbellt(4, 0.8_dp, 1.3_dp), 2.0e-13_dp)

   call check_close('ZIP zero mass', dzip(0, 0.3_dp, 1.2_dp), 0.3_dp + 0.7_dp * exp(-1.2_dp), 1.0e-14_dp)
   call check_close('ZOIP zero mass', dzoip(0, 0.2_dp, 0.1_dp, 1.2_dp), &
      0.2_dp + 0.7_dp * exp(-1.2_dp), 1.0e-14_dp)
   call check_close('ZOIP one mass', dzoip(1, 0.2_dp, 0.1_dp, 1.2_dp), &
      0.1_dp + 0.7_dp * 1.2_dp * exp(-1.2_dp), 1.0e-14_dp)

   call get_data_criminal(x)
   call check_true('criminal length', size(x) == 4301)
   call check_true('criminal sum', sum(x) == 334)
   call get_data_sbirth(x)
   call check_true('sbirth length', size(x) == 402)
   call check_true('sbirth sum', sum(x) == 175)

   call random_seed(size=i)
   allocate(seed(i))
   seed = 1729
   call random_seed(put=seed)
   allocate(y(200))
   call rbellt(size(y), 0.7_dp, 1.1_dp, y)
   call check_true('Bell-Touchard RNG support', all(y >= 0))
   call rzibellt(size(y), 0.7_dp, 1.1_dp, 0.3_dp, y)
   call check_true('ZIBT RNG support', all(y >= 0))

   if (failures == 0) then
      print '(a)', 'test_core: PASS'
   else
      print '(a,i0)', 'test_core: FAIL ', failures
      error stop 1
   end if

contains

   subroutine check_close(name, got, expected, tol)
      character(*), intent(in) :: name
      real(dp), intent(in) :: got, expected, tol
      if (abs(got - expected) > tol * max(1.0_dp, abs(expected))) then
         failures = failures + 1
         print '(a,2(1x,es24.16))', 'failure: '//trim(name), got, expected
      end if
   end subroutine check_close

   subroutine check_true(name, condition)
      character(*), intent(in) :: name
      logical, intent(in) :: condition
      if (.not. condition) then
         failures = failures + 1
         print '(a)', 'failure: '//trim(name)
      end if
   end subroutine check_true

end program test_core
