! SPDX-License-Identifier: GPL-2.0-only
program test_cross_validation
   use fincovregularization
   implicit none
   real(dp) :: data(30,5)
   type(cv_result) :: band_result, band_repeat, taper_result, hard_result, soft_result, operator_result
   integer :: i, j

   do i = 1, size(data,1)
      do j = 1, size(data,2)
         data(i,j) = 0.01_dp*sin(0.17_dp*real(i*j,dp)) + &
            0.004_dp*cos(0.11_dp*real(i*(j+1),dp)) + 0.001_dp*real(j,dp)
      end do
   end do

   band_result = banding_cv(data, n_cv=4, norm='F', seed=1234)
   call assert_true(band_result%status == fincov_ok, 'banding CV status')
   call assert_true(size(band_result%cv_error) == 5, 'banding CV grid size')
   call assert_true(band_result%parameter_opt >= 0.0_dp .and. band_result%parameter_opt <= 4.0_dp, 'banding CV optimum')

   band_repeat = banding_cv(data, n_cv=4, norm='f', seed=1234)
   call assert_true(maxval(abs(band_repeat%cv_error-band_result%cv_error)) < 1.0e-18_dp, 'CV reproducibility')

   taper_result = tapering_cv(data, h=0.5_dp, n_cv=3, norm='F', seed=77)
   call assert_true(taper_result%status == fincov_ok, 'tapering CV status')
   call assert_true(size(taper_result%parameter_grid) == 5, 'tapering CV grid size')

   hard_result = threshold_cv(data, method='hard', thresh_len=8, n_cv=3, norm='F', seed=99)
   call assert_true(hard_result%status == fincov_ok, 'hard threshold CV status')
   call assert_true(size(hard_result%parameter_grid) == 8, 'hard threshold CV grid size')
   call assert_true(hard_result%parameter_index >= 1 .and. hard_result%parameter_index <= 8, 'hard threshold CV index')

   soft_result = threshold_cv(data, method='soft', thresh_len=8, n_cv=3, norm='F', seed=99)
   call assert_true(soft_result%status == fincov_ok, 'soft threshold CV status')

   operator_result = banding_cv(data, n_cv=2, norm='O', seed=55)
   call assert_true(operator_result%status == fincov_ok, 'operator-norm CV status')
   call assert_true(all(operator_result%cv_error >= 0.0_dp), 'operator-norm CV errors')

   print '(a)', 'test_cross_validation: PASS'
contains
   subroutine assert_true(condition, label)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: label
      if (.not. condition) then
         print '(a)', trim(label)//' failed'
         error stop 1
      end if
   end subroutine assert_true
end program test_cross_validation
