program test_pearson
   use ppcor, only : dp, ppcor_result, pcor, spcor, ppcor_pearson, ppcor_success
   implicit none
   real(dp) :: x(30,4)
   type(ppcor_result) :: r
   integer :: i

   call make_data(x)
   call pcor(x, r, ppcor_pearson)
   call assert_true(r%status == ppcor_success, 'Pearson pcor status')
   call assert_close(r%estimate(1,2), 0.36423313554331355_dp, 2.0e-11_dp, &
                     'Pearson partial correlation')
   call assert_close(r%statistic(1,2), 1.9942191622079826_dp, 2.0e-10_dp, &
                     'Pearson t statistic')
   call assert_close(r%p_value(1,2), 0.05671348995699716_dp, 2.0e-10_dp, &
                     'Pearson p value')
   call assert_true(maxval(abs(r%estimate-transpose(r%estimate))) < 1.0e-12_dp, &
                    'partial correlation symmetry')
   do i = 1, 4
      call assert_close(r%estimate(i,i), 1.0_dp, 0.0_dp, 'partial diagonal')
      call assert_close(r%statistic(i,i), 0.0_dp, 0.0_dp, 'statistic diagonal')
      call assert_close(r%p_value(i,i), 0.0_dp, 0.0_dp, 'p-value diagonal')
   end do

   call spcor(x, r, ppcor_pearson)
   call assert_true(r%status == ppcor_success, 'Pearson spcor status')
   call assert_close(r%estimate(1,2), 0.24877665718245207_dp, 2.0e-11_dp, &
                     'Pearson semi-partial 1,2')
   call assert_close(r%estimate(2,1), 0.21816653633311450_dp, 2.0e-11_dp, &
                     'Pearson semi-partial 2,1')
   call assert_true(abs(r%estimate(1,2)-r%estimate(2,1)) > 1.0e-3_dp, &
                    'semi-partial matrix is directional')

   print '(a)', 'test_pearson: PASS'

contains

   subroutine make_data(a)
      real(dp), intent(out) :: a(:,:)
      integer :: k
      real(dp) :: t, z1, z2, xx, yy
      do k = 1, size(a,1)
         t = real(k,dp)
         z1 = sin(0.3_dp*t)
         z2 = cos(0.17_dp*t) + 0.03_dp*t
         xx = 0.8_dp*z1 - 0.4_dp*z2 + 0.5_dp*sin(1.1_dp*t)
         yy = 0.2_dp*z1 + 0.6_dp*z2 + 0.45_dp*xx + 0.6_dp*cos(0.7_dp*t)
         a(k,:) = [xx, yy, z1, z2]
      end do
   end subroutine make_data

   subroutine assert_close(actual, expected, tol, label)
      real(dp), intent(in) :: actual, expected, tol
      character(len=*), intent(in) :: label
      if (abs(actual-expected) > tol) then
         write(*,'(a,2(1x,es24.16))') trim(label)//' failed:', actual, expected
         error stop 1
      end if
   end subroutine assert_close

   subroutine assert_true(condition, label)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: label
      if (.not. condition) then
         write(*,'(a)') trim(label)//' failed'
         error stop 1
      end if
   end subroutine assert_true

end program test_pearson
