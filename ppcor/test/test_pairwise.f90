program test_pairwise
   use ppcor, only : dp, ppcor_result, ppcor_test_result, pcor, spcor, &
                     pcor_test, spcor_test, ppcor_pearson, ppcor_success
   implicit none
   real(dp) :: x(30,4)
   type(ppcor_result) :: full
   type(ppcor_test_result) :: one

   call make_data(x)
   call pcor(x, full, ppcor_pearson)
   call pcor_test(x(:,1), x(:,2), x(:,3:4), one, ppcor_pearson)
   call assert_true(one%status == ppcor_success, 'matrix-z pcor.test status')
   call assert_close(one%estimate, full%estimate(1,2), 2.0e-12_dp, &
                     'matrix-z pcor.test estimate')
   call assert_close(one%p_value, full%p_value(1,2), 2.0e-12_dp, &
                     'matrix-z pcor.test p-value')

   call pcor(x(:,1:3), full, ppcor_pearson)
   call pcor_test(x(:,1), x(:,2), x(:,3), one, ppcor_pearson)
   call assert_close(one%estimate, full%estimate(1,2), 2.0e-12_dp, &
                     'vector-z pcor.test estimate')

   call spcor(x, full, ppcor_pearson)
   call spcor_test(x(:,1), x(:,2), x(:,3:4), one, ppcor_pearson)
   call assert_close(one%estimate, full%estimate(1,2), 2.0e-12_dp, &
                     'matrix-z spcor.test estimate')
   call assert_true(one%n == 30 .and. one%gp == 2, 'pairwise metadata')

   print '(a)', 'test_pairwise: PASS'

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

end program test_pairwise
