program test_rank_methods
   use ppcor, only : dp, ppcor_result, pcor, spcor, ppcor_spearman, &
                     ppcor_kendall, ppcor_success
   implicit none
   real(dp) :: x(30,4), tied(12,4)
   type(ppcor_result) :: r

   call make_data(x)
   call pcor(x, r, ppcor_spearman)
   call assert_true(r%status == ppcor_success, 'Spearman status')
   call assert_close(r%estimate(1,2), 0.3623558901908250_dp, 3.0e-11_dp, &
                     'Spearman partial')
   call assert_close(r%p_value(1,2), 0.05809742714549789_dp, 3.0e-10_dp, &
                     'Spearman p value')
   call spcor(x, r, ppcor_spearman)
   call assert_close(r%estimate(1,2), 0.2639695792881723_dp, 3.0e-11_dp, &
                     'Spearman semi-partial')

   call pcor(x, r, ppcor_kendall)
   call assert_true(r%status == ppcor_success, 'Kendall status')
   call assert_close(r%estimate(1,2), 0.24556553978755083_dp, 4.0e-11_dp, &
                     'Kendall partial')
   call assert_close(r%statistic(1,2), 1.8338747011677814_dp, 4.0e-10_dp, &
                     'Kendall z statistic')
   call assert_close(r%p_value(1,2), 0.06667259574123358_dp, 4.0e-10_dp, &
                     'Kendall p value')

   call make_tied_data(tied)
   call pcor(tied, r, ppcor_spearman)
   call assert_close(r%estimate(1,2), 0.12067041220203246_dp, 5.0e-11_dp, &
                     'tied Spearman partial')
   call pcor(tied, r, ppcor_kendall)
   call assert_close(r%estimate(1,2), 0.18247971747205952_dp, 5.0e-11_dp, &
                     'tied Kendall tau-b partial')

   print '(a)', 'test_rank_methods: PASS'

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

   subroutine make_tied_data(a)
      real(dp), intent(out) :: a(:,:)
      a(:,1) = real([1,1,2,2,3,3,4,5,5,6,7,8],dp)
      a(:,2) = real([4,3,3,2,2,1,1,2,4,5,6,6],dp)
      a(:,3) = real([2,1,2,3,3,4,4,5,6,5,7,8],dp)
      a(:,4) = real([1,2,2,3,4,4,5,5,6,7,7,8],dp)
   end subroutine make_tied_data

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

end program test_rank_methods
