! SPDX-License-Identifier: MIT
program test_inference
   use vasicekfit, only : dp, vasicek_fit_result, covariance_result, &
      confidence_interval_result, fit_vasicek, vasicek_covariance, &
      vasicek_confidence_intervals
   implicit none
   real(dp), parameter :: y(20) = [ &
      0.009600134840871976_dp, 0.029502890965674730_dp, 0.037198232908477340_dp, &
      0.020785442980663090_dp, 0.050685547743244040_dp, 0.017640200610933460_dp, &
      0.032230717502823025_dp, 0.072262483067218770_dp, 0.016096229432628024_dp, &
      0.048090550856388070_dp, 0.031519799988685336_dp, 0.065329288703111180_dp, &
      0.015700650584264047_dp, 0.060401782241646530_dp, 0.039279530963784316_dp, &
      0.033014580456703310_dp, 0.081703730226032780_dp, 0.045890279005243326_dp, &
      0.071152806570663060_dp, 0.033199582822105970_dp ]
   real(dp), parameter :: covariance_ref(4,4) = reshape([ &
      2.162765901062822e-05_dp, 1.382545612941461e-05_dp, -2.989830290235622e-07_dp, 1.385924403179396e-05_dp, &
      1.382545612941461e-05_dp, 1.754521830589704e-04_dp, -1.181251841231018e-05_dp, 2.289572557427135e-06_dp, &
     -2.989830290235619e-07_dp,-1.181251841231018e-05_dp,  3.409060525317405e-03_dp, 4.924406692691221e-04_dp, &
      1.385924403179396e-05_dp, 2.289572557427135e-06_dp,  4.924406692691221e-04_dp, 1.066445375875629e-02_dp ], [4,4])
   real(dp) :: x(20,2)
   integer :: i
   type(vasicek_fit_result) :: fit
   type(covariance_result) :: iid, hac
   type(confidence_interval_result) :: ci95, ci99

   do i = 1, 20
      x(i,1) = -1.5_dp + 3.0_dp * real(i-1,dp) / 19.0_dp
      x(i,2) = (-1.0_dp)**(i-1) * (0.2_dp + 0.03_dp * real(i-1,dp))
   end do
   fit = fit_vasicek(y, x)
   iid = vasicek_covariance(fit)
   hac = vasicek_covariance(fit, use_hac=.true., hac_lag=3)
   ci95 = vasicek_confidence_intervals(fit, 0.95_dp)
   ci99 = vasicek_confidence_intervals(fit, 0.99_dp)

   call assert_true(iid%ok .and. hac%ok .and. ci95%ok .and. ci99%ok)
   call assert_close(maxval(abs(iid%covariance - covariance_ref)), 0.0_dp, 3.0e-12_dp)
   call assert_close(maxval(abs(iid%covariance - transpose(iid%covariance))), 0.0_dp, 0.0_dp)
   call assert_close(maxval(abs(hac%covariance - transpose(hac%covariance))), 0.0_dp, 0.0_dp)
   call assert_true(all(diagonal(hac%covariance) >= 0.0_dp))
   call assert_true(all(ci95%lower <= ci95%estimate) .and. all(ci95%estimate <= ci95%upper))
   call assert_true(all((ci99%upper - ci99%lower) > (ci95%upper - ci95%lower)))

   print '(a)', 'test_inference: PASS'

contains

   function diagonal(a) result(d)
      real(dp), intent(in) :: a(:,:)
      real(dp), allocatable :: d(:)
      integer :: j
      allocate(d(min(size(a,1),size(a,2))))
      do j = 1, size(d)
         d(j) = a(j,j)
      end do
   end function diagonal

   subroutine assert_close(actual, expected, tolerance)
      real(dp), intent(in) :: actual, expected, tolerance
      if (abs(actual - expected) > tolerance + tolerance * abs(expected)) then
         print '(a,3es25.16)', 'mismatch: ', actual, expected, abs(actual - expected)
         error stop 1
      end if
   end subroutine assert_close

   subroutine assert_true(condition)
      logical, intent(in) :: condition
      if (.not. condition) error stop 1
   end subroutine assert_true

end program test_inference
