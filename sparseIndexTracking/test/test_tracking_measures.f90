program test_tracking_measures
   use sparse_index_tracking, only : dp, sp_index_track, sit_success
   implicit none

   integer, parameter :: m = 120, n = 6
   character(len=4), parameter :: measures(4) = [character(len=4) :: 'ete ', 'dr  ', 'hete', 'hdr ']
   real(dp) :: x(m, n), contaminated(m), clean_target(m), t
   real(dp), allocatable :: weights(:)
   real(dp) :: ete_weights(n), hete_weights(n)
   integer :: i, info, j
   character(len=256) :: message

   do i = 1, m
      t = real(i, dp)
      x(i, 1) = 0.010_dp * sin(0.13_dp * t) + 0.002_dp * cos(0.07_dp * t)
      x(i, 2) = 0.008_dp * cos(0.11_dp * t) - 0.001_dp * sin(0.03_dp * t)
      x(i, 3) = 0.006_dp * sin(0.17_dp * t + 0.4_dp)
      x(i, 4) = 0.007_dp * cos(0.19_dp * t + 0.2_dp)
      x(i, 5) = 0.005_dp * sin(0.23_dp * t) + 0.003_dp * cos(0.05_dp * t)
      x(i, 6) = 0.004_dp * cos(0.29_dp * t) - 0.002_dp * sin(0.09_dp * t)
   end do
   clean_target = 0.7_dp * x(:, 1) + 0.3_dp * x(:, 2)
   contaminated = clean_target
   contaminated(30) = contaminated(30) + 0.2_dp

   do j = 1, size(measures)
      call sp_index_track(x, contaminated, 1.0e-6_dp, weights, info, message, &
                          upper_bound=0.8_dp, measure=trim(measures(j)), huber_parameter=0.02_dp)
      call assert_true(info == sit_success, trim(measures(j)) // ' fit failed: ' // trim(message))
      call assert_true(abs(sum(weights) - 1.0_dp) < 2.0e-11_dp, trim(measures(j)) // ' sum')
      call assert_true(all(weights >= 0.0_dp), trim(measures(j)) // ' nonnegativity')
      call assert_true(all(weights <= 0.8_dp + 2.0e-12_dp), trim(measures(j)) // ' cap')
      if (trim(measures(j)) == 'ete') then
         ete_weights = weights
      else if (trim(measures(j)) == 'hete') then
         hete_weights = weights
      end if
   end do

   call assert_true(norm2(hete_weights - [0.7_dp, 0.3_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp]) < &
                    norm2(ete_weights - [0.7_dp, 0.3_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp]), &
                    'Huber ETE was not more robust to the seeded outlier')
   call assert_true(norm2(matmul(x, hete_weights) - clean_target) < &
                    norm2(matmul(x, ete_weights) - clean_target), &
                    'Huber ETE did not improve clean-target tracking')

   print '(a)', 'test_tracking_measures: PASS'

contains

   subroutine assert_true(condition, message)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: message

      if (.not. condition) then
         print '(a)', 'FAIL: ' // trim(message)
         error stop 1
      end if
   end subroutine assert_true

end program test_tracking_measures
