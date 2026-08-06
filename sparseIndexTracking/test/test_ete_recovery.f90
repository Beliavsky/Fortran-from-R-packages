program test_ete_recovery
   use sparse_index_tracking, only : dp, sparse_index_fit, fit_sparse_index_tracking, &
                                     spIndexTrack, sit_success
   implicit none

   integer, parameter :: m = 120, n = 6
   real(dp) :: x(m, n), index_returns(m)
   real(dp), allocatable :: weights(:), alias_weights(:)
   type(sparse_index_fit) :: fit
   integer :: info
   character(len=256) :: message

   call make_data(x, index_returns)

   call fit_sparse_index_tracking(x, index_returns, 1.0e-7_dp, fit, &
                                  upper_bound=0.8_dp, measure='ete')
   call assert_true(fit%info == sit_success, 'typed ETE fit failed: ' // trim(fit%message))
   call assert_true(fit%converged, 'typed ETE fit did not converge')
   call assert_close(sum(fit%weights), 1.0_dp, 2.0e-12_dp, 'weight sum')
   call assert_true(all(fit%weights >= 0.0_dp) .and. all(fit%weights <= 0.8_dp), 'weight bounds')
   call assert_true(abs(fit%weights(1) - 0.7_dp) < 2.0e-3_dp, 'first weight recovery')
   call assert_true(abs(fit%weights(2) - 0.3_dp) < 2.0e-3_dp, 'second weight recovery')
   call assert_true(maxval(abs(fit%weights(3:n))) < 1.0e-10_dp, 'irrelevant assets were not removed')
   call assert_true(fit%cardinality == 2, 'unexpected ETE cardinality')

   call spIndexTrack(x, index_returns, 1.0e-7_dp, weights, info, message, &
                     upper_bound=0.8_dp, measure='ETE')
   call assert_true(info == sit_success, 'compatibility ETE fit failed: ' // trim(message))
   call assert_close(maxval(abs(weights - fit%weights)), 0.0_dp, 1.0e-12_dp, 'compatibility weights')

   call spIndexTrack(x, index_returns, 1.0e-7_dp, alias_weights, info, message, &
                     upper_bound=0.8_dp, measure='ete', threshold=0.31_dp)
   call assert_true(info == sit_success, 'thresholded fit failed')
   call assert_true(count(alias_weights > 0.0_dp) == 1, 'threshold did not reduce cardinality')
   call assert_close(alias_weights(1), 1.0_dp, 2.0e-12_dp, 'threshold renormalization')

   print '(a)', 'test_ete_recovery: PASS'

contains

   subroutine make_data(data, target)
      real(dp), intent(out) :: data(:, :), target(:)
      integer :: row
      real(dp) :: time

      do row = 1, size(data, 1)
         time = real(row, dp)
         data(row, 1) = 0.010_dp * sin(0.13_dp * time) + 0.002_dp * cos(0.07_dp * time)
         data(row, 2) = 0.008_dp * cos(0.11_dp * time) - 0.001_dp * sin(0.03_dp * time)
         data(row, 3) = 0.006_dp * sin(0.17_dp * time + 0.4_dp)
         data(row, 4) = 0.007_dp * cos(0.19_dp * time + 0.2_dp)
         data(row, 5) = 0.005_dp * sin(0.23_dp * time) + 0.003_dp * cos(0.05_dp * time)
         data(row, 6) = 0.004_dp * cos(0.29_dp * time) - 0.002_dp * sin(0.09_dp * time)
      end do
      target = 0.7_dp * data(:, 1) + 0.3_dp * data(:, 2)
   end subroutine make_data


   subroutine assert_true(condition, message)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: message

      if (.not. condition) then
         print '(a)', 'FAIL: ' // trim(message)
         error stop 1
      end if
   end subroutine assert_true


   subroutine assert_close(actual, expected_value, tolerance, message)
      real(dp), intent(in) :: actual, expected_value, tolerance
      character(len=*), intent(in) :: message

      call assert_true(abs(actual - expected_value) <= tolerance, message)
   end subroutine assert_close

end program test_ete_recovery
