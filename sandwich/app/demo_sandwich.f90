! SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only
program demo_sandwich
   use sandwich, only : dp, ols_model, fit_ols, vcov_hc, vcov_hac, &
      vcov_cluster, newey_west_weights, SANDWICH_SUCCESS
   implicit none

   integer, parameter :: n = 12
   real(dp) :: x(n, 2), y(n)
   integer :: cluster(n, 1), i, status
   real(dp), allocatable :: hc3(:, :), hac(:, :), clustered(:, :), lag_weights(:)
   type(ols_model) :: model

   do i = 1, n
      x(i, 1) = 1.0_dp
      x(i, 2) = real(i - 1, dp)
      y(i) = 1.5_dp + 0.35_dp * x(i, 2) + 0.25_dp * sin(0.9_dp * real(i, dp))
      cluster(i, 1) = (i + 2) / 3
   end do

   call fit_ols(x, y, model, status)
   if (status /= SANDWICH_SUCCESS) error stop 'OLS fit failed'

   call vcov_hc(x, model%residuals, model%bread, 'HC3', hc3, status, model%hat)
   if (status /= SANDWICH_SUCCESS) error stop 'HC3 failed'

   call newey_west_weights(model%scores, lag_weights, status, lag = 2, prewhite_order = 0)
   if (status /= SANDWICH_SUCCESS) error stop 'Newey-West weights failed'
   call vcov_hac(model%scores, model%bread, lag_weights, hac, status, &
      adjust = .false., prewhite_order = 0)
   if (status /= SANDWICH_SUCCESS) error stop 'HAC failed'

   call vcov_cluster(model%scores, cluster, model%bread, clustered, status, &
      type = 'HC1', cadjust = .true.)
   if (status /= SANDWICH_SUCCESS) error stop 'cluster covariance failed'

   print '(a,2f12.6)', 'coefficients:       ', model%coefficients
   print '(a,2f12.6)', 'classical std. err.:', sqrt(diagonal(model%covariance))
   print '(a,2f12.6)', 'HC3 std. err.:      ', sqrt(diagonal(hc3))
   print '(a,2f12.6)', 'HAC std. err.:      ', sqrt(diagonal(hac))
   print '(a,2f12.6)', 'cluster std. err.:  ', sqrt(diagonal(clustered))

contains

   pure function diagonal(a) result(d)
      real(dp), intent(in) :: a(:, :)
      real(dp) :: d(min(size(a, 1), size(a, 2)))
      integer :: j
      do j = 1, size(d)
         d(j) = a(j, j)
      end do
   end function diagonal

end program demo_sandwich
