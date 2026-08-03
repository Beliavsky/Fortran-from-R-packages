! SPDX-License-Identifier: GPL-2.0-or-later
program demo_icsnp
   use icsnp, only : dp, icsnp_ok, spatial_median, tyler_shape, HP_loc_test, test_result
   implicit none
   integer, parameter :: n = 80, p = 3
   real(dp) :: x(n, p), zero(p)
   real(dp), allocatable :: center(:), shape(:,:)
   type(test_result) :: test
   integer :: i, status, iterations

   do i = 1, n
      x(i, 1) = sin(0.19_dp * real(i, dp))
      x(i, 2) = 0.5_dp * x(i, 1) + cos(0.27_dp * real(i, dp))
      x(i, 3) = -0.2_dp * x(i, 1) + sin(0.37_dp * real(i, dp))
   end do
   zero = 0.0_dp

   call spatial_median(x, center, status, iterations)
   if (status /= icsnp_ok) error stop 'spatial median failed'
   call tyler_shape(x, shape, status, iterations, location=center)
   if (status /= icsnp_ok) error stop 'Tyler shape failed'
   call HP_loc_test(x, test, mu=zero, scores='rank')
   if (test%status /= icsnp_ok) error stop 'HP location test failed'

   write(*, '(a,3f11.5)') 'Spatial median: ', center
   write(*, '(a,f11.5)') 'Tyler shape trace: ', shape(1,1) + shape(2,2) + shape(3,3)
   write(*, '(a,f11.5,a,f11.5)') 'HP rank statistic: ', test%statistic, &
      ', p-value: ', test%p_value
end program demo_icsnp
