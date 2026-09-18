! SPDX-License-Identifier: GPL-3.0-only
program basic_kza
   use kza_api, only : dp, kz, kza, periodogram, rlv
   implicit none

   real(dp) :: signal(80)
   real(dp), allocatable :: baseline(:)
   real(dp), allocatable :: adaptive(:)
   real(dp), allocatable :: local_variance(:)
   real(dp), allocatable :: pg(:, :)
   integer :: i

   do i = 1, size(signal)
      signal(i) = sin(2.0_dp * acos(-1.0_dp) * real(i, dp) / 20.0_dp)
      if (i > 40) signal(i) = signal(i) + 0.75_dp
      signal(i) = signal(i) + 0.15_dp * sin(2.0_dp * acos(-1.0_dp) * real(7 * i, dp) / 19.0_dp)
   end do

   baseline = kz(signal, 9, k=3)
   adaptive = kza(signal, 9, y=baseline, k=3, min_size=1, impute_tails=.true.)
   local_variance = rlv(adaptive, 5)
   pg = periodogram(signal)

   print '(a,f10.5)', "KZ center value:        ", baseline(40)
   print '(a,f10.5)', "KZA center value:       ", adaptive(40)
   print '(a,f10.5)', "RLV center variance:    ", local_variance(40)
   print '(a,f10.5)', "First spectral magnitude:", pg(1, 2)

end program basic_kza
