! SPDX-License-Identifier: GPL-2.0-or-later
program demo_moments
   use moments, only : dp, all_moments, all_cumulants, skewness, kurtosis, jarque_test, &
      moments_test_result
   implicit none

   real(dp) :: returns(12)
   real(dp), allocatable :: raw(:), cumulants(:)
   type(moments_test_result) :: jb

   returns = [0.012_dp, -0.008_dp, 0.004_dp, 0.019_dp, -0.015_dp, 0.006_dp, &
      0.003_dp, -0.004_dp, 0.027_dp, -0.011_dp, 0.008_dp, 0.002_dp]
   raw = all_moments(returns, 4)
   cumulants = all_cumulants(raw)
   jb = jarque_test(returns)

   write(*, '(a,f12.8)') 'mean: ', raw(2)
   write(*, '(a,f12.8)') 'variance cumulant: ', cumulants(3)
   write(*, '(a,f12.8)') 'skewness: ', skewness(returns)
   write(*, '(a,f12.8)') 'Pearson kurtosis: ', kurtosis(returns)
   write(*, '(a,f12.8)') 'Jarque-Bera p-value: ', jb%p_value
end program demo_moments
