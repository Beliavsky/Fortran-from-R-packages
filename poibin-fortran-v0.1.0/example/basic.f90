! SPDX-License-Identifier: GPL-2.0-only
program basic
   use poibin, only : dp, dpoibin_vec, ppoibin_vec, qpoibin_vec, rpoibin_sample, poibin_seed
   implicit none
   real(dp), parameter :: pp(5) = [0.1_dp, 0.2_dp, 0.3_dp, 0.4_dp, 0.5_dp]
   integer, parameter :: wts(5) = [2, 2, 2, 2, 2]
   integer :: k(11), qx(5), draws(5), i
   real(dp) :: pmf(11), cdf(11)
   real(dp), parameter :: q(5) = [0.1_dp, 0.25_dp, 0.5_dp, 0.75_dp, 0.9_dp]

   k = [(i, i=0,10)]
   call dpoibin_vec(k, pp, pmf, wts)
   call ppoibin_vec(k, pp, cdf, 'DFT-CF', wts)
   call qpoibin_vec(q, pp, qx, wts)
   call poibin_seed(20260816)
   call rpoibin_sample(5, pp, draws, wts)

   print '(a)', ' k          pmf          cdf'
   do i = 1, size(k)
      print '(i2,2(2x,f12.8))', k(i), pmf(i), cdf(i)
   end do
   print '(a,5(1x,i0))', 'quantiles:', qx
   print '(a,5(1x,i0))', 'draws:    ', draws
end program basic
