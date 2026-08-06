! SPDX-License-Identifier: GPL-3.0-only
program wqc_demo
   use wqc, only : dp, quantile_correlation_analysis, wqc_pair_result
   implicit none

   integer, parameter :: n = 256
   real(dp) :: x(n), y(n)
   real(dp), parameter :: quantiles(3) = [0.1_dp, 0.5_dp, 0.9_dp]
   type(wqc_pair_result) :: result
   integer :: i, j, stat
   character(len=:), allocatable :: errmsg

   do i = 1, n
      x(i) = sin(0.045_dp * real(i, dp)) + 0.3_dp * cos(0.24_dp * real(i, dp))
      y(i) = 0.7_dp * x(i) + 0.3_dp * sin(0.17_dp * real(i, dp) + 0.5_dp)
   end do

   allocate(character(len=1) :: errmsg)
   errmsg = ''
   result = quantile_correlation_analysis(x, y, quantiles, wf='la8', j_levels=4, &
      n_sim=250, seed=20250612, stat=stat, errmsg=errmsg)
   if (stat /= 0) then
      write(*, '(a)') 'WQC failed: '//errmsg
      error stop 1
   end if

   write(*, '(a)') ' level  quantile     estimated      ci_lower      ci_upper'
   do j = 1, result%levels
      do i = 1, size(result%quantiles)
         write(*, '(i6,1x,f9.3,3(1x,f13.6))') j, result%quantiles(i), &
            result%estimated_qc(j, i), result%ci_lower(j, i), result%ci_upper(j, i)
      end do
   end do
end program wqc_demo
