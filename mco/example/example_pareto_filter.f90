! SPDX-License-Identifier: GPL-2.0-only
program example_pareto_filter
   use mco, only : dp, pareto_filter
   implicit none
   real(dp) :: points(2,5)
   real(dp),allocatable :: front(:,:)
   integer :: i
   points=reshape([1.0_dp,4.0_dp,2.0_dp,3.0_dp,3.0_dp,2.0_dp,4.0_dp,1.0_dp,4.0_dp,4.0_dp],[2,5])
   front=pareto_filter(points)
   do i=1,size(front,2); print '(2f10.4)',front(:,i); end do
end program
