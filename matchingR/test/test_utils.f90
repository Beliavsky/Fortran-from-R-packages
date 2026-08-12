program test_utils
   use matchingr
   implicit none
   real(dp) :: u(4,3)
   integer, allocatable :: s(:,:), r(:,:)
   u=reshape([3.0_dp,2.0_dp,8.0_dp,1.0_dp, &
              12.0_dp,2.0_dp,9.0_dp,2.0_dp, &
              13.0_dp,5.0_dp,3.1_dp,2.1_dp],[4,3])
   s=sort_index(u)
   if(any(s /= reshape([3,1,2,4, 1,3,2,4, 1,2,3,4],[4,3]))) error stop "sort_index"
   r=rank_index(s)
   if(any(sort_index(-real(r,dp)) /= s)) error stop "rank inverse"
   print *, "test_utils: PASS"
end program
