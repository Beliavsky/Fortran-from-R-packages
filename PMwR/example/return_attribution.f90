program return_attribution_example
   use pmwr, only : dp, attribution_result, return_attribution
   implicit none
   type(attribution_result) :: result
   real(dp) :: r(3,2), w(3,2), rb(3,2), wb(3,2)

   r(1,:) = [0.03_dp, 0.01_dp]
   r(2,:) = [-0.02_dp, 0.02_dp]
   r(3,:) = [0.01_dp, 0.04_dp]
   rb(1,:) = [0.02_dp, 0.00_dp]
   rb(2,:) = [-0.01_dp, 0.01_dp]
   rb(3,:) = [0.00_dp, 0.03_dp]
   w = spread([0.60_dp, 0.40_dp], dim=1, ncopies=3)
   wb = 0.50_dp

   call return_attribution(r, w, rb, wb, result)
   print '(a,3(f10.6,1x))', 'Linked allocation/selection/interaction: ', result%linked_total
end program return_attribution_example
