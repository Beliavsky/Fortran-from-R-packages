program pcoa_example
   use ape, only : dp, pcoa, pcoa_cailliez, pcoa_result
   implicit none

   real(dp) :: distance(4, 4)
   type(pcoa_result) :: result
   integer :: info

   distance = reshape([ &
      0.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, &
      1.0_dp, 0.0_dp, 1.0_dp, 3.0_dp, &
      1.0_dp, 1.0_dp, 0.0_dp, 3.0_dp, &
      1.0_dp, 3.0_dp, 3.0_dp, 0.0_dp], [4, 4])

   call pcoa(distance, result, info, 'cailliez')
   if (info /= 0) error stop 'PCoA failed'
   if (result%correction /= pcoa_cailliez) error stop 'unexpected correction mode'

   print '(a,f12.7)', 'Cailliez constant: ', result%correction_constant
   print '(a,*(1x,f12.7))', 'Corrected eigenvalues:', result%corrected_eigenvalues
   print '(a)', 'First two corrected principal coordinates:'
   print '(*(1x,f12.7))', result%corrected_vectors(:, 1)
   print '(*(1x,f12.7))', result%corrected_vectors(:, 2)
end program pcoa_example
