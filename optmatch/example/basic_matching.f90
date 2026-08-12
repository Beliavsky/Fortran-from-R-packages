program basic_matching
   use optmatch, only : dp, distance_spec, match_result
   use optmatch, only : dense_distance, fullmatch, effective_sample_size
   implicit none
   type(distance_spec) :: d
   type(match_result) :: m
   integer :: i

   ! Rows are treatment units; columns are controls.
   d = dense_distance(reshape([ &
      1.0_dp, 4.0_dp, &
      3.0_dp, 1.0_dp, &
      2.0_dp, 5.0_dp  &
   ], [2,3]))

   m = fullmatch(d, min_controls=1.0_dp, max_controls=2.0_dp)

   write(*,'(a,l1)') 'feasible: ', m%feasible
   write(*,'(a,f8.3)') 'objective: ', m%objective
   write(*,'(a,f8.3)') 'effective sample size: ', effective_sample_size(m)
   do i = 1, size(m%treatment_group)
      write(*,'(a,i0,a,i0)') 'treatment ', i, ' -> set ', m%treatment_group(i)
   end do
   do i = 1, size(m%control_group)
      write(*,'(a,i0,a,i0)') 'control   ', i, ' -> set ', m%control_group(i)
   end do
end program basic_matching
