program basic_subset_sum
  use flsss_api
  implicit none
  type(subset_solutions) :: r
  real(dp) :: v(10)
  integer :: i
  v=[(real(i,dp),i=1,10)]
  r=flsss(3,v,target=15.0_dp,me=0.0_dp,solution_need=5)
  print '(a,i0)', 'solutions found: ', r%size()
  do i=1,r%size()
    print '(*(i0,1x))',r%sol(i)%idx
  end do
end program basic_subset_sum
