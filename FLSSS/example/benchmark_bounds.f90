program benchmark_bounds
  use flsss_mod
  implicit none
  type(subset_solutions) :: r
  real(dp) :: v(26,3), target(3), me(3)
  integer(i8) :: c0,c1,rate
  integer :: i
  real(dp) :: seconds

  do i=1,26
    v(i,1)=real(mod(17*i+11,97),dp)
    v(i,2)=real(mod(29*i+7,113),dp)
    v(i,3)=real(mod(41*i+3,127),dp)
  end do
  target=10000.0_dp
  me=0.0_dp
  call system_clock(c0,rate)
  r=mflsss_par(8,v,target,me,solution_need=1)
  call system_clock(c1)
  seconds=real(c1-c0,dp)/real(rate,dp)
  print '(a,f10.6)', 'elapsed seconds: ',seconds
  print '(a,i0)', 'nodes: ',r%nodes
  print '(a,i0)', 'pruned states: ',r%pruned
end program benchmark_bounds
