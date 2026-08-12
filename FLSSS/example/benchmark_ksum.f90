program benchmark_ksum
  use flsss_mod
  implicit none
  character(len=32) :: v(20,2), target(2)
  type(subset_solutions) :: plain, fast
  type(ksum_table) :: tab
  integer(i8) :: c0,c1,rate
  integer :: i
  real(dp) :: tplain,tfast

  v='0'
  do i=1,20
    write(v(i,1),'(i0)') i
    write(v(i,2),'(i0)') i*i
  end do
  target=[character(len=32) :: '132','2220']
  tab=build_ksum_hash(4,v,max_entries=10000)
  call system_clock(c0,rate)
  plain=arb_flsss(8,v,target,solution_need=1)
  call system_clock(c1)
  tplain=real(c1-c0,dp)/real(rate,dp)
  call system_clock(c0)
  fast=arb_flsss(8,v,target,solution_need=1,given_ksum=tab)
  call system_clock(c1)
  tfast=real(c1-c0,dp)/real(rate,dp)
  print '(a,f10.6,a,i0)', 'plain seconds: ',tplain,' nodes: ',plain%nodes
  print '(a,f10.6,a,i0,a,i0)', 'hash seconds: ',tfast,' nodes: ',fast%nodes, &
    ' lookups: ',fast%hash_lookups
  if(tfast>0.0_dp) print '(a,f10.2)', 'speed ratio: ',tplain/tfast
end program benchmark_ksum
