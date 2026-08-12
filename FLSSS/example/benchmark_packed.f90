program benchmark_packed
  use flsss_kinds, only : dp, i8
  use flsss_types, only : subset_solutions
  use flsss_search, only : search_md_i8
  use flsss_fast_search, only : search_md_i8_packed
  implicit none
  integer, parameter :: n=80, d=24, len=8
  integer(i8) :: v(n,d), target(d), me(d)
  integer :: i, j, k, idx(len)
  real(dp) :: t0, t1, td, tp
  type(subset_solutions) :: ref, packed

  do j=1,d
    do i=1,n
      v(i,j)=int(i*(2+mod(j,3))+i*i/(20+mod(j,7))+j*7,i8)
    end do
  end do
  idx=[7,13,21,32,44,57,68,77]
  target=0_i8
  do k=1,len
    target=target+v(idx(k),:)
  end do
  me=0_i8

  call cpu_time(t0)
  call search_md_i8(v,len,target,me,ref,1,30.0_dp)
  call cpu_time(t1)
  td=t1-t0
  call cpu_time(t0)
  call search_md_i8_packed(v,len,target,me,packed,1,30.0_dp)
  call cpu_time(t1)
  tp=t1-t0

  print '(a,f10.6,2x,a,i0)', 'v0.2 DFS seconds: ',td,'nodes: ',ref%nodes
  print '(a,f10.6,2x,a,i0)', 'packed seconds:   ',tp,'nodes: ',packed%nodes
  print '(a,i0)', 'packed lanes:     ',packed%packed_lanes
  print '(a,f8.2)', 'speed ratio:      ',td/max(tp,1.0e-12_dp)
end program benchmark_packed
