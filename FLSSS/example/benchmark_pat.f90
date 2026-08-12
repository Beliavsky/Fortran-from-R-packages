program benchmark_pat
  use flsss_kinds, only : dp, i8
  use flsss_types, only : subset_solutions
  use flsss_fast_search, only : search_md_i8_packed, search_md_i8_pat
  implicit none
  integer, parameter :: n=80, d=8, len=8
  integer(i8) :: v(n,d), target(d), me(d)
  integer :: i, j, k, idx(len)
  real(dp) :: t0, t1, td, tp
  type(subset_solutions) :: dfs, pat

  do j=1,d
    do i=1,n
      v(i,j)=int(i*(j+2)+i*i/(20+j)+j*7,i8)
    end do
  end do
  idx=[70,71,72,73,74,75,76,77]
  target=0_i8
  do k=1,len
    target=target+v(idx(k),:)
  end do
  me=0_i8

  call cpu_time(t0)
  call search_md_i8_packed(v,len,target,me,dfs,1,30.0_dp)
  call cpu_time(t1)
  td=t1-t0
  call cpu_time(t0)
  call search_md_i8_pat(v,len,target,me,pat,1,30.0_dp)
  call cpu_time(t1)
  tp=t1-t0

  print '(a,f10.6,2x,a,i0)', 'packed DFS seconds: ',td,'states: ',dfs%nodes
  print '(a,f10.6,2x,a,i0)', 'PAT seconds:        ',tp,'states: ',pat%nodes
  print '(a,f8.2)', 'speed ratio:        ',td/max(tp,1.0e-12_dp)
end program benchmark_pat
