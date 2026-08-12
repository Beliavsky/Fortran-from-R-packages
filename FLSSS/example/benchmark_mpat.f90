program benchmark_mpat
  use flsss_kinds, only : dp, i8
  use flsss_types, only : subset_solutions
  use flsss_fast_search, only : search_md_i8_packed, search_md_i8_pat
  use flsss_mpat, only : search_md_i8_mpat
  implicit none
  integer, parameter :: n=45, d=24, len=10
  integer(i8) :: v(n,d), target(d), me(d)
  integer :: i, j, k, idx(len)
  real(dp) :: t0, t1, tp, tt, tm
  type(subset_solutions) :: packed, pat, mpat

  do j=1,d
    do i=1,n
      v(i,j)=int(i*(2+mod(j,5))+i*i/(18+mod(j,9))+j*3,i8)
    end do
  end do
  idx=[5,9,13,17,21,25,29,33,37,41]
  target=0_i8
  do k=1,len
    target=target+v(idx(k),:)
  end do
  me=0_i8

  call cpu_time(t0)
  call search_md_i8_packed(v,len,target,me,packed,1,30.0_dp)
  call cpu_time(t1)
  tp=t1-t0

  call cpu_time(t0)
  call search_md_i8_pat(v,len,target,me,pat,1,30.0_dp)
  call cpu_time(t1)
  tt=t1-t0

  call cpu_time(t0)
  call search_md_i8_mpat(v,len,target,me,mpat,1,30.0_dp)
  call cpu_time(t1)
  tm=t1-t0

  print '(a,f10.6,2x,a,i0)', 'packed DFS seconds: ',tp,'states: ',packed%nodes
  print '(a,f10.6,2x,a,i0)', 'v0.3 PAT seconds:   ',tt,'states: ',pat%nodes
  print '(a,f10.6,2x,a,i0)', 'v0.4 mPAT seconds:  ',tm,'states: ',mpat%nodes
  print '(a,i0,2x,a,i0)', 'triangular entries: ',mpat%tri_entries,'lookups: ',mpat%tri_lookups
  print '(a,f8.2)', 'packed/mPAT ratio:  ',tp/max(tm,1.0e-12_dp)
  print '(a,f8.2)', 'PAT/mPAT ratio:     ',tt/max(tm,1.0e-12_dp)
end program benchmark_mpat
