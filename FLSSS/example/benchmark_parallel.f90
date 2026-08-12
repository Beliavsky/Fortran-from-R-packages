program benchmark_parallel
  use flsss_mod
  implicit none
  integer, parameter :: n=65, d=6, len=8
  real(dp) :: v(n,d), target(d), me(d), ts, tp
  integer :: i, j, k, idx(len)
  integer(i8) :: c0, c1, rate
  type(mflsss_decomposition) :: dec
  type(subset_solutions) :: serial, parallel

  do j=1,d
    do i=1,n
      v(i,j)=real(i*(j+2)+mod(i*i+3*j,23),dp)
    end do
  end do
  idx=[5,12,20,29,37,46,55,63]
  target=0.0_dp
  do k=1,len
    target=target+v(idx(k),:)
  end do
  me=0.0_dp
  dec=decompose_mflsss(len,v,target,me,approx_ninstance=24)

  call system_clock(c0,rate)
  serial=mflsss_decomp_run(dec,10000,30.0_dp,.false.,4)
  call system_clock(c1)
  ts=real(c1-c0,dp)/real(rate,dp)
  call system_clock(c0)
  parallel=mflsss_decomp_run(dec,10000,30.0_dp,.true.,4)
  call system_clock(c1)
  tp=real(c1-c0,dp)/real(rate,dp)

  print '(a,f10.5,2x,a,i0)', 'serial wall seconds:   ',ts,'nodes: ',serial%nodes
  print '(a,f10.5,2x,a,i0)', 'parallel wall seconds: ',tp,'nodes: ',parallel%nodes
  print '(a,a)', 'parallel engine:       ',trim(parallel%engine)
  print '(a,f8.2)', 'wall speed ratio:      ',ts/max(tp,1.0e-12_dp)
end program benchmark_parallel
