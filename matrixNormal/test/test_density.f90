program test_density
  use matrixNormal
  use test_support
  implicit none
  real(dp)::a(2,2),m(2,2),u(2,2),v(2,2),ld,d
  logical::ok
  character(len=256)::msg
  a=reshape([0.4_dp,1.2_dp,-0.7_dp,0.3_dp],[2,2])
  m=reshape([0.1_dp,0.5_dp,-0.2_dp,0.4_dp],[2,2])
  u=reshape([2.0_dp,0.3_dp,0.3_dp,1.5_dp],[2,2])
  v=reshape([1.2_dp,0.2_dp,0.2_dp,0.8_dp],[2,2])
  ld=dmatnorm(a,m,u,v,ok=ok,message=msg)
  call assert_true(ok,'density status: '//trim(msg))
  call assert_close(ld,-4.908844744462413_dp,2e-12_dp,'independent scipy log-density reference')
  d=dmatnorm(a,m,u,v,log_density=.false.)
  call assert_close(log(d),ld,2e-14_dp,'density/log density consistency')
  call assert_close(dmatnorm(reshape([0.5_dp],[1,1]),reshape([0.0_dp],[1,1]), &
    reshape([4.0_dp],[1,1]),reshape([9.0_dp],[1,1])), &
    -0.5_dp*log(2.0_dp*acos(-1.0_dp)*36.0_dp)-0.5_dp*0.25_dp/36.0_dp,1e-13_dp,'univariate normal reduction')
  write(*,'(a)') 'test_density: PASS'
end program test_density
