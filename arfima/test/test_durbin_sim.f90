program test_durbin_sim
  use arfima
  use test_support
  implicit none
  real(dp),allocatable::r(:),z(:),inv(:,:),cov(:,:),x(:)
  real(dp)::ll,s,ldet,expected
  integer::i,j,info
  type(arfima_error)::err

  call tacvf_arma([0.4_dp],[real(dp)::],4,1.0_dp,r,err)
  x=[0.2_dp,-0.1_dp,0.5_dp,0.3_dp]
  ll=dl_loglikelihood(r,x,err)
  allocate(cov(4,4))
  do i=1,4; do j=1,4; cov(i,j)=r(abs(i-j)+1); end do; end do
  call invert_matrix(cov,inv,info)
  s=dot_product(x,matmul(inv,x))
  call logdet_positive_definite(cov,ldet,info)
  expected=-2.0_dp*log(s/4.0_dp)-0.5_dp*ldet
  call assert_close(ll,expected,2.0e-11_dp,'Durbin likelihood equals direct covariance likelihood')

  call dl_simulate(r,[1.0_dp,2.0_dp,-1.0_dp,0.5_dp],z,err)
  call assert_close(z(1),sqrt(1.0_dp/(1.0_dp-0.16_dp)),1.0e-12_dp,'AR simulation first value')
  call assert_close(z(2),0.4_dp*z(1)+2.0_dp,1.0e-12_dp,'AR simulation second value')
  call assert_close(z(3),0.4_dp*z(2)-1.0_dp,1.0e-12_dp,'AR simulation third value')
  write(*,'(a)') 'test_durbin_sim: PASS'
end program test_durbin_sim
