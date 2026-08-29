program test_random
  use matrixNormal
  use test_support
  implicit none
  integer,parameter::s=6000
  real(dp)::m(2,2),u(2,2),v(2,2)
  real(dp),allocatable::x(:,:),xs(:,:,:),z(:,:)
  real(dp)::mu(4),cov(4,4),d(4)
  integer::i,j,k
  m=reshape([0.5_dp,-0.25_dp,1.0_dp,2.0_dp],[2,2])
  u=reshape([1.0_dp,0.2_dp,0.2_dp,2.0_dp],[2,2])
  v=reshape([3.0_dp,0.4_dp,0.4_dp,1.5_dp],[2,2])
  x=rmatnorm(m,u,v,777)
  xs=rmatnorm(1,m,u,v,777)
  call assert_true(maxval(abs(x-xs(:,:,1)))<1e-14_dp,'one/many seeded reproducibility')
  xs=rmatnorm(s,m,u,v,2468)
  allocate(z(s,4))
  do k=1,s
    z(k,:)=vec(xs(:,:,k))
  end do
  mu=sum(z,dim=1)/real(s,dp)
  cov=0.0_dp
  do k=1,s
    do i=1,4
      do j=1,4
        cov(i,j)=cov(i,j)+(z(k,i)-mu(i))*(z(k,j)-mu(j))
      end do
    end do
  end do
  cov=cov/real(s-1,dp)
  d=vec(m)
  call assert_true(maxval(abs(mu-d))<0.12_dp,'random sample mean')
  d=[3.0_dp,6.0_dp,1.5_dp,3.0_dp]
  call assert_true(maxval(abs([(cov(i,i),i=1,4)]-d)/d)<0.10_dp,'random sample variances')
  call assert_close(cov(1,2),0.6_dp,0.18_dp,'row covariance')
  call assert_close(cov(1,3),0.4_dp,0.16_dp,'column covariance')
  write(*,'(a)') 'test_random: PASS'
end program test_random
