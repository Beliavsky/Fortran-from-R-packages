program test_probability
  use gpareto, only : dp, probability_nondomination
  implicit none
  real(dp)::front(2,2),mean(3,2),sd(3,2),f3(1,3),m3(1,3),s3(1,3),expect
  real(dp),allocatable::p(:)
  front=reshape([0.2_dp,0.8_dp, 0.8_dp,0.2_dp],[2,2])
  mean=reshape([0.1_dp,0.5_dp,0.9_dp, 0.1_dp,0.5_dp,0.9_dp],[3,2]);sd=0.2_dp
  call probability_nondomination(front,mean,sd,p)
  if(any(p<0.0_dp).or.any(p>1.0_dp)) error stop 'probability bounds failed'
  if(.not.(p(1)>p(2).and.p(2)>p(3))) error stop 'probability ordering failed'
  f3(1,:)=[0.2_dp,0.3_dp,0.4_dp];m3(1,:)=[0.1_dp,0.2_dp,0.3_dp];s3=0.5_dp
  call probability_nondomination(f3,m3,s3,p)
  expect=1.0_dp-(1.0_dp-0.5_dp*erfc(-(0.2_dp-0.1_dp)/0.5_dp/sqrt(2.0_dp)))* &
    (1.0_dp-0.5_dp*erfc(-(0.3_dp-0.2_dp)/0.5_dp/sqrt(2.0_dp)))* &
    (1.0_dp-0.5_dp*erfc(-(0.4_dp-0.3_dp)/0.5_dp/sqrt(2.0_dp)))
  if(abs(p(1)-expect)>1.0e-12_dp) error stop '3D probability failed'
  print *, 'test_probability PASS'
end program test_probability
