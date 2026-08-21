program test_distributions
  use compositions
  implicit none
  real(dp) :: alpha(3),x0(3),m(3),covclr(3,3),den
  real(dp), allocatable :: x(:,:),lnx(:,:)
  integer :: i
  alpha=[2.0_dp,3.0_dp,5.0_dp]; x0=[0.2_dp,0.3_dp,0.5_dp]
  call rng_seed(12345)
  den=dirichlet_pdf(x0,[1.0_dp,1.0_dp,1.0_dp])
  if(abs(den-2.0_dp)>1e-12_dp) error stop 'Dirichlet uniform density'
  x=rdirichlet(10000,alpha); m=sum(x,dim=1)/real(size(x,1),dp)
  if(maxval(abs(m-alpha/sum(alpha)))>0.015_dp) then; print *,m; error stop 'Dirichlet RNG mean'; end if
  covclr=0.0_dp
  covclr(1,1)=0.2_dp; covclr(2,2)=0.2_dp; covclr(3,3)=0.2_dp
  covclr=covclr-0.2_dp/3.0_dp
  lnx=rlogistic_normal(3000,[0.25_dp,0.35_dp,0.40_dp],covclr)
  if(maxval(abs(sum(lnx,dim=2)-1.0_dp))>1e-12_dp) error stop 'logistic normal closure'
  if(.not.(logistic_normal_pdf(x0,[0.25_dp,0.35_dp,0.40_dp],covclr)>0.0_dp)) error stop 'lognormal density'
  print *, 'test_distributions: PASS'
end program
