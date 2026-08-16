program test_rng
  use mbbefd, only : dp, rmbbefd, mmbbefd, tlmbbefd, roibeta
  implicit none
  integer,parameter :: n=40000
  real(dp),allocatable :: x(:)
  real(dp)::a,b,m,p1
  allocate(x(n));call seed_fixed()
  a=0.5_dp;b=0.3_dp
  call rmbbefd(x,a,b);m=sum(x)/real(n,dp)
  if(abs(m-mmbbefd(1.0_dp,a,b))>0.015_dp) error stop 'mbbefd RNG mean'
  if(abs(real(count(x==1.0_dp),dp)/real(n,dp)-tlmbbefd(a,b))>0.015_dp) error stop 'mbbefd RNG atom'
  p1=0.25_dp;call roibeta(x,2.0_dp,5.0_dp,p1)
  if(abs(real(count(x==1.0_dp),dp)/real(n,dp)-p1)>0.015_dp) error stop 'oibeta RNG atom'
  print '(a)', 'test_rng: PASS'
contains
  subroutine seed_fixed()
    integer,allocatable::s(:);integer::i,nseed
    call random_seed(size=nseed);allocate(s(nseed));do i=1,nseed;s(i)=7919+37*i;end do;call random_seed(put=s)
  end subroutine seed_fixed
end program test_rng
