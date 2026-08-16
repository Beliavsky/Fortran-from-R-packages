program test_fit
  use mbbefd, only : dp, dr_fit_result, fit_dr, roiunif, roibeta, rmbbefd, rmbbefd_gb
  implicit none
  type(dr_fit_result)::fit
  real(dp),allocatable::x(:)
  integer,parameter::n=1200
  allocate(x(n));call seed_fixed()
  call roiunif(x,0.27_dp);call fit_dr(x,'oiunif',fit)
  if(fit%convergence/=0.or.abs(fit%estimate(1)-0.27_dp)>0.05_dp) error stop 'oiunif fit'
  call roibeta(x,2.0_dp,4.0_dp,0.18_dp);call fit_dr(x,'oibeta',fit)
  if(fit%convergence/=0) error stop 'oibeta fit convergence'
  if(abs(fit%estimate(3)-0.18_dp)>0.05_dp) error stop 'oibeta p1'
  call rmbbefd(x,0.6_dp,0.25_dp);call fit_dr(x,'mbbefd',fit,method='tlmme')
  if(fit%convergence/=0) error stop 'mbbefd tlmme'
  call fit_dr(x,'mbbefd',fit,method='mle')
  if(fit%convergence/=0.or.fit%estimate(1)<=0.0_dp.or.fit%estimate(2)>=1.0_dp) error stop 'mbbefd mle'
  call rmbbefd_gb(x,2.5_dp,0.2_dp);call fit_dr(x,'MBBEFD',fit,method='tlmme')
  if(fit%convergence/=0) error stop 'MBBEFD tlmme'
  call fit_dr(x,'MBBEFD',fit,method='mle')
  if(fit%convergence/=0.or.fit%estimate(1)<=1.0_dp) error stop 'MBBEFD mle'
  print '(a)', 'test_fit: PASS'
contains
  subroutine seed_fixed()
    integer,allocatable::s(:);integer::i,nseed
    call random_seed(size=nseed);allocate(s(nseed));do i=1,nseed;s(i)=12347+53*i;end do;call random_seed(put=s)
  end subroutine seed_fixed
end program test_fit
