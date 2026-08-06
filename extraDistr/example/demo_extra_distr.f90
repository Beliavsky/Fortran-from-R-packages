program demo_extra_distr
  use extra_distr
  implicit none
  integer, parameter :: n = 20000
  real(dp), allocatable :: z(:,:), tail(:)
  real(dp) :: mean1, mean2, cov12, q99
  integer, allocatable :: category(:)

  call seed_rng(20260804)

  q99 = qgpd(0.99_dp, mu=0.0_dp, sigma=1.5_dp, xi=0.2_dp)
  tail = rgpd(n, mu=0.0_dp, sigma=1.5_dp, xi=0.2_dp)

  z = rbvnorm(n, mean1=0.0_dp, mean2=0.0_dp, sd1=1.0_dp, sd2=2.0_dp, cor=0.6_dp)
  mean1 = sum(z(:,1))/real(n,dp)
  mean2 = sum(z(:,2))/real(n,dp)
  cov12 = sum((z(:,1)-mean1)*(z(:,2)-mean2))/real(n-1,dp)

  category = rcat(12,[0.1_dp,0.2_dp,0.7_dp])

  write(*,'(a,f12.6)') 'GPD 99% quantile:           ', q99
  write(*,'(a,f12.6)') 'Simulated GPD mean:         ', sum(tail)/real(n,dp)
  write(*,'(a,2f12.6)') 'Bivariate-normal means:     ', mean1, mean2
  write(*,'(a,f12.6)') 'Bivariate-normal covariance:', cov12
  write(*,'(a,12(1x,i0))') 'Categorical draws:          ', category
end program demo_extra_distr
