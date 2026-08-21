program test_energy_gof
  use compositions
  implicit none
  real(dp) :: x(8,3),same(8,3),z(20,3),s0,s1
  integer :: sizes(2),i
  type(energy_test_result) :: e,norm
  x=reshape([ &
    0.20_dp,0.22_dp,0.18_dp,0.21_dp,0.55_dp,0.52_dp,0.58_dp,0.54_dp, &
    0.30_dp,0.31_dp,0.29_dp,0.32_dp,0.25_dp,0.27_dp,0.23_dp,0.26_dp, &
    0.50_dp,0.47_dp,0.53_dp,0.47_dp,0.20_dp,0.21_dp,0.19_dp,0.20_dp],[8,3])
  sizes=[4,4]
  same(1:4,:)=x(1:4,:); same(5:8,:)=x(1:4,:)
  s0=energy_ksample_statistic(same,sizes)
  s1=energy_ksample_statistic(x,sizes)
  if(abs(s0)>1.0e-12_dp) error stop 'identical samples energy statistic'
  if(s1<=s0) error stop 'separated samples not detected'
  e=acomp_energy_test(x,sizes,reps=39,seed=1234)
  if(e%p_value<0.0_dp.or.e%p_value>1.0_dp) error stop 'energy p-value range'
  do i=1,20
    z(i,1)=0.2_dp+0.015_dp*real(mod(3*i,7),dp)
    z(i,2)=0.3_dp+0.010_dp*real(mod(5*i,9),dp)
    z(i,3)=1.0_dp-z(i,1)-z(i,2)
  end do
  norm=acomp_normal_energy_test(z,reps=19,seed=44)
  if(norm%statistic/=norm%statistic.or.norm%statistic<0.0_dp) error stop 'normal energy statistic'
  if(norm%p_value<0.0_dp.or.norm%p_value>1.0_dp) error stop 'normal energy p-value'
  print *, 'test_energy_gof: PASS'
end program test_energy_gof
