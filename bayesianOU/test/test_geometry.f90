! SPDX-License-Identifier: MIT
module geometry_test_target
  use bayesianou, only : dp
  implicit none
contains
  function lp_gauss(x) result(v)
    real(dp),intent(in)::x(:)
    real(dp)::v
    v=-0.5_dp*(x(1)**2+x(2)**2/4.0_dp)
  end function lp_gauss
  subroutine grad_gauss(x,g)
    real(dp),intent(in)::x(:)
    real(dp),intent(out)::g(size(x))
    g=[-x(1),-x(2)/4.0_dp]
  end subroutine grad_gauss
  subroutine hess_gauss(x,h)
    real(dp),intent(in)::x(:)
    real(dp),intent(out)::h(size(x),size(x))
    h=0.0_dp;h(1,1)=-1.0_dp;h(2,2)=-0.25_dp
  end subroutine hess_gauss
end module geometry_test_target
program test_geometry
  use bayesianou
  use geometry_test_target
  implicit none
  type(ou_geom_target_type)::target
  type(ou_geom_metric_type)::metric,rm
  type(ou_geom_hmc_result)::out
  procedure(lp_gauss),pointer::lp
  procedure(grad_gauss),pointer::gr
  procedure(hess_gauss),pointer::he
  real(dp)::init(2),mass(2,2),inv(2,2),ld
  integer::status
  lp=>lp_gauss;gr=>grad_gauss;he=>hess_gauss
  target=ou_geom_target(2,lp,gr,he)
  metric=ou_geom_metric_euclidean(2)
  init=0.0_dp
  call ou_geom_hmc(target,metric,0.18_dp,8,1200,300,init,42,out)
  if(out%status/=status_ok)error stop 'HMC failed'
  if(abs(sum(out%draws(:,1))/size(out%draws,1))>0.2_dp)error stop 'bad HMC mean'
  if(out%accept_rate<0.5_dp)error stop 'low HMC acceptance'
  rm=ou_geom_metric_riemannian(target,alpha=10.0_dp)
  call ou_geom_mass(target,rm,[0.0_dp,0.0_dp],mass,inv,ld,status)
  if(status/=status_ok.or.any([(mass(1,1)<=0.0_dp),(mass(2,2)<=0.0_dp)]))error stop 'bad Riemannian mass'
  print *, 'test_geometry: PASS'
end program test_geometry
