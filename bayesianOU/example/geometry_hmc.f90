! SPDX-License-Identifier: MIT
module geometry_example_target
  use bayesianou, only : dp
  implicit none
contains
  function logp(x) result(v)
    real(dp),intent(in)::x(:)
    real(dp)::v
    v=-0.5_dp*(x(1)**2+x(2)**2/9.0_dp)
  end function logp
  subroutine grad(x,g)
    real(dp),intent(in)::x(:)
    real(dp),intent(out)::g(size(x))
    g=[-x(1),-x(2)/9.0_dp]
  end subroutine grad
  subroutine hess(x,h)
    real(dp),intent(in)::x(:)
    real(dp),intent(out)::h(size(x),size(x))
    h=0.0_dp;h(1,1)=-1.0_dp;h(2,2)=-1.0_dp/9.0_dp
  end subroutine hess
end module geometry_example_target
program geometry_hmc
  use bayesianou
  use geometry_example_target
  implicit none
  type(ou_geom_target_type)::target
  type(ou_geom_metric_type)::metric
  type(ou_geom_hmc_result)::sample
  procedure(logp),pointer::lp_ptr
  procedure(grad),pointer::g_ptr
  procedure(hess),pointer::h_ptr
  lp_ptr=>logp;g_ptr=>grad;h_ptr=>hess
  target=ou_geom_target(2,lp_ptr,g_ptr,h_ptr)
  metric=ou_geom_metric_riemannian(target,alpha=20.0_dp)
  call ou_geom_hmc(target,metric,0.12_dp,6,800,200,[0.0_dp,0.0_dp],123,sample)
  print '(a,f8.3)','acceptance: ',sample%accept_rate
  print '(a,2f10.4)','posterior mean: ',sum(sample%draws,dim=1)/real(size(sample%draws,1),dp)
  print '(a,f8.3)','E-BFMI: ',sample%ebfmi
end program geometry_hmc
