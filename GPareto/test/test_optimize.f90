program test_optimize
  use gpareto, only : dp, differential_evolution_max, optim_control
  implicit none
  real(dp),allocatable::x(:);real(dp)::v
  type(optim_control)::ctl
  ctl%population=24;ctl%generations=80;ctl%seed=77
  call differential_evolution_max(obj,[-2.0_dp,-2.0_dp],[2.0_dp,2.0_dp],x,v,ctl)
  if(sum((x-[0.3_dp,-0.4_dp])**2)>1.0e-5_dp) then
    print *,x,v; error stop 'DE optimizer failed'
  end if
  print *, 'test_optimize PASS'
contains
  function obj(z) result(f)
    real(dp),intent(in)::z(:);real(dp)::f
    f=-sum((z-[0.3_dp,-0.4_dp])**2)
  end function obj
end program test_optimize
