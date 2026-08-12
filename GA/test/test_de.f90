program test_de
  use ga
  implicit none
  type(de_control_type)::c
  type(ga_real_result)::r
  real(dp)::lo(4),up(4)
  lo=-5.0_dp;up=5.0_dp;c%pop_size=50;c%max_iter=300;c%run=80;c%seed=333
  c%stepsize=0.8_dp;c%pcrossover=0.7_dp
  call de_real(sphere_fit,lo,up,c,r)
  if(r%fitness_value < -1.0e-10_dp) error stop "DE fitness"
  if(maxval(abs(r%solution))>1.0e-4_dp) error stop "DE solution"
  print *, "test_de: PASS",r%fitness_value
contains
  function sphere_fit(x) result(f)
    real(dp),intent(in)::x(:);real(dp)::f
    f=-sum(x*x)
  end function sphere_fit
end program test_de
