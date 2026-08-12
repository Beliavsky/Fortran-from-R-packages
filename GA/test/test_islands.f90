program test_islands
  use ga
  implicit none
  type(island_control_type)::c
  type(island_real_result)::rr
  type(island_int_result)::rb,rp
  real(dp)::lo(2),up(2)

  lo=-5.0_dp;up=5.0_dp
  c%pop_size=120;c%num_islands=4;c%migration_rate=0.1_dp
  c%migration_interval=20;c%max_iter=200;c%seed=7654
  c%ga%pcrossover=0.9_dp;c%ga%pmutation=0.2_dp;c%ga%elitism=2
  call ga_islands_real(sphere_fit,lo,up,c,rr)
  if(rr%fitness_value < -1.0e-5_dp) error stop "real island GA fitness"
  if(maxval(abs(rr%solution))>0.01_dp) error stop "real island GA solution"

  c%pop_size=80;c%migration_interval=10;c%max_iter=100;c%seed=7655
  call ga_islands_binary(onemax,16,c,rb)
  if(abs(rb%fitness_value-16.0_dp)>1.0e-12_dp) error stop "binary island GA"

  c%pop_size=100;c%migration_interval=15;c%max_iter=150;c%seed=7656
  c%ga%pmutation=0.3_dp
  call ga_islands_permutation(perm_fit,1,6,c,rp)
  if(abs(rp%fitness_value-6.0_dp)>1.0e-12_dp) error stop "permutation island GA"

  print *, "test_islands: PASS",rr%fitness_value,rb%fitness_value,rp%fitness_value
contains
  function sphere_fit(x) result(f)
    real(dp),intent(in)::x(:)
    real(dp)::f
    f=-sum(x*x)
  end function sphere_fit

  function onemax(x) result(f)
    integer,intent(in)::x(:)
    real(dp)::f
    f=real(sum(x),dp)
  end function onemax

  function perm_fit(x) result(f)
    integer,intent(in)::x(:)
    integer,parameter::target(6)=[3,1,6,2,5,4]
    real(dp)::f
    f=real(count(x==target),dp)
  end function perm_fit
end program test_islands
