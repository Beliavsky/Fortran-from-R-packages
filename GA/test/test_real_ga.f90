program test_real_ga
  use ga
  implicit none
  type(ga_control_type) :: c
  type(ga_real_result) :: r
  real(dp) :: lo(3),up(3),suggestion(1,3)
  lo=-5.0_dp
  up=5.0_dp
  c%pop_size=80
  c%max_iter=300
  c%run=100
  c%seed=2026
  c%pcrossover=0.9_dp;c%pmutation=0.15_dp;c%elitism=4
  call ga_real(rosenbrock_fit,lo,up,c,r)
  if(r%fitness_value < -1.0e-3_dp) error stop "real GA fitness"
  if(maxval(abs(r%solution-[1.0_dp,-2.0_dp,0.5_dp]))>0.05_dp) error stop "real GA solution"
  if(r%evaluations<=c%pop_size) error stop "real GA evaluations"
  suggestion(1,:)=[1.0_dp,-2.0_dp,0.5_dp]
  c%max_iter=50
  c%max_fitness=0.0_dp
  c%seed=2027
  call ga_real(rosenbrock_fit,lo,up,c,r,suggestion)
  if(r%iter/=1) error stop "maxFitness stopping"
  if(abs(r%fitness_value)>1.0e-14_dp) error stop "suggestion handling"
  print *, "test_real_ga: PASS",r%fitness_value,r%solution
contains
  function rosenbrock_fit(x) result(f)
    real(dp),intent(in)::x(:)
    real(dp)::f
    f=-sum((x-[1.0_dp,-2.0_dp,0.5_dp])**2)
  end function rosenbrock_fit
end program test_real_ga
