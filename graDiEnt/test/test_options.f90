program test_options
  use gradient
  implicit none
  type(sqgde_options) :: opt
  integer :: status
  character(len=160) :: message

  call get_algo_params(3,opt)
  if (opt%n_particles /= 9) error stop 'options: default n_particles'
  if (abs(opt%step_size - 2.38_dp/sqrt(6.0_dp)) > 1.0e-14_dp) error stop 'options: step_size'
  call validate_options(opt,status,message)
  if (status /= 0) error stop 'options: valid defaults rejected'
  opt%n_diff = 5
  call validate_options(opt,status,message)
  if (status == 0) error stop 'options: invalid parent count accepted'
  print *, 'PASS test_options'
end program test_options
