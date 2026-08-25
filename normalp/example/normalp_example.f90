program normalp_example
  use normalp
  implicit none
  real(dp) :: x(5)
  call rnormp(x, mu=0.0_dp, sigmap=1.0_dp, p=2.0_dp)
  print '(a,f10.6)', 'density at zero = ', dnormp(0.0_dp)
  print '(a,5f10.5)', 'sample = ', x
end program
