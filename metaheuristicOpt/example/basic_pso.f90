program basic_pso
   use metaheuristic_opt, only : dp, mh_control, mh_result, metaopt, sphere
   implicit none
   type(mh_control) :: control
   type(mh_result) :: result
   real(dp) :: lower(5), upper(5)

   lower = -10.0_dp
   upper = 10.0_dp
   control%num_population = 40
   control%max_iter = 300
   control%seed = 12345
   control%legacy_quirks = .false.

   call metaopt('PSO', sphere, lower, upper, result, control)
   print '(a,es14.6)', 'best objective = ', result%value
   print '(a,*(f11.6,1x))', 'best point     = ', result%par
end program basic_pso
