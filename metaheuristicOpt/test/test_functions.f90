program test_functions
   use metaheuristic_opt, only : dp, sphere, schwefel, rastrigin, &
      cantilever_beam, centilever_beam, mae
   implicit none
   real(dp) :: x(5), data(2,4), coef(4)

   x = 0.0_dp
   if (abs(sphere(x)) > tiny(1.0_dp)) error stop 'sphere zero failed'
   if (abs(schwefel(x)) > tiny(1.0_dp)) error stop 'schwefel zero failed'
   if (abs(rastrigin(x)) > tiny(1.0_dp)) error stop 'rastrigin zero failed'
   x = 10.0_dp
   if (cantilever_beam(x) <= 0.0_dp) error stop 'cantilever benchmark failed'
   if (abs(centilever_beam(x)-cantilever_beam(x)) > tiny(1.0_dp)) error stop 'legacy beam alias failed'

   data(1,:) = [1.0_dp, 2.0_dp, 3.0_dp, 7.0_dp]
   data(2,:) = [2.0_dp, 1.0_dp, 0.0_dp, 4.0_dp]
   coef = [1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp]
   if (abs(mae(data,coef)) > tiny(1.0_dp)) error stop 'MAE helper failed'
end program test_functions
