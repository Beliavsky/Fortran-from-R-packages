! SPDX-License-Identifier: GPL-3.0-only
program optimization_example
   use nmof
   implicit none
   type(optimization_result) :: result

   call de_opt(trefethen_objective,[-1.0_dp,-1.0_dp],[1.0_dp,1.0_dp],result, &
      n_population=50,n_generations=250,minmax_constraint=.true.,seed=987654_i8)
   write(*,'(a,2f14.8)') 'xbest: ',result%xbest
   write(*,'(a,f14.8)') 'objective: ',result%ofvalue
contains
   function trefethen_objective(x,context) result(value)
      real(dp),intent(in)::x(:)
      class(*),intent(in),optional::context
      real(dp)::value
      value=exp(sin(50.0_dp*x(1)))+sin(60.0_dp*exp(x(2)))+ &
         sin(70.0_dp*sin(x(1)))+sin(sin(80.0_dp*x(2)))- &
         sin(10.0_dp*(x(1)+x(2)))+(x(1)**2+x(2)**2)/4.0_dp
   end function trefethen_objective
end program optimization_example
