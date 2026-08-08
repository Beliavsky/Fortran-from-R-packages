! SPDX-License-Identifier: GPL-3.0-only
program demo_nmof
   use nmof
   implicit none
   type(option_result) :: option
   type(optimization_result) :: opt
   real(dp) :: covariance(3,3), weights(3)
   integer :: status

   option=vanilla_option_european(100.0_dp,100.0_dp,1.0_dp,0.02_dp,0.01_dp,0.20_dp**2,'call')
   write(*,'(a,f12.6)') 'European call value: ',option%value
   write(*,'(a,f12.6)') 'Delta:               ',option%delta

   covariance=reshape([0.04_dp,0.006_dp,0.004_dp, &
                       0.006_dp,0.09_dp,0.012_dp, &
                       0.004_dp,0.012_dp,0.16_dp],[3,3])
   call minimum_variance(covariance,weights,status=status)
   write(*,'(a,3f11.6)') 'Minimum-variance weights: ',weights

   call de_opt(sphere,[-5.0_dp,-5.0_dp],[5.0_dp,5.0_dp],opt, &
      n_population=30,n_generations=100,minmax_constraint=.true.,seed=20260728_i8)
   write(*,'(a,es13.5)') 'DE sphere objective: ',opt%ofvalue
contains
   function sphere(x,context) result(value)
      real(dp),intent(in)::x(:)
      class(*),intent(in),optional::context
      real(dp)::value
      value=dot_product(x,x)
   end function sphere
end program demo_nmof
