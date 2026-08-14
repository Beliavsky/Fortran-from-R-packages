program basic
   use, intrinsic :: iso_fortran_env, only : int64, real64
   use randtoolbox, only : sfmt_generate, halton, frequency_test, rng_test_result
   implicit none
   real(real64),allocatable::u(:,:),q(:,:)
   type(rng_test_result)::test
   u=sfmt_generate(1000,1,mexp=19937,seed=12345_int64,use_parameter_sets=.false.)
   q=halton(5,2,start=1_int64)
   test=frequency_test(u(:,1),10)
   print '(a,f10.6)', 'SFMT mean: ',sum(u(:,1))/real(size(u,1),real64)
   print '(a,f10.6)', 'frequency-test p-value: ',test%p_value
   print '(a,2f10.6)', 'first Halton point: ',q(1,:)
end program basic
