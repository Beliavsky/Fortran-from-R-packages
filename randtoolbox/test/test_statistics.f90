program test_statistics
   use, intrinsic :: iso_fortran_env, only : real64
   use randtoolbox_tests
   implicit none
   real(real64),allocatable::s(:)
   integer,allocatable::p(:,:)
   type(rng_test_result)::r
   real(real64)::u(16)
   integer::i
   s=stirling_second(4)
   if(any(abs(s-[0.0_real64,1.0_real64,7.0_real64,6.0_real64,1.0_real64])>1e-14_real64)) error stop 'Stirling mismatch'
   p=permutations(3)
   if(size(p,1)/=6 .or. size(p,2)/=3) error stop 'permutation shape mismatch'
   if(collision_count([0,1,0,2,1],4)/=2) error stop 'collision count mismatch'
   do i=1,16; u(i)=(real(i,real64)-0.5_real64)/16.0_real64; end do
   r=frequency_test(u,16)
   if(abs(r%statistic)>1e-14_real64) error stop 'frequency test mismatch'
   if(any(abs(r%observed-1.0_real64)>1e-14_real64)) error stop 'frequency counts mismatch'
   print '(a)', 'statistical tests: PASS'
end program test_statistics
