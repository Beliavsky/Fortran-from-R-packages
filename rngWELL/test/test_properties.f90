program test_properties
   use, intrinsic :: iso_fortran_env, only : int32, int64, real64
   use rngwell, only : well_rng, well_state_size, well_variant_supported, init_mt2002
   implicit none
   character(len=6), parameter :: variants(17) = [character(len=6) :: &
      '512a','521a','521b','607a','607b','800a','800b','1024a','1024b', &
      '19937a','19937b','19937c','21701a','23209a','23209b','44497a','44497b']
   type(well_rng) :: r
   integer(int32), allocatable :: s(:)
   real(real64) :: x(10000), mean, var
   integer :: i

   do i=1,size(variants)
      if (.not. well_variant_supported(variants(i))) error stop 'supported variant rejected'
      if (well_state_size(variants(i)) <= 0) error stop 'bad state size'
      call r%init(variants(i), seed=987654321_int64)
      call r%fill(x)
      if (any(x < 0.0_real64) .or. any(x >= 1.0_real64)) error stop 'uniform outside [0,1)'
      mean=sum(x)/real(size(x),real64)
      var=sum((x-mean)**2)/real(size(x)-1,real64)
      ! Loose smoke bounds only; exact sequence tests provide correctness.
      if (abs(mean-0.5_real64) > 0.03_real64) error stop 'uniform mean smoke test failed'
      if (abs(var-1.0_real64/12.0_real64) > 0.01_real64) error stop 'uniform variance smoke test failed'
   end do

   allocate(s(16))
   call init_mt2002(1_int64,s)
   if (s(1) /= 1_int32) error stop 'MT2002 first state failed'
   if (s(2) /= int(z'6c078966',int32)) error stop 'MT2002 second state failed'
   deallocate(s)

   print '(a)', 'rngWELL property tests: PASS'
end program test_properties
