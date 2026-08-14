program test_sequences
   use, intrinsic :: iso_fortran_env, only : int32, int64, real64
   use randtoolbox_sfmt, only : sfmt_rng
   use randtoolbox_mt19937, only : mt19937_rng
   use randtoolbox_knuth, only : knuth_rng
   use rngwell, only : well_rng
   implicit none
   type(sfmt_rng) :: sf
   type(mt19937_rng) :: mt
   type(knuth_rng) :: kg
   type(well_rng) :: wg
   integer(int64), parameter :: sfexp(8)=[1196421539_int64,2865311212_int64,3866479472_int64,2692900087_int64, &
      3838928621_int64,3188765817_int64,2751632982_int64,3712143069_int64]
   integer(int64), parameter :: mtexp(10)=[3499211612_int64,581869302_int64,3890346734_int64,3586334585_int64, &
      545404204_int64,4161255391_int64,3922919429_int64,949333985_int64,2715962298_int64,1323567403_int64]
   integer(int32), parameter :: wexp(4)=[int(z'3e099644',int32),int(z'1ed900fc',int32), &
      int(z'95171007',int32),int(z'0941445e',int32)]
   real(real64) :: a(2001), b(2001)
   integer :: i

   call sf%init(607,1234_int64,1)
   do i=1,size(sfexp)
      if(sf%next_uint32()/=sfexp(i)) error stop 'SFMT reference sequence mismatch'
   end do
   call mt%init(5489_int64)
   do i=1,size(mtexp)
      if(mt%next_uint32()/=mtexp(i)) error stop 'MT19937 reference sequence mismatch'
   end do
   call wg%init('512a',seed=123456789_int64)
   do i=1,size(wexp)
      if(wg%next_uint32()/=wexp(i)) error stop 'WELL reference sequence mismatch'
   end do
   call kg%seed(1302_int64); call kg%fill(a)
   call kg%seed(1302_int64); call kg%fill(b)
   if(maxval(abs(a-b))>0.0_real64) error stop 'Knuth reproducibility mismatch'
   print '(a)', 'sequence tests: PASS'
end program test_sequences
