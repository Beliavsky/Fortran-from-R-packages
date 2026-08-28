program test_spacefillr
use, intrinsic :: iso_fortran_env, only: int32, int64, real64
use spacefillr
use spacefillr_rng, only: pcg32_state, u32_to_i64
implicit none

type(pcg32_state) :: r
real(real64) :: x(256,2), x16(16,2), s(8,4), o(8,4), h(8,6)
real(real64) :: hr1(32,6), hr2(32,6), hr3(32,6)
integer(int64), parameter :: pcg_ref(8) = [ &
   1273465047_int64, 4201302492_int64, 1760530922_int64, &
   3811196712_int64, 629196892_int64, 2353405979_int64, &
   1215268356_int64, 1930523812_int64 ]
real(real64), parameter :: pj_ref(4,2) = reshape([ &
   0.2965016961d0, 0.7049527764d0, 0.7739725113d0, &
   0.2247425467d0, 0.9781919718d0, 0.4436816871d0, &
   0.6414758563d0, 0.4866460562d0 ], [4,2])
integer :: i

call r%init(7_int64)
do i = 1, 8
   if (u32_to_i64(r%next_u32()) /= pcg_ref(i)) then
      error stop 'PCG reference mismatch'
   end if
end do

call generate_sobol_set(8,4,s,7_int32)
if (abs(s(1,1)-0.3196877837d0) > 8d-8 .or. &
    abs(s(8,4)-0.9915646911d0) > 8d-8) then
   error stop 'Sobol mismatch'
end if
if (abs(sobol_single(4294967295_int64,1023,123_int32) - &
    0.50342857837677002d0) > 1d-15) then
   error stop 'high-dimensional Sobol mismatch'
end if

call generate_sobol_owen_set(8,4,o,7_int32)
if (abs(o(1,1)-0.5223080516d0) > 8d-8 .or. &
    abs(o(8,4)-0.4919772744d0) > 8d-8) then
   error stop 'Owen mismatch'
end if
if (abs(sobol_owen_single(4294967295_int64,21200,123_int32) - &
    0.05553311109542847d0) > 1d-15) then
   error stop 'high-dimensional Owen mismatch'
end if

call generate_halton_faure_set(8,6,h)
if (maxval(abs(h(:,1) - [0d0,.5d0,.25d0,.75d0,.125d0, &
    .625d0,.375d0,.875d0])) > 1d-7) then
   error stop 'Halton2 mismatch'
end if
if (abs(generate_halton_faure_single(4294967295_int64,255) - &
    0.77262938022613525d0) > 1d-15) then
   error stop 'high-dimensional Halton mismatch'
end if
call generate_halton_random_set(32,6,hr1,7_int64)
call generate_halton_random_set(32,6,hr2,7_int64)
call generate_halton_random_set(32,6,hr3,8_int64)
if (maxval(abs(hr1-hr2)) > 1d-15) then
   error stop 'random Halton reproducibility mismatch'
end if
if (maxval(abs(hr1-hr3)) < 1d-12) then
   error stop 'random Halton seed has no effect'
end if
if (minval(hr1) < 0d0 .or. maxval(hr1) >= 1d0) then
   error stop 'random Halton range mismatch'
end if

call generate_pj_set(16,x16,7_int64)
if (maxval(abs(x16(1:4,:)-pj_ref)) > 2d-7) error stop 'PJ mismatch'

call generate_pmj_set(256,x,123_int64)
if (abs(x(256,2)-0.7667365055531263d0) > 1d-15) then
   error stop 'PMJ mismatch'
end if
call generate_pmjbn_set(256,x,123_int64)
if (abs(x(256,1)-0.3555436320602894d0) > 1d-15 .or. &
    abs(x(256,2)-0.8196288119070232d0) > 1d-15) then
   error stop 'PMJBN mismatch'
end if
call generate_pmj02_set(256,x,123_int64)
if (abs(x(256,2)-0.9589003358269110d0) > 1d-15) then
   error stop 'PMJ02 mismatch'
end if
call generate_pmj02bn_set(256,x,123_int64)
if (abs(x(256,1)-0.32413161429576576d0) > 1d-15 .or. &
    abs(x(256,2)-0.8640170855214819d0) > 1d-15) then
   error stop 'PMJ02BN mismatch'
end if

print *, 'spacefillr tests passed'
end program test_spacefillr
