program basic
use, intrinsic :: iso_fortran_env, only: real64,int32
use spacefillr
implicit none
real(real64)::x(16,2)
integer::i
call generate_sobol_owen_set(16,2,x,7_int32)
do i=1,5
 write(*,'(i3,2f14.9)') i,x(i,1),x(i,2)
end do
end program basic
