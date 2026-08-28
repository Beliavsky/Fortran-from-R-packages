module spacefillr_rng
use spacefillr_kinds, only: int32, int64, int128, real32, real64
implicit none
private
public :: pcg32_state, u32_to_i64
integer(int64), parameter :: pcg_mult = int(z'5851F42D4C957F2D', int64)
integer(int64), parameter :: pcg_inc  = int(z'14057B7EF767814F', int64)
type :: pcg32_state
   integer(int64) :: state = 0_int64
contains
   procedure :: init => pcg_init
   procedure :: next_u32 => pcg_next_u32
   procedure :: uniform32 => pcg_uniform32
   procedure :: uniform => pcg_uniform
   procedure :: uniform_int => pcg_uniform_int
end type
contains
pure integer(int128) function u64_to_i128(x) result(y)
integer(int64), intent(in) :: x
integer(int128), parameter :: two64 = 2_int128**64
if (x >= 0_int64) then
   y = int(x,int128)
else
   y = int(x,int128) + two64
end if
end function
pure integer(int64) function i128_to_u64(x) result(y)
integer(int128), intent(in) :: x
integer(int128), parameter :: two64 = 2_int128**64
integer(int128), parameter :: max64 = 2_int128**63-1_int128
integer(int128) :: z
z = modulo(x,two64)
if (z <= max64) then
   y = int(z,int64)
else
   y = int(z-two64,int64)
end if
end function
pure integer(int64) function add64(a,b) result(c)
integer(int64),intent(in)::a,b
c=i128_to_u64(u64_to_i128(a)+u64_to_i128(b))
end function
pure integer(int64) function mul64(a,b) result(c)
integer(int64),intent(in)::a,b
c=i128_to_u64(u64_to_i128(a)*u64_to_i128(b))
end function
pure integer(int64) function bump64(s) result(t)
integer(int64),intent(in)::s
t=add64(mul64(s,pcg_mult),pcg_inc)
end function
pure integer(int64) function u32_to_i64(x) result(y)
integer(int32),intent(in)::x
if(x>=0_int32) then
 y=int(x,int64)
else
 y=int(x,int64)+4294967296_int64
end if
end function
pure integer(int32) function i64_to_u32(x) result(y)
integer(int64),intent(in)::x
integer(int64)::z
z=modulo(x,4294967296_int64)
if(z<=2147483647_int64) then
 y=int(z,int32)
else
 y=int(z-4294967296_int64,int32)
end if
end function
pure integer(int32) function rotr32(x,r) result(y)
integer(int32),intent(in)::x
integer,intent(in)::r
integer::rr
rr=iand(r,31)
if(rr==0) then
 y=x
else
 y=ior(shiftr(x,rr),shiftl(x,32-rr))
end if
end function
subroutine pcg_init(self,seed)
class(pcg32_state),intent(inout)::self
integer(int64),intent(in)::seed
integer(int64)::s
s=i128_to_u64(int(modulo(seed,4294967296_int64),int128)+u64_to_i128(pcg_inc))
self%state=bump64(s)
end subroutine
integer(int32) function pcg_next_u32(self) result(out)
class(pcg32_state),intent(inout)::self
integer(int64)::old,xs
integer(int32)::x32
integer::rot
old=self%state
self%state=bump64(self%state)
rot=int(iand(shiftr(old,59),31_int64))
xs=ieor(old,shiftr(old,18))
x32=i64_to_u32(shiftr(xs,27))
out=rotr32(x32,rot)
end function
real(real32) function pcg_uniform32(self) result(u)
class(pcg32_state),intent(inout)::self
integer(int32)::x
x=self%next_u32()
u=real(u32_to_i64(x),real32)*2.3283064365386962890625e-10_real32
end function
real(real64) function pcg_uniform(self) result(u)
class(pcg32_state),intent(inout)::self
u=real(self%uniform32(),real64)
end function
integer function pcg_uniform_int(self,lo,hi) result(v)
class(pcg32_state),intent(inout)::self
integer,intent(in)::lo,hi
if(hi<lo) error stop 'pcg_uniform_int: invalid range'
v=int(self%uniform32()*real(hi-lo+1,real32))+lo
if(v>hi)v=hi
end function
end module spacefillr_rng
