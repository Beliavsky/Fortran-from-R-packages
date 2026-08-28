module spacefillr_sobol
use spacefillr_kinds, only: int32,int64,int128,real32,real64
use spacefillr_rng, only: u32_to_i64
use spacefillr_sobol_matrices, only: sobol_matrices
use spacefillr_sobol_directions, only: sobol_directions
implicit none
private
public :: sobol_single, sobol_owen_single, generate_sobol_set, generate_sobol_owen_set
integer(int128), parameter :: two32_128=2_int128**32
contains
pure integer(int32) function low32(x) result(y)
integer(int128),intent(in)::x
integer(int128)::z
z=modulo(x,two32_128)
if(z<=2147483647_int128)then
 y=int(z,int32)
else
 y=int(z-two32_128,int32)
end if
end function
pure integer(int128) function u32i(x) result(y)
integer(int32),intent(in)::x
y=int(u32_to_i64(x),int128)
end function
pure integer(int32) function add32(a,b) result(c)
integer(int32),intent(in)::a,b
c=low32(u32i(a)+u32i(b))
end function
pure integer(int32) function mul32(a,b) result(c)
integer(int32),intent(in)::a,b
c=low32(u32i(a)*u32i(b))
end function
pure integer(int32) function hash_combine(seed,v) result(h)
integer(int32),intent(in)::seed,v
integer(int32)::t
t=add32(v,add32(shiftl(seed,6),shiftr(seed,2)))
h=ieor(seed,t)
end function
pure integer(int32) function hash_u32(n0,seed) result(n)
integer(int32),intent(in)::n0,seed
n=ieor(int(z'6217c6e1',int32),add32(n0,mul32(seed,int(z'9e3779b9',int32))))
n=ieor(n,shiftr(n,17)); n=mul32(n,int(z'ed5ad4bb',int32))
n=ieor(n,shiftr(n,11)); n=mul32(n,int(z'ac4c1b51',int32))
n=ieor(n,shiftr(n,15)); n=mul32(n,int(z'31848bab',int32))
n=ieor(n,shiftr(n,14))
end function
pure integer(int32) function reverse_bits(x0) result(x)
integer(int32),intent(in)::x0
x=x0
x=ior(shiftr(iand(x,int(z'aaaaaaaa',int32)),1),shiftl(iand(x,int(z'55555555',int32)),1))
x=ior(shiftr(iand(x,int(z'cccccccc',int32)),2),shiftl(iand(x,int(z'33333333',int32)),2))
x=ior(shiftr(iand(x,int(z'f0f0f0f0',int32)),4),shiftl(iand(x,int(z'0f0f0f0f',int32)),4))
x=ior(shiftr(iand(x,int(z'ff00ff00',int32)),8),shiftl(iand(x,int(z'00ff00ff',int32)),8))
x=ior(shiftr(x,16),shiftl(x,16))
end function
pure integer(int32) function owen_scramble(x0,seed0) result(x)
integer(int32),intent(in)::x0,seed0
integer(int32)::seed
x=reverse_bits(x0)
seed=hash_u32(seed0,int(z'a14a177d',int32))
x=ieor(x,mul32(x,int(z'3d20adea',int32)))
x=add32(x,seed)
x=mul32(x,ior(shiftr(seed,16),1_int32))
x=ieor(x,mul32(x,int(z'05526c56',int32)))
x=ieor(x,mul32(x,int(z'53a22864',int32)))
x=reverse_bits(x)
end function
pure real(real32) function u32_to_f32(x) result(v)
integer(int32),intent(in)::x
v=min(real(u32_to_i64(x),real32)*2.3283064365386962890625e-10_real32,0.99999994_real32)
end function
pure integer(int32) function sobol_matrix_u32(index,dim,scramble) result(v)
integer(int32),intent(in)::index,scramble
integer,intent(in)::dim
integer(int32)::idx
integer::b
idx=owen_scramble(index,scramble); v=0_int32
if(dim<0 .or. dim>=1024) return
do b=0,31
 if(btest(idx,b)) v=ieor(v,sobol_matrices(b+1,dim+1))
end do
end function
pure integer(int32) function sobol_dir_u32(index,dim) result(v)
integer(int32),intent(in)::index
integer,intent(in)::dim
integer::b
v=0_int32
if(dim<0 .or. dim>=21201) return
do b=0,31
 if(btest(index,b))v=ieor(v,sobol_directions(b+1,dim+1))
end do
end function
pure real(real64) function sobol_single(i,dim,seed) result(x)
integer(int64),intent(in)::i
integer,intent(in)::dim
integer(int32),intent(in),optional::seed
integer(int32)::s,idx
if(dim<0 .or. dim>=1024) error stop 'sobol_single: dimension must be 0..1023'
s=0_int32;if(present(seed))s=seed
idx=low32(int(modulo(i,4294967296_int64),int128))
x=real(u32_to_f32(sobol_matrix_u32(idx,dim,s)),real64)
end function
pure real(real64) function sobol_owen_single(i,dim,seed) result(x)
integer(int64),intent(in)::i
integer,intent(in)::dim
integer(int32),intent(in),optional::seed
integer(int32)::s,idx,a,b
if(dim<0 .or. dim>=21201) error stop 'sobol_owen_single: dimension must be 0..21200'
s=0_int32;if(present(seed))s=seed
idx=low32(int(modulo(i,4294967296_int64),int128))
a=owen_scramble(idx,s)
b=sobol_dir_u32(a,dim)
b=owen_scramble(b,hash_combine(s,int(dim,int32)))
x=real(u32_to_f32(b),real64)
end function
subroutine generate_sobol_set(n,dim,x,seed)
integer,intent(in)::n,dim
real(real64),intent(out)::x(n,dim)
integer(int32),intent(in),optional::seed
integer::i,j
if(dim>1024)error stop 'generate_sobol_set: dim > 1024'
do j=1,dim;do i=1,n;x(i,j)=sobol_single(int(i-1,int64),j-1,seed);end do;end do
end subroutine
subroutine generate_sobol_owen_set(n,dim,x,seed)
integer,intent(in)::n,dim
real(real64),intent(out)::x(n,dim)
integer(int32),intent(in),optional::seed
integer::i,j
if(dim>21201)error stop 'generate_sobol_owen_set: dim > 21201'
do j=1,dim;do i=1,n;x(i,j)=sobol_owen_single(int(i-1,int64),j-1,seed);end do;end do
end subroutine
end module spacefillr_sobol
