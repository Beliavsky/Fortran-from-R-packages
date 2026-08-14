module spacefillr_halton
use spacefillr_kinds, only: int32,int64,real32,real64
use spacefillr_rng, only: pcg32_state
implicit none
private
public :: halton_sampler, generate_halton_faure_set, generate_halton_random_set, &
          generate_halton_faure_single, generate_halton_random_single
integer, parameter :: primes(256)=[ &
  2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, &
  59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 131, &
  137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193, 197, 199, 211, 223, &
  227, 229, 233, 239, 241, 251, 257, 263, 269, 271, 277, 281, 283, 293, 307, 311, &
  313, 317, 331, 337, 347, 349, 353, 359, 367, 373, 379, 383, 389, 397, 401, 409, &
  419, 421, 431, 433, 439, 443, 449, 457, 461, 463, 467, 479, 487, 491, 499, 503, &
  509, 521, 523, 541, 547, 557, 563, 569, 571, 577, 587, 593, 599, 601, 607, 613, &
  617, 619, 631, 641, 643, 647, 653, 659, 661, 673, 677, 683, 691, 701, 709, 719, &
  727, 733, 739, 743, 751, 757, 761, 769, 773, 787, 797, 809, 811, 821, 823, 827, &
  829, 839, 853, 857, 859, 863, 877, 881, 883, 887, 907, 911, 919, 929, 937, 941, &
  947, 953, 967, 971, 977, 983, 991, 997, 1009, 1013, 1019, 1021, 1031, 1033, 1039, 1049, &
  1051, 1061, 1063, 1069, 1087, 1091, 1093, 1097, 1103, 1109, 1117, 1123, 1129, 1151, 1153, 1163, &
  1171, 1181, 1187, 1193, 1201, 1213, 1217, 1223, 1229, 1231, 1237, 1249, 1259, 1277, 1279, 1283, &
  1289, 1291, 1297, 1301, 1303, 1307, 1319, 1321, 1327, 1361, 1367, 1373, 1381, 1399, 1409, 1423, &
  1427, 1429, 1433, 1439, 1447, 1451, 1453, 1459, 1471, 1481, 1483, 1487, 1489, 1493, 1499, 1511, &
  1523, 1531, 1543, 1549, 1553, 1559, 1567, 1571, 1579, 1583, 1597, 1601, 1607, 1609, 1613, 1619 ]
integer, parameter :: halton_digits(256)=[ &
  0, 20, 12, 9, 8, 8, 6, 6, 7, 6, 6, 6, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 4, &
  4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, &
  4, 4, 4, 4, 4, 4, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, &
  3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, &
  3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, &
  3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, &
  3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, &
  3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, &
  3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, &
  3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, &
  3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3 ]
type :: halton_sampler
   logical :: random_mode=.false.
   integer(int64) :: seed=0_int64
   integer, allocatable :: offsets(:), perms(:)
contains
   procedure :: init_faure => halton_init_faure
   procedure :: init_random => halton_init_random
   procedure :: sample => halton_sample
end type
contains
recursive subroutine faure_perm(base,p)
integer,intent(in)::base
integer,allocatable,intent(out)::p(:)
integer,allocatable::q(:)
integer::b,i
allocate(p(0:base-1))
if(base<=3)then
 do i=0,base-1;p(i)=i;end do;return
end if
b=base/2
if(iand(base,1)==1)then
 call faure_perm(base-1,q)
 do i=0,base-2;p(i+merge(1,0,i>=b))=q(i)+merge(1,0,q(i)>=b);end do
 p(b)=b
else
 call faure_perm(b,q)
 do i=0,b-1;p(i)=2*q(i);p(b+i)=2*q(i)+1;end do
end if
end subroutine
subroutine halton_init_faure(self)
class(halton_sampler),intent(inout)::self
if(allocated(self%perms))deallocate(self%perms)
if(allocated(self%offsets))deallocate(self%offsets)
self%random_mode=.false.;self%seed=0_int64
end subroutine
subroutine shuffle_perm(p,rng)
integer,intent(inout)::p(0:)
type(pcg32_state),intent(inout)::rng
integer::i,j,t
do i=ubound(p,1),1,-1
 j=rng%uniform_int(0,i)
 t=p(i);p(i)=p(j);p(j)=t
end do
end subroutine
subroutine halton_init_random(self,seed)
class(halton_sampler),intent(inout)::self
integer(int64),intent(in)::seed
type(pcg32_state)::rng
integer,allocatable::p(:)
integer::base,i,total,off
if(allocated(self%perms))deallocate(self%perms)
if(allocated(self%offsets))deallocate(self%offsets)
allocate(self%offsets(0:1619));self%offsets=-1
total=sum([(base,base=4,1619)])
allocate(self%perms(0:total-1));off=0
call rng%init(seed)
do base=4,1619
 self%offsets(base)=off
 allocate(p(0:base-1));do i=0,base-1;p(i)=i;end do
 call shuffle_perm(p,rng)
 self%perms(off:off+base-1)=p;off=off+base;deallocate(p)
end do
self%random_mode=.true.;self%seed=seed
end subroutine
subroutine get_perm(self,base,p)
class(halton_sampler),intent(in)::self
integer,intent(in)::base
integer,allocatable,intent(out)::p(:)
integer::off
if(.not.self%random_mode .or. base<=3)then
 call faure_perm(base,p)
else
 allocate(p(0:base-1));off=self%offsets(base);p=self%perms(off:off+base-1)
end if
end subroutine
pure integer(int32) function reverse32(x0) result(x)
integer(int32),intent(in)::x0
x=x0
x=ior(shiftr(iand(x,int(z'aaaaaaaa',int32)),1),shiftl(iand(x,int(z'55555555',int32)),1))
x=ior(shiftr(iand(x,int(z'cccccccc',int32)),2),shiftl(iand(x,int(z'33333333',int32)),2))
x=ior(shiftr(iand(x,int(z'f0f0f0f0',int32)),4),shiftl(iand(x,int(z'0f0f0f0f',int32)),4))
x=ior(shiftr(iand(x,int(z'ff00ff00',int32)),8),shiftl(iand(x,int(z'00ff00ff',int32)),8))
x=ior(shiftr(x,16),shiftl(x,16))
end function
real(real64) function halton_sample(self,dimension,index) result(v)
class(halton_sampler),intent(in)::self
integer,intent(in)::dimension
integer(int64),intent(in)::index
integer,allocatable::p(:)
integer::base,digit,k
integer(int64)::n
integer(int32)::ix
real(real32)::vf
if(dimension<0 .or. dimension>=256)error stop 'halton_sample: dimension must be 0..255'
base=primes(dimension+1)
if(base==2)then
 ix=int(modulo(index,4294967296_int64),int32)
 ix=reverse32(ix)
 ! Match the generated C++ sampler's 23-bit float mantissa construction.
 vf=real(shiftr(ix,9),real32)/8388608.0_real32
 v=real(vf,real64);return
end if
call get_perm(self,base,p)
n=modulo(index,4294967296_int64)
block
 integer(int64) :: numer,denom
 integer :: nd
 numer=0_int64
 denom=1_int64
 nd=halton_digits(dimension+1)
 do k=1,nd
  digit=int(modulo(n,int(base,int64)))
  numer=numer*int(base,int64)+int(p(digit),int64)
  n=n/int(base,int64)
  denom=denom*int(base,int64)
 end do
 vf=real(numer,real32)*real((1.0_real64-2.0_real64**(-23))/real(denom,real64),real32)
end block
v=real(vf,real64)
end function
real(real64) function generate_halton_faure_single(i,dim) result(v)
integer(int64),intent(in)::i
integer,intent(in)::dim
type(halton_sampler)::h
call h%init_faure();v=h%sample(dim,i)
end function
real(real64) function generate_halton_random_single(i,dim,seed) result(v)
integer(int64),intent(in)::i,seed
integer,intent(in)::dim
type(halton_sampler)::h
call h%init_random(seed);v=h%sample(dim,i)
end function
subroutine generate_halton_faure_set(n,dim,x)
integer,intent(in)::n,dim
real(real64),intent(out)::x(n,dim)
type(halton_sampler)::h
integer::i,j
call h%init_faure();do j=1,dim;do i=1,n;x(i,j)=h%sample(j-1,int(i-1,int64));end do;end do
end subroutine
subroutine generate_halton_random_set(n,dim,x,seed)
integer,intent(in)::n,dim
real(real64),intent(out)::x(n,dim)
integer(int64),intent(in),optional::seed
type(halton_sampler)::h
integer::i,j
integer(int64)::s
s=0_int64;if(present(seed))s=seed
call h%init_random(s);do j=1,dim;do i=1,n;x(i,j)=h%sample(j-1,int(i-1,int64));end do;end do
end subroutine
end module spacefillr_halton
