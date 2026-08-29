module spam_distance
use spam_kinds, only: dp
use spam_types, only: csr_matrix
use spam_csr, only: csr_from_triplet
implicit none
private
public :: nearest_dist, rdist, rdist_earth
contains
function nearest_dist(x,y,method,delta,upper,p,radius,full) result(a)
real(dp),intent(in)::x(:,:)
real(dp),intent(in),optional::y(:,:)
character(len=*),intent(in),optional::method
real(dp),intent(in),optional::delta,p,radius
logical,intent(in),optional::upper,full
type(csr_matrix)::a
real(dp),allocatable::yy(:,:),val(:)
integer,allocatable::ir(:),jc(:)
integer::n1,n2,nd,i,j,k,nz,part,cap
real(dp)::del,pp,r,dd,tmp,rad,dotv
character(len=16)::meth
logical::same
n1=size(x,1)
nd=size(x,2)
same=.not.present(y)
if(same) then
   allocate(yy(n1,nd))
   yy=x
   n2=n1
else
   if(size(y,2)/=nd) error stop 'nearest_dist: x/y column mismatch'
   n2=size(y,1)
   allocate(yy(n2,nd))
   yy=y
end if
meth='euclidean'
if(present(method))meth=lower(trim(method))
del=1.0_dp
if(present(delta))del=abs(delta)
pp=2.0_dp
if(present(p))pp=p
part=0
if(same) then
   part=-1
   if(present(upper)) then
      if(upper) part=1
   end if
else if(present(upper)) then
   if(upper) part=1
end if
if(present(full))then
   if(full)part=0
end if
r=3963.34_dp
if(present(radius))r=radius
if(r<=0.0_dp) error stop 'nearest_dist: radius must be positive'
cap=max(16,min(n1*n2,1024))
allocate(ir(cap),jc(cap),val(cap))
nz=0
rad=acos(-1.0_dp)/180.0_dp
do i=1,n1
   do j=1,n2
      if(same .and. part<0 .and. j>i) cycle
      if(same .and. part>0 .and. j<i) cycle
      select case(trim(meth))
      case('euclidean','euclidian')
         tmp=0.0_dp
         do k=1,nd
         tmp=tmp+(x(i,k)-yy(j,k))**2
         if(tmp>del*del)exit
         end do
         if(tmp>del*del)cycle
         dd=sqrt(tmp)
      case('maximum')
         tmp=0.0_dp
         do k=1,nd
         tmp=max(tmp,abs(x(i,k)-yy(j,k)))
         if(tmp>del)exit
         end do
         if(tmp>del)cycle
         dd=tmp
      case('minkowski')
         if(pp<=0.0_dp) error stop 'nearest_dist: p must be positive'
         tmp=0.0_dp
         do k=1,nd
         tmp=tmp+abs(x(i,k)-yy(j,k))**pp
         if(tmp>del**pp)exit
         end do
         if(tmp>del**pp)cycle
         dd=tmp**(1.0_dp/pp)
      case('greatcircle')
         if(nd/=2) error stop 'nearest_dist: greatcircle requires two columns'
         if(del>180.1_dp) error stop 'nearest_dist: greatcircle delta must be <=180 degrees'
         dotv=cos(x(i,1)*rad)*cos(x(i,2)*rad)*cos(yy(j,1)*rad)*cos(yy(j,2)*rad) + &
              sin(x(i,1)*rad)*cos(x(i,2)*rad)*sin(yy(j,1)*rad)*cos(yy(j,2)*rad) + &
              sin(x(i,2)*rad)*sin(yy(j,2)*rad)
         dotv=max(-1.0_dp,min(1.0_dp,dotv))
         if(dotv < cos(del*rad))cycle
         if(dotv>=0.99999999999_dp)then
         dd=0.0_dp
         else
         dd=acos(dotv)*r
         end if
      case default
         error stop 'nearest_dist: invalid method'
      end select
      nz=nz+1
      call grow_triplet(ir,jc,val,nz)
      ir(nz)=i
      jc(nz)=j
      val(nz)=dd
   end do
end do
a=csr_from_triplet(n1,n2,ir(:nz),jc(:nz),val(:nz),eps=-1.0_dp)
end function nearest_dist

function rdist(x,y,delta) result(a)
real(dp),intent(in)::x(:,:),y(:,:)
real(dp),intent(in),optional::delta
type(csr_matrix)::a
if(present(delta))then
a=nearest_dist(x,y,delta=delta)
else
a=nearest_dist(x,y)
end if
end function rdist

function rdist_earth(x,y,delta,miles,radius) result(a)
real(dp),intent(in)::x(:,:),y(:,:)
real(dp),intent(in),optional::delta,radius
logical,intent(in),optional::miles
type(csr_matrix)::a
real(dp)::d,r
logical::mi
d=1.0_dp
if(present(delta))d=delta
mi=.true.
if(present(miles))mi=miles
if(mi)then
r=3963.34_dp
else
r=6378.388_dp
end if
if(present(radius))r=radius
a=nearest_dist(x,y,'greatcircle',d,p=2.0_dp,radius=r)
end function rdist_earth

subroutine grow_triplet(ir,jc,val,need)
integer,allocatable,intent(inout)::ir(:),jc(:)
real(dp),allocatable,intent(inout)::val(:)
integer,intent(in)::need
integer,allocatable::ii(:),jj(:)
real(dp),allocatable::vv(:)
integer::n
if(need<=size(ir))return
n=max(need,2*size(ir))
allocate(ii(n),jj(n),vv(n))
ii=0
jj=0
vv=0.0_dp
ii(:size(ir))=ir
jj(:size(jc))=jc
vv(:size(val))=val
call move_alloc(ii,ir)
call move_alloc(jj,jc)
call move_alloc(vv,val)
end subroutine grow_triplet

pure function lower(s) result(t)
character(len=*),intent(in)::s
character(len=len(s))::t
integer::i,c
t=s
do i=1,len(s)
c=iachar(t(i:i))
if(c>=65.and.c<=90)t(i:i)=achar(c+32)
end do
end function lower
end module spam_distance
