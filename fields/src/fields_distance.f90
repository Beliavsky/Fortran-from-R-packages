! GPL-2.0-or-later. Translation of fields distance kernels.
module fields_distance
use fields_kinds, only: dp
use spam_types, only: csr_matrix
use spam_distance, only: nearest_dist
implicit none
private
public :: fields_rdist, fields_rdist_self, fields_rdist_vec, fields_rdist_earth, &
          fields_rdist_earth_vec, fields_rdist_near, compact_to_matrix, add_to_diagonal

contains

function fields_rdist(x1,x2) result(d)
real(dp), intent(in) :: x1(:,:),x2(:,:)
real(dp), allocatable :: d(:,:)
integer :: i,j,k,n1,n2,nd
real(dp) :: s
n1=size(x1,1); n2=size(x2,1); nd=size(x1,2)
if(size(x2,2)/=nd) error stop 'fields_rdist: coordinate dimensions differ'
allocate(d(n1,n2))
do j=1,n2
   do i=1,n1
      s=0.0_dp
      do k=1,nd
         s=s+(x1(i,k)-x2(j,k))**2
      end do
      d(i,j)=sqrt(s)
   end do
end do
end function fields_rdist

function fields_rdist_self(x) result(d)
real(dp), intent(in) :: x(:,:)
real(dp), allocatable :: d(:,:)
integer :: i,j,k,n,nd
real(dp) :: s
n=size(x,1); nd=size(x,2); allocate(d(n,n)); d=0.0_dp
do j=1,n
   do i=j+1,n
      s=0.0_dp
      do k=1,nd; s=s+(x(i,k)-x(j,k))**2; end do
      d(i,j)=sqrt(s); d(j,i)=d(i,j)
   end do
end do
end function fields_rdist_self

function fields_rdist_vec(x1,x2) result(d)
real(dp), intent(in) :: x1(:,:),x2(:,:)
real(dp), allocatable :: d(:)
integer :: i,k,n,nd
if(any(shape(x1)/=shape(x2))) error stop 'fields_rdist_vec: dimensions differ'
n=size(x1,1); nd=size(x1,2); allocate(d(n)); d=0.0_dp
do i=1,n
   do k=1,nd; d(i)=d(i)+(x1(i,k)-x2(i,k))**2; end do
   d(i)=sqrt(d(i))
end do
end function fields_rdist_vec

function fields_rdist_earth(x1,x2,miles,radius) result(d)
real(dp), intent(in) :: x1(:,:),x2(:,:)
logical, intent(in), optional :: miles
real(dp), intent(in), optional :: radius
real(dp), allocatable :: d(:,:)
real(dp) :: r,torad,lon1,lat1,lon2,lat2,a
logical :: mi
integer :: i,j
if(size(x1,2)/=2 .or. size(x2,2)/=2) error stop 'fields_rdist_earth: lon/lat matrices require two columns'
mi=.true.; if(present(miles)) mi=miles
if(mi) then; r=3963.34_dp; else; r=6378.388_dp; end if
if(present(radius)) r=radius
if(r<=0.0_dp) error stop 'fields_rdist_earth: radius must be positive'
torad=acos(-1.0_dp)/180.0_dp
allocate(d(size(x1,1),size(x2,1)))
do j=1,size(x2,1)
   lon2=x2(j,1)*torad; lat2=x2(j,2)*torad
   do i=1,size(x1,1)
      lon1=x1(i,1)*torad; lat1=x1(i,2)*torad
      a=sin(0.5_dp*(lat2-lat1))**2 + cos(lat1)*cos(lat2)*sin(0.5_dp*(lon2-lon1))**2
      a=max(0.0_dp,min(1.0_dp,a))
      d(i,j)=2.0_dp*atan2(sqrt(a),sqrt(max(0.0_dp,1.0_dp-a)))*r
   end do
end do
end function fields_rdist_earth

function fields_rdist_earth_vec(x1,x2,miles,radius) result(d)
real(dp), intent(in) :: x1(:,:),x2(:,:)
logical, intent(in), optional :: miles
real(dp), intent(in), optional :: radius
real(dp), allocatable :: d(:),dm(:,:)
integer :: i,n
if(any(shape(x1)/=shape(x2))) error stop 'fields_rdist_earth_vec: dimensions differ'
n=size(x1,1); dm=fields_rdist_earth(x1,x2,miles,radius); allocate(d(n))
do i=1,n; d(i)=dm(i,i); end do
end function fields_rdist_earth_vec

function fields_rdist_near(x1,x2,delta,method,p,radius,upper,full) result(a)
real(dp), intent(in) :: x1(:,:),x2(:,:),delta
character(len=*), intent(in), optional :: method
real(dp), intent(in), optional :: p,radius
logical, intent(in), optional :: upper,full
type(csr_matrix) :: a
character(len=16) :: meth
real(dp) :: pp,rr
logical :: up,fu
meth='euclidean'; if(present(method)) meth=method
pp=2.0_dp; if(present(p)) pp=p
rr=3963.34_dp; if(present(radius)) rr=radius
up=.false.; if(present(upper)) up=upper
fu=.false.; if(present(full)) fu=full
a=nearest_dist(x1,x2,method=meth,delta=delta,p=pp,radius=rr,upper=up,full=fu)
end function fields_rdist_near

function compact_to_matrix(v,n,diag_value,lower,upper) result(a)
real(dp), intent(in) :: v(:)
integer, intent(in) :: n
real(dp), intent(in), optional :: diag_value
logical, intent(in), optional :: lower,upper
real(dp), allocatable :: a(:,:)
real(dp) :: dv
logical :: lo,up
integer :: i,j,k
if(size(v)/=n*(n-1)/2) error stop 'compact_to_matrix: wrong compact vector length'
dv=0.0_dp; if(present(diag_value)) dv=diag_value
lo=.false.; if(present(lower)) lo=lower
up=.true.; if(present(upper)) up=upper
allocate(a(n,n)); a=0.0_dp; k=0
! R dist order: columns of lower triangle: (2,1),(3,1),...,(n,1),(3,2),...
do j=1,n-1
   do i=j+1,n
      k=k+1
      if(up) a(j,i)=v(k)
      if(lo) a(i,j)=v(k)
   end do
end do
do i=1,n; a(i,i)=dv; end do
end function compact_to_matrix

subroutine add_to_diagonal(a,v)
real(dp), intent(inout) :: a(:,:)
real(dp), intent(in) :: v(:)
integer :: i,n
n=size(a,1)
if(size(a,2)/=n .or. size(v)/=n) error stop 'add_to_diagonal: dimension mismatch'
do i=1,n; a(i,i)=a(i,i)+v(i); end do
end subroutine add_to_diagonal

end module fields_distance
