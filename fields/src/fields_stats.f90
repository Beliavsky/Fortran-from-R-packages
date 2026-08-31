! GPL-2.0-or-later. Descriptive/binning utilities translated from fields.
module fields_stats
use fields_kinds, only: dp
use, intrinsic :: ieee_arithmetic, only: ieee_is_nan, ieee_value, ieee_quiet_nan
implicit none
private
public :: describe_values, stats_bin_result, stats_bin, quantile_type7
public :: fast_oneway_result, fast_oneway

type :: stats_bin_result
   real(dp), allocatable :: centers(:),breaks(:),stats(:,:)
end type

type :: fast_oneway_result
   integer :: n=0, ngroups=0
   integer, allocatable :: tags(:), group_index(:)
   real(dp), allocatable :: means(:,:), weight_sums(:), sse(:), mse(:)
end type fast_oneway_result

contains

function describe_values(x) result(s)
real(dp), intent(in) :: x(:)
real(dp) :: s(9)
real(dp), allocatable :: z(:)
integer :: n,i,k
real(dp) :: nanv,meanv
nanv=ieee_value(0.0_dp,ieee_quiet_nan)
n=count(.not.ieee_is_nan(x)); allocate(z(n)); k=0
do i=1,size(x)
   if(.not.ieee_is_nan(x(i))) then; k=k+1; z(k)=x(i); end if
end do
s=nanv; s(1)=real(n,dp); s(9)=real(size(x)-n,dp)
if(n==0) return
call sort_real(z)
meanv=sum(z)/real(n,dp); s(2)=meanv
if(n>1) s(3)=sqrt(sum((z-meanv)**2)/real(n-1,dp))
s(4)=z(1); s(5)=quantile_type7_sorted(z,0.25_dp); s(6)=quantile_type7_sorted(z,0.5_dp)
s(7)=quantile_type7_sorted(z,0.75_dp); s(8)=z(n)
end function describe_values

real(dp) function quantile_type7(x,p) result(q)
real(dp), intent(in) :: x(:),p
real(dp), allocatable :: z(:)
integer :: i,k,n
n=count(.not.ieee_is_nan(x)); if(n==0) then; q=ieee_value(0.0_dp,ieee_quiet_nan); return; end if
allocate(z(n)); k=0
do i=1,size(x); if(.not.ieee_is_nan(x(i))) then; k=k+1; z(k)=x(i); end if; end do
call sort_real(z); q=quantile_type7_sorted(z,p)
end function quantile_type7

pure real(dp) function quantile_type7_sorted(z,p) result(q)
real(dp), intent(in) :: z(:),p
real(dp) :: h,g
integer :: j,n
n=size(z); if(n==1) then; q=z(1); return; end if
if(p<=0.0_dp) then; q=z(1); return; end if
if(p>=1.0_dp) then; q=z(n); return; end if
h=1.0_dp+real(n-1,dp)*p; j=int(floor(h)); g=h-real(j,dp)
if(j>=n) then; q=z(n); else; q=(1.0_dp-g)*z(j)+g*z(j+1); end if
end function quantile_type7_sorted

function stats_bin(x,y,n_breaks,breaks) result(out)
real(dp), intent(in) :: x(:),y(:)
integer, intent(in), optional :: n_breaks
real(dp), intent(in), optional :: breaks(:)
type(stats_bin_result) :: out
real(dp), allocatable :: z(:)
integer :: nb,k,i,nz,nbr
if(size(x)/=size(y)) error stop 'stats_bin: x/y size mismatch'
if(present(breaks)) then
   if(size(breaks)<2) error stop 'stats_bin: need at least two break points'
   out%breaks=breaks
else
   nbr=10; if(present(n_breaks)) nbr=n_breaks
   if(nbr<2) error stop 'stats_bin: n_breaks must be >=2'
   allocate(out%breaks(nbr))
   do k=1,nbr; out%breaks(k)=minval(x,mask=.not.ieee_is_nan(x))+real(k-1,dp)/real(nbr-1,dp)* &
      (maxval(x,mask=.not.ieee_is_nan(x))-minval(x,mask=.not.ieee_is_nan(x))); end do
end if
nb=size(out%breaks)-1; allocate(out%centers(nb),out%stats(9,nb))
out%centers=0.5_dp*(out%breaks(:nb)+out%breaks(2:nb+1))
do k=1,nb
   nz=0
   do i=1,size(x)
      if(ieee_is_nan(x(i))) cycle
      if((k==1 .and. x(i)>=out%breaks(k) .and. x(i)<=out%breaks(k+1)) .or. &
         (k>1 .and. x(i)>out%breaks(k) .and. x(i)<=out%breaks(k+1))) nz=nz+1
   end do
   allocate(z(nz)); nz=0
   do i=1,size(x)
      if(ieee_is_nan(x(i))) cycle
      if((k==1 .and. x(i)>=out%breaks(k) .and. x(i)<=out%breaks(k+1)) .or. &
         (k>1 .and. x(i)>out%breaks(k) .and. x(i)<=out%breaks(k+1))) then
         nz=nz+1; z(nz)=y(i)
      end if
   end do
   out%stats(:,k)=describe_values(z); deallocate(z)
end do
end function stats_bin

function fast_oneway(level,y,weights) result(out)
integer,intent(in)::level(:)
real(dp),intent(in)::y(:,:)
real(dp),intent(in),optional::weights(:)
type(fast_oneway_result)::out
real(dp),allocatable::w(:)
integer::n,r,i,j,g,ng
logical::found
n=size(level);if(size(y,1)/=n)error stop 'fast_oneway: level/y mismatch'
r=size(y,2);allocate(w(n));w=1.0_dp
if(present(weights))then
 if(size(weights)/=n .or. any(weights<0.0_dp))error stop 'fast_oneway: invalid weights'
 w=weights
end if
allocate(out%tags(n),out%group_index(n));ng=0
do i=1,n
 found=.false.
 do g=1,ng
  if(level(i)==out%tags(g))then;out%group_index(i)=g;found=.true.;exit;end if
 end do
 if(.not.found)then;ng=ng+1;out%tags(ng)=level(i);out%group_index(i)=ng;end if
end do
out%ngroups=ng;out%n=n;out%tags=out%tags(:ng)
allocate(out%means(ng,r),out%weight_sums(ng),out%sse(r),out%mse(r))
out%means=0.0_dp;out%weight_sums=0.0_dp
do i=1,n
 g=out%group_index(i);out%weight_sums(g)=out%weight_sums(g)+w(i)
 out%means(g,:)=out%means(g,:)+w(i)*y(i,:)
end do
do g=1,ng
 if(out%weight_sums(g)>0.0_dp)out%means(g,:)=out%means(g,:)/out%weight_sums(g)
end do
out%sse=0.0_dp
do i=1,n
 g=out%group_index(i);out%sse=out%sse+w(i)*(y(i,:)-out%means(g,:))**2
end do
if(n>ng)then;out%mse=out%sse/real(n-ng,dp);else;out%mse=huge(1.0_dp);end if
end function fast_oneway

subroutine sort_real(a)
real(dp), intent(inout) :: a(:)
integer :: i,j
real(dp) :: key
do i=2,size(a)
   key=a(i); j=i-1
   do while(j>=1)
      if(a(j)<=key) exit
      a(j+1)=a(j); j=j-1
   end do
   a(j+1)=key
end do
end subroutine sort_real

end module fields_stats
