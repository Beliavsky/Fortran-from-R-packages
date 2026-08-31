! GPL-2.0-or-later. Grid variogram translated from fields::vgram.matrix.
module fields_variogram
use fields_kinds, only: dp
use, intrinsic :: ieee_arithmetic, only: ieee_is_nan, ieee_value, ieee_quiet_nan
implicit none
private
public :: variogram_grid_result, variogram_grid

type :: variogram_grid_result
   real(dp) :: dx=1.0_dp,dy=1.0_dp
   integer, allocatable :: shift(:,:)
   real(dp), allocatable :: distance_full(:),variogram_full(:),robust_full(:)
   integer, allocatable :: count_full(:)
   real(dp), allocatable :: distance(:),variogram(:)
end type

contains

function variogram_grid(dat,dx,dy,radius) result(out)
real(dp), intent(in) :: dat(:,:)
real(dp), intent(in), optional :: dx,dy,radius
type(variogram_grid_result) :: out
real(dp) :: ddx,ddy,rr,d,nanv,sum2,sumroot
integer :: m,n,imax,jmax,di,dj,i,j,ii,jj,k,nb,ng,g,cnt
integer, allocatable :: si(:),sj(:),ord(:)
real(dp), allocatable :: dist(:),vg(:),rvg(:),uniq(:),num(:),den(:)
nanv=ieee_value(0.0_dp,ieee_quiet_nan)
ddx=1.0_dp; if(present(dx)) ddx=dx
ddy=1.0_dp; if(present(dy)) ddy=dy
rr=5.0_dp*max(ddx,ddy); if(present(radius)) rr=radius
m=size(dat,1); n=size(dat,2); imax=min(nint(rr/ddx),m); jmax=min(nint(rr/ddy),n)
nb=0
do dj=0,jmax
   do di=-imax,imax
      if(di==0 .and. dj==0) cycle
      if(dj==0 .and. di<=0) cycle
      d=sqrt((ddx*real(di,dp))**2+(ddy*real(dj,dp))**2)
      if(d>0.0_dp .and. d<=rr) nb=nb+1
   end do
end do
allocate(si(nb),sj(nb),dist(nb),vg(nb),rvg(nb),out%count_full(nb),ord(nb)); k=0
do dj=0,jmax
   do di=-imax,imax
      if(di==0 .and. dj==0) cycle
      if(dj==0 .and. di<=0) cycle
      d=sqrt((ddx*real(di,dp))**2+(ddy*real(dj,dp))**2)
      if(d<=0.0_dp .or. d>rr) cycle
      k=k+1; si(k)=di; sj(k)=dj; dist(k)=d; ord(k)=k
      sum2=0.0_dp; sumroot=0.0_dp; cnt=0
      do j=1,n
         jj=j+dj; if(jj<1 .or. jj>n) cycle
         do i=1,m
            ii=i+di; if(ii<1 .or. ii>m) cycle
            if(ieee_is_nan(dat(i,j)) .or. ieee_is_nan(dat(ii,jj))) cycle
            d=dat(i,j)-dat(ii,jj); sum2=sum2+0.5_dp*d*d; sumroot=sumroot+sqrt(abs(d)); cnt=cnt+1
         end do
      end do
      out%count_full(k)=cnt
      if(cnt>0) then
         vg(k)=sum2/real(cnt,dp); rvg(k)=0.5_dp*(sumroot/real(cnt,dp))**4/(0.457_dp+0.494_dp*real(cnt,dp))
      else
         vg(k)=nanv; rvg(k)=nanv
      end if
   end do
end do
call sort_by_distance(dist,si,sj,vg,rvg,out%count_full)
allocate(out%shift(nb,2)); out%shift(:,1)=si; out%shift(:,2)=sj
out%distance_full=dist; out%variogram_full=vg; out%robust_full=rvg; out%dx=ddx; out%dy=ddy
! Collapse exactly/common-to-roundoff equal distances using weighted means.
allocate(uniq(nb),num(nb),den(nb)); ng=0
k=1
do while(k<=nb)
   ng=ng+1; uniq(ng)=dist(k); num(ng)=0.0_dp; den(ng)=0.0_dp
   g=k
   do while(g<=nb)
      if(abs(dist(g)-dist(k))>64.0_dp*epsilon(1.0_dp)*max(1.0_dp,dist(k))) exit
      if(.not.ieee_is_nan(vg(g))) then
         num(ng)=num(ng)+vg(g)*real(out%count_full(g),dp); den(ng)=den(ng)+real(out%count_full(g),dp)
      end if
      g=g+1
   end do
   k=g
end do
allocate(out%distance(ng),out%variogram(ng)); out%distance=uniq(:ng)
do k=1,ng
   if(den(k)>0.0_dp) then; out%variogram(k)=num(k)/den(k); else; out%variogram(k)=nanv; end if
end do
end function variogram_grid

subroutine sort_by_distance(d,si,sj,vg,rvg,cnt)
real(dp), intent(inout) :: d(:),vg(:),rvg(:)
integer, intent(inout) :: si(:),sj(:),cnt(:)
integer :: i,j,ki,kj,kc
real(dp) :: kd,kv,kr
do i=2,size(d)
   kd=d(i); kv=vg(i); kr=rvg(i); ki=si(i); kj=sj(i); kc=cnt(i); j=i-1
   do while(j>=1)
      if(d(j)<=kd) exit
      d(j+1)=d(j); vg(j+1)=vg(j); rvg(j+1)=rvg(j); si(j+1)=si(j); sj(j+1)=sj(j); cnt(j+1)=cnt(j); j=j-1
   end do
   d(j+1)=kd; vg(j+1)=kv; rvg(j+1)=kr; si(j+1)=ki; sj(j+1)=kj; cnt(j+1)=kc
end do
end subroutine sort_by_distance

end module fields_variogram
