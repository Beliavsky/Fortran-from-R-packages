! GPL-2.0-or-later. Geometry helpers from fields.
module fields_geometry
use fields_kinds, only: dp
use fields_distance, only: fields_rdist
implicit none
private
public :: points_in_polygon, minimax_criterion, greedy_cover_design

contains

function points_in_polygon(points,polygon) result(inside)
real(dp), intent(in) :: points(:,:),polygon(:,:)
logical, allocatable :: inside(:)
integer :: i,j,k,nv
real(dp) :: x,y,x1,y1,x2,y2,xint
logical :: c
if(size(points,2)/=2 .or. size(polygon,2)/=2 .or. size(polygon,1)<3) error stop 'points_in_polygon: requires 2D polygon'
nv=size(polygon,1); allocate(inside(size(points,1)))
do k=1,size(points,1)
   x=points(k,1); y=points(k,2); c=.false.; j=nv
   do i=1,nv
      x1=polygon(i,1); y1=polygon(i,2); x2=polygon(j,1); y2=polygon(j,2)
      if(point_on_segment(x,y,x1,y1,x2,y2)) then; c=.true.; exit; end if
      if((y1>y) .neqv. (y2>y)) then
         xint=(x2-x1)*(y-y1)/(y2-y1)+x1
         if(x<xint) c=.not.c
      end if
      j=i
   end do
   inside(k)=c
end do
end function points_in_polygon

pure logical function point_on_segment(x,y,x1,y1,x2,y2) result(on)
real(dp),intent(in)::x,y,x1,y1,x2,y2
real(dp)::cross,tol
tol=64.0_dp*epsilon(1.0_dp)*max(1.0_dp,abs(x),abs(y),abs(x1),abs(y1),abs(x2),abs(y2))
cross=(x-x1)*(y2-y1)-(y-y1)*(x2-x1)
on=abs(cross)<=tol .and. x>=min(x1,x2)-tol .and. x<=max(x1,x2)+tol .and. y>=min(y1,y2)-tol .and. y<=max(y1,y2)+tol
end function point_on_segment

real(dp) function minimax_criterion(candidates,design_index) result(v)
real(dp), intent(in) :: candidates(:,:)
integer, intent(in) :: design_index(:)
logical, allocatable :: isdes(:)
real(dp), allocatable :: des(:,:),cand(:,:),d(:,:)
integer :: i,j,nc
if(any(design_index<1) .or. any(design_index>size(candidates,1))) error stop 'minimax_criterion: bad index'
allocate(isdes(size(candidates,1))); isdes=.false.; isdes(design_index)=.true.
allocate(des(size(design_index),size(candidates,2))); des=candidates(design_index,:)
nc=count(.not.isdes); if(nc==0) then; v=0.0_dp; return; end if
allocate(cand(nc,size(candidates,2))); j=0
do i=1,size(candidates,1); if(.not.isdes(i)) then; j=j+1; cand(j,:)=candidates(i,:); end if; end do
d=fields_rdist(cand,des); v=0.0_dp
do i=1,nc; v=max(v,minval(d(i,:))); end do
end function minimax_criterion

function greedy_cover_design(candidates,nd,start_index) result(index)
! A deterministic farthest-point design: the same minimax objective used by cover.design,
! but without its stochastic interchange refinement.
real(dp), intent(in) :: candidates(:,:)
integer, intent(in) :: nd
integer, intent(in), optional :: start_index
integer, allocatable :: index(:)
real(dp), allocatable :: dmin(:),d(:,:)
logical, allocatable :: chosen(:)
integer :: n,k,j,seed
n=size(candidates,1); if(nd<1 .or. nd>n) error stop 'greedy_cover_design: invalid nd'
allocate(index(nd),chosen(n),dmin(n)); chosen=.false.; dmin=huge(1.0_dp)
seed=1; if(present(start_index)) seed=start_index
if(seed<1 .or. seed>n) error stop 'greedy_cover_design: invalid start index'
index(1)=seed; chosen(seed)=.true.
do k=2,nd
   d=fields_rdist(candidates,candidates(index(k-1:k-1),:)); dmin=min(dmin,d(:,1)); where(chosen) dmin=-1.0_dp
   j=maxloc(dmin,dim=1); index(k)=j; chosen(j)=.true.
end do
end function greedy_cover_design

end module fields_geometry
