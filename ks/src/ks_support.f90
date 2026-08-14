! SPDX-License-Identifier: GPL-2.0-only
module ks_support
  use ks_kinds, only: dp
  implicit none
  private
  public :: convex_hull_2d, contour_support_points
contains
  pure real(dp) function cross_xy(o1,o2,a1,a2,b1,b2) result(v)
    real(dp),intent(in)::o1,o2,a1,a2,b1,b2
    v=(a1-o1)*(b2-o2)-(a2-o2)*(b1-o1)
  end function

  subroutine convex_hull_2d(points,hull)
    real(dp),intent(in)::points(:,:)
    real(dp),allocatable,intent(out)::hull(:,:)
    integer,allocatable::idx(:),stack(:)
    integer::n,i,j,k,t
    if(size(points,2)/=2)error stop 'convex_hull_2d: requires 2D'
    n=size(points,1)
    if(n<=2)then;allocate(hull(n,2));hull=points;return;end if
    allocate(idx(n));idx=[(i,i=1,n)]
    do i=2,n
      k=idx(i);j=i-1
      do while(j>=1)
        if(points(idx(j),1)<points(k,1))exit
        if(points(idx(j),1)<=points(k,1).and.points(idx(j),2)<=points(k,2))exit
        idx(j+1)=idx(j);j=j-1
      end do
      idx(j+1)=k
    end do
    allocate(stack(2*n));t=0
    do i=1,n
      do while(t>=2)
        if(cross_xy(points(stack(t-1),1),points(stack(t-1),2), &
                     points(stack(t),1),points(stack(t),2),points(idx(i),1),points(idx(i),2))>0.0_dp)exit
        t=t-1
      end do
      t=t+1;stack(t)=idx(i)
    end do
    k=t
    do i=n-1,1,-1
      do while(t>=k+1)
        if(cross_xy(points(stack(t-1),1),points(stack(t-1),2), &
                     points(stack(t),1),points(stack(t),2),points(idx(i),1),points(idx(i),2))>0.0_dp)exit
        t=t-1
      end do
      t=t+1;stack(t)=idx(i)
    end do
    if(t>1)t=t-1
    allocate(hull(t,2));do i=1,t;hull(i,:)=points(stack(i),:);end do
  end subroutine

  subroutine contour_support_points(grid,density,level,points,convex_hull)
    real(dp),intent(in)::grid(:,:),density(:),level
    real(dp),allocatable,intent(out)::points(:,:)
    logical,intent(in),optional::convex_hull
    real(dp),allocatable::raw(:,:),h(:,:)
    logical::ch
    integer::i,k,n
    if(size(grid,1)/=size(density))error stop 'contour_support_points: shape'
    n=count(density>level);allocate(raw(n,size(grid,2)));k=0
    do i=1,size(density);if(density(i)>level)then;k=k+1;raw(k,:)=grid(i,:);end if;end do
    ch=.true.;if(present(convex_hull))ch=convex_hull
    if(ch.and.size(grid,2)==2.and.n>=3)then;call convex_hull_2d(raw,h);call move_alloc(h,points)
    else;call move_alloc(raw,points);end if
  end subroutine
end module ks_support
