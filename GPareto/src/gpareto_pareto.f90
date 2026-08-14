! Modern Fortran translation of the computational core of GPareto 1.1.9.
! GPareto is GPL-3.0-only; see LICENSE and UPSTREAM.md.
module gpareto_pareto
  use gpareto_kinds, only : dp
  implicit none
  private
  public :: dominates, nondominated_mask, nondominated_points, nondominated_indices
  public :: nondom_set, squared_distances, hypervolume, hypervolume_improvement
contains
  pure logical function dominates(a,b) result(ok)
    real(dp), intent(in) :: a(:), b(:)
    ok = all(a <= b)
  end function dominates

  subroutine nondominated_mask(points, mask)
    real(dp), intent(in) :: points(:,:)
    logical, allocatable, intent(out) :: mask(:)
    integer :: i,j,n
    n=size(points,1)
    allocate(mask(n))
    mask=.true.
    do i=1,n
      do j=1,n
        if (j==i) cycle
        if (dominates(points(j,:),points(i,:)) .and. any(points(j,:) < points(i,:))) then
          mask(i)=.false.
          exit
        end if
      end do
    end do
  end subroutine nondominated_mask

  subroutine nondominated_indices(points, idx)
    real(dp), intent(in) :: points(:,:)
    integer, allocatable, intent(out) :: idx(:)
    logical, allocatable :: mask(:)
    integer :: i,k
    call nondominated_mask(points,mask)
    allocate(idx(count(mask)))
    k=0
    do i=1,size(mask)
      if(mask(i)) then
      k=k+1
      idx(k)=i
      end if
    end do
  end subroutine nondominated_indices

  subroutine nondominated_points(points, front)
    real(dp), intent(in) :: points(:,:)
    real(dp), allocatable, intent(out) :: front(:,:)
    integer, allocatable :: idx(:)
    integer :: i
    call nondominated_indices(points,idx)
    allocate(front(size(idx),size(points,2)))
    do i=1,size(idx)
    front(i,:)=points(idx(i),:)
    end do
  end subroutine nondominated_points

  subroutine nondom_set(points,ref,keep)
    real(dp), intent(in) :: points(:,:),ref(:,:)
    logical, allocatable, intent(out) :: keep(:)
    integer :: i,j
    allocate(keep(size(points,1)))
    keep=.true.
    do i=1,size(points,1)
      do j=1,size(ref,1)
        if (all(points(i,:) > ref(j,:))) then
          keep(i)=.false.
          exit
        end if
      end do
    end do
  end subroutine nondom_set

  subroutine squared_distances(x1,x2,d2)
    real(dp),intent(in)::x1(:,:),x2(:,:)
    real(dp),allocatable,intent(out)::d2(:,:)
    integer::i,j
    allocate(d2(size(x1,1),size(x2,1)))
    do j=1,size(x2,1)
    do i=1,size(x1,1)
      d2(i,j)=sum((x1(i,:)-x2(j,:))**2)
    end do
    end do
  end subroutine squared_distances

  recursive real(dp) function hypervolume(points,ref) result(hv)
    real(dp), intent(in) :: points(:,:),ref(:)
    real(dp), allocatable :: p(:,:),f(:,:),active(:,:)
    integer, allocatable :: ord(:)
    integer :: i,j,n,m,k,tmpi
    real(dp) :: z,nextz
    if(size(points,2)/=size(ref)) error stop 'hypervolume: dimension mismatch'
    m=size(ref)
    n=count([(all(points(i,:) < ref),i=1,size(points,1))])
    if(n==0) then
    hv=0.0_dp
    return
    end if
    allocate(p(n,m))
    k=0
    do i=1,size(points,1)
      if(all(points(i,:) < ref)) then
      k=k+1
      p(k,:)=points(i,:)
      end if
    end do
    call nondominated_points(p,f)
    n=size(f,1)
    if(m==1) then
      hv=max(0.0_dp,ref(1)-minval(f(:,1)))
      return
    end if
    allocate(ord(n))
    ord=[(i,i=1,n)]
    do i=2,n
      tmpi=ord(i)
      j=i-1
      do while(j>=1)
        if(f(ord(j),m)<=f(tmpi,m)) exit
        ord(j+1)=ord(j)
        j=j-1
      end do
      ord(j+1)=tmpi
    end do
    hv=0.0_dp
    do i=1,n
      z=f(ord(i),m)
      if(i<n) then
      nextz=min(f(ord(i+1),m),ref(m))
      else
      nextz=ref(m)
      end if
      if(nextz<=z) cycle
      allocate(active(i,m-1))
      do j=1,i
      active(j,:)=f(ord(j),1:m-1)
      end do
      hv=hv+(nextz-z)*hypervolume(active,ref(1:m-1))
      deallocate(active)
    end do
  end function hypervolume

  real(dp) function hypervolume_improvement(point,front,ref) result(v)
    real(dp),intent(in)::point(:),front(:,:),ref(:)
    real(dp),allocatable::aug(:,:)
    integer::n
    n=size(front,1)
    allocate(aug(n+1,size(front,2)))
    aug(1:n,:)=front
    aug(n+1,:)=point
    v=max(0.0_dp,hypervolume(aug,ref)-hypervolume(front,ref))
  end function hypervolume_improvement
end module gpareto_pareto
