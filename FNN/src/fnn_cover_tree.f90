! FNN-fortran: modern Fortran translation of computational code from FNN 1.1.4.1.
! Modified/translated 2026 by the FNN-fortran contributors.
! SPDX-License-Identifier: GPL-2.0-or-later
! See UPSTREAM.md and upstream/FNN-1.1.4.1 for original authorship and notices.
module fnn_cover_tree
  use fnn_kinds, only : dp
  implicit none
  private
  public :: cover_search_self, cover_search_query
  real(dp), parameter :: cover_base = 1.3_dp

  type :: cover_tree
    real(dp), pointer :: x(:,:) => null()
    integer, allocatable :: parent(:), first_child(:), next_sibling(:), level(:)
    real(dp), allocatable :: radius(:)
    integer :: root = 0
  end type cover_tree
contains
  subroutine cover_search_self(data,k,index,distance)
    real(dp), target, intent(in) :: data(:,:)
    integer, intent(in) :: k
    integer, intent(out) :: index(:,:)
    real(dp), intent(out) :: distance(:,:)
    type(cover_tree) :: tree
    integer :: i
    call build_cover_tree(tree,data)
    do i=1,size(data,1)
      call query_cover(tree,data(i,:),k,index(i,:),distance(i,:),i)
    end do
  end subroutine cover_search_self

  subroutine cover_search_query(data,query,k,index,distance)
    real(dp), target, intent(in) :: data(:,:),query(:,:)
    integer, intent(in) :: k
    integer, intent(out) :: index(:,:)
    real(dp), intent(out) :: distance(:,:)
    type(cover_tree) :: tree
    integer :: i
    call build_cover_tree(tree,data)
    do i=1,size(query,1)
      call query_cover(tree,query(i,:),k,index(i,:),distance(i,:),0)
    end do
  end subroutine cover_search_query

  subroutine build_cover_tree(tree,data)
    type(cover_tree), intent(out) :: tree
    real(dp), target, intent(in) :: data(:,:)
    integer :: n,i,lvl,p,a
    real(dp) :: maxd
    n=size(data,1); tree%x=>data
    allocate(tree%parent(max(1,n)),tree%first_child(max(1,n)),tree%next_sibling(max(1,n)))
    allocate(tree%level(max(1,n)),tree%radius(max(1,n)))
    tree%parent=0; tree%first_child=0; tree%next_sibling=0; tree%level=0; tree%radius=0.0_dp
    if(n==0) return
    tree%root=1
    if(n==1) return
    maxd=0.0_dp
    do i=2,n
      maxd=max(maxd,euclidean(data(1,:),data(i,:)))
    end do
    if(maxd>0.0_dp) then
      lvl=ceiling(log(maxd)/log(cover_base))
      do while(cover_base**real(lvl,dp)<maxd)
        lvl=lvl+1
      end do
    else
      lvl=0
    end if
    tree%level(1)=lvl
    do i=2,n
      call insert_cover(tree,1,i,lvl)
    end do
    do i=1,n
      p=tree%parent(i)
      do while(p/=0)
        a=p
        tree%radius(a)=max(tree%radius(a),euclidean(data(a,:),data(i,:)))
        p=tree%parent(p)
      end do
    end do
  end subroutine build_cover_tree

  recursive subroutine insert_cover(tree,node,newpoint,level)
    type(cover_tree), intent(inout) :: tree
    integer, intent(in) :: node,newpoint,level
    integer :: child,best_child
    real(dp) :: threshold,d,best
    threshold=cover_base**real(level-1,dp)
    child=tree%first_child(node); best_child=0; best=huge(1.0_dp)
    do while(child/=0)
      d=euclidean(tree%x(child,:),tree%x(newpoint,:))
      if(d<=threshold .and. d<best) then
        best=d; best_child=child
      end if
      child=tree%next_sibling(child)
    end do
    if(best_child/=0) then
      call insert_cover(tree,best_child,newpoint,level-1)
    else
      tree%parent(newpoint)=node
      tree%level(newpoint)=level-1
      tree%next_sibling(newpoint)=tree%first_child(node)
      tree%first_child(node)=newpoint
    end if
  end subroutine insert_cover

  subroutine query_cover(tree,q,k,index,distance,exclude)
    type(cover_tree), intent(in) :: tree
    real(dp), intent(in) :: q(:)
    integer, intent(in) :: k,exclude
    integer, intent(out) :: index(:)
    real(dp), intent(out) :: distance(:)
    if(size(index)/=k .or. size(distance)/=k) error stop "query_cover: result size mismatch"
    index=-1; distance=huge(1.0_dp)
    if(tree%root/=0) call search_cover_node(tree,tree%root,q,exclude,index,distance)
  end subroutine query_cover

  recursive subroutine search_cover_node(tree,node,q,exclude,index,distance)
    type(cover_tree), intent(in) :: tree
    integer, intent(in) :: node,exclude
    real(dp), intent(in) :: q(:)
    integer, intent(inout) :: index(:)
    real(dp), intent(inout) :: distance(:)
    integer :: child
    real(dp) :: d,dc,lower
    d=euclidean(tree%x(node,:),q)
    if(node/=exclude) call best_insert(node,d,index,distance)
    child=tree%first_child(node)
    do while(child/=0)
      dc=euclidean(tree%x(child,:),q)
      lower=max(0.0_dp,dc-tree%radius(child))
      if(lower<=distance(size(distance))) then
        call search_cover_node(tree,child,q,exclude,index,distance)
      end if
      child=tree%next_sibling(child)
    end do
  end subroutine search_cover_node

  pure real(dp) function euclidean(a,b) result(d)
    real(dp), intent(in) :: a(:),b(:)
    d=sqrt(max(0.0_dp,sum((a-b)**2)))
  end function euclidean

  subroutine best_insert(id,d,index,distance)
    integer, intent(in) :: id
    real(dp), intent(in) :: d
    integer, intent(inout) :: index(:)
    real(dp), intent(inout) :: distance(:)
    integer :: i,j,pos
    pos=size(distance)+1
    do i=1,size(distance)
      if(d<distance(i)) then
        pos=i; exit
      else if(d>distance(i)) then
        cycle
      else if(index(i)<0 .or. id<index(i)) then
        pos=i; exit
      end if
    end do
    if(pos>size(distance)) return
    do j=size(distance),pos+1,-1
      distance(j)=distance(j-1); index(j)=index(j-1)
    end do
    distance(pos)=d; index(pos)=id
  end subroutine best_insert
end module fnn_cover_tree
