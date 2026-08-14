! FNN-fortran: modern Fortran translation of computational code from FNN 1.1.4.1.
! Modified/translated 2026 by the FNN-fortran contributors.
! SPDX-License-Identifier: GPL-2.0-or-later
! See UPSTREAM.md and upstream/FNN-1.1.4.1 for original authorship and notices.
module fnn_neighbors
  use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
  use fnn_kinds, only : dp
  use fnn_types, only : knn_result
  use fnn_cover_tree, only : cover_search_self, cover_search_query
  implicit none
  private
  public :: get_knn, get_knnx, knn_index, knn_dist, knnx_index, knnx_dist
  public :: mean_log_knn_distance

  type :: kd_tree
    real(dp), pointer :: x(:,:) => null()
    integer, allocatable :: point(:), axis(:), left(:), right(:)
    integer :: root = 0
    integer :: used = 0
  end type kd_tree
contains
  function get_knn(data, k, algorithm) result(res)
    real(dp), target, intent(in) :: data(:,:)
    integer, intent(in) :: k
    character(len=*), intent(in), optional :: algorithm
    type(knn_result) :: res
    character(len=24) :: alg
    integer :: n, keff, i
    type(kd_tree) :: tree

    n = size(data,1)
    if (k < 1) error stop "get_knn: k must be positive"
    allocate(res%index(n,k), res%distance(n,k))
    res%index = -1
    res%distance = ieee_value(0.0_dp, ieee_quiet_nan)
    if (n <= 1) return
    keff = min(k,n-1)
    alg = normalized_algorithm(algorithm)
    select case(trim(alg))
    case("cr")
      do i=1,n
        call query_correlation(data, data(i,:), keff, res%index(i,1:keff), res%distance(i,1:keff), i)
      end do
    case("brute")
      do i=1,n
        call query_brute(data, data(i,:), keff, res%index(i,1:keff), res%distance(i,1:keff), i)
      end do
    case("kd_tree")
      call build_kd_tree(tree,data)
      do i=1,n
        call query_kd(tree,data(i,:),keff,res%index(i,1:keff),res%distance(i,1:keff),i)
      end do
    case("cover_tree")
      call cover_search_self(data,keff,res%index(:,1:keff),res%distance(:,1:keff))
    case default
      error stop "get_knn: unknown algorithm"
    end select
  end function get_knn

  function get_knnx(data, query, k, algorithm) result(res)
    real(dp), target, intent(in) :: data(:,:), query(:,:)
    integer, intent(in) :: k
    character(len=*), intent(in), optional :: algorithm
    type(knn_result) :: res
    character(len=24) :: alg
    integer :: n, m, keff, i
    type(kd_tree) :: tree

    if (size(data,2) /= size(query,2)) error stop "get_knnx: dimensions differ"
    if (k < 1) error stop "get_knnx: k must be positive"
    n=size(data,1); m=size(query,1); keff=min(k,n)
    allocate(res%index(m,k),res%distance(m,k))
    res%index=-1
    res%distance=ieee_value(0.0_dp,ieee_quiet_nan)
    if (n == 0) return
    alg=normalized_algorithm(algorithm)
    select case(trim(alg))
    case("cr")
      do i=1,m
        call query_correlation(data,query(i,:),keff,res%index(i,1:keff),res%distance(i,1:keff),0)
      end do
    case("brute")
      do i=1,m
        call query_brute(data,query(i,:),keff,res%index(i,1:keff),res%distance(i,1:keff),0)
      end do
    case("kd_tree")
      call build_kd_tree(tree,data)
      do i=1,m
        call query_kd(tree,query(i,:),keff,res%index(i,1:keff),res%distance(i,1:keff),0)
      end do
    case("cover_tree")
      call cover_search_query(data,query,keff,res%index(:,1:keff),res%distance(:,1:keff))
    case default
      error stop "get_knnx: unknown algorithm"
    end select
  end function get_knnx

  function knn_index(data,k,algorithm) result(idx)
    real(dp), target, intent(in) :: data(:,:)
    integer, intent(in) :: k
    character(len=*), intent(in), optional :: algorithm
    integer, allocatable :: idx(:,:)
    type(knn_result) :: z
    z=get_knn(data,k,algorithm); idx=z%index
  end function knn_index

  function knn_dist(data,k,algorithm) result(dist)
    real(dp), target, intent(in) :: data(:,:)
    integer, intent(in) :: k
    character(len=*), intent(in), optional :: algorithm
    real(dp), allocatable :: dist(:,:)
    type(knn_result) :: z
    z=get_knn(data,k,algorithm); dist=z%distance
  end function knn_dist

  function knnx_index(data,query,k,algorithm) result(idx)
    real(dp), target, intent(in) :: data(:,:), query(:,:)
    integer, intent(in) :: k
    character(len=*), intent(in), optional :: algorithm
    integer, allocatable :: idx(:,:)
    type(knn_result) :: z
    z=get_knnx(data,query,k,algorithm); idx=z%index
  end function knnx_index

  function knnx_dist(data,query,k,algorithm) result(dist)
    real(dp), target, intent(in) :: data(:,:), query(:,:)
    integer, intent(in) :: k
    character(len=*), intent(in), optional :: algorithm
    real(dp), allocatable :: dist(:,:)
    type(knn_result) :: z
    z=get_knnx(data,query,k,algorithm); dist=z%distance
  end function knnx_dist

  function mean_log_knn_distance(data,k,algorithm) result(mld)
    real(dp), target, intent(in) :: data(:,:)
    integer, intent(in) :: k
    character(len=*), intent(in), optional :: algorithm
    real(dp), allocatable :: mld(:)
    type(knn_result) :: z
    integer :: j
    if (k >= size(data,1)) error stop "mean_log_knn_distance: k must be < n"
    z=get_knn(data,k,algorithm)
    allocate(mld(k))
    do j=1,k
      mld(j)=sum(log(z%distance(:,j)))/real(size(data,1),dp)
    end do
  end function mean_log_knn_distance

  function normalized_algorithm(algorithm) result(alg)
    character(len=*), intent(in), optional :: algorithm
    character(len=24) :: alg
    integer :: i,c
    alg="kd_tree"
    if (present(algorithm)) then
      alg=adjustl(algorithm)
      do i=1,len_trim(alg)
        c=iachar(alg(i:i)); if(c>=65 .and. c<=90) alg(i:i)=achar(c+32)
      end do
    end if
  end function normalized_algorithm

  subroutine query_brute(data,q,k,idx,dist,exclude)
    real(dp), intent(in) :: data(:,:), q(:)
    integer, intent(in) :: k, exclude
    integer, intent(out) :: idx(:)
    real(dp), intent(out) :: dist(:)
    integer :: j
    real(dp) :: d2
    if(size(idx)/=k .or. size(dist)/=k) error stop "query_brute: result size mismatch"
    call init_best(idx,dist)
    do j=1,size(data,1)
      if (j == exclude) cycle
      d2=sum((data(j,:)-q)**2)
      call best_insert(j,sqrt(max(0.0_dp,d2)),idx,dist)
    end do
  end subroutine query_brute

  subroutine query_correlation(data,q,k,idx,dist,exclude)
    real(dp), intent(in) :: data(:,:), q(:)
    integer, intent(in) :: k, exclude
    integer, intent(out) :: idx(:)
    real(dp), intent(out) :: dist(:)
    integer :: j
    real(dp) :: d
    if(size(idx)/=k .or. size(dist)/=k) error stop "query_correlation: result size mismatch"
    call init_best(idx,dist)
    do j=1,size(data,1)
      if(j==exclude) cycle
      d=1.0_dp-dot_product(data(j,:),q)
      call best_insert(j,d,idx,dist)
    end do
  end subroutine query_correlation

  subroutine build_kd_tree(tree,data)
    type(kd_tree), intent(out) :: tree
    real(dp), target, intent(in) :: data(:,:)
    integer, allocatable :: ids(:)
    integer :: i,n
    n=size(data,1); tree%x=>data
    allocate(tree%point(max(1,n)),tree%axis(max(1,n)),tree%left(max(1,n)),tree%right(max(1,n)),ids(n))
    tree%point=0; tree%axis=0; tree%left=0; tree%right=0; tree%used=0
    do i=1,n; ids(i)=i; end do
    if(n>0) tree%root=build_node(tree,ids)
  end subroutine build_kd_tree

  recursive integer function build_node(tree,ids) result(node)
    type(kd_tree), intent(inout) :: tree
    integer, intent(in) :: ids(:)
    integer :: axis,m
    integer, allocatable :: work(:)
    if(size(ids)==0) then; node=0; return; end if
    axis=widest_axis(tree%x,ids)
    work=ids
    call sort_ids_axis(tree%x,work,axis)
    m=(size(work)+1)/2
    tree%used=tree%used+1; node=tree%used
    tree%point(node)=work(m); tree%axis(node)=axis
    if(m>1) tree%left(node)=build_node(tree,work(:m-1))
    if(m<size(work)) tree%right(node)=build_node(tree,work(m+1:))
  end function build_node

  integer function widest_axis(x,ids) result(axis)
    real(dp), intent(in) :: x(:,:)
    integer, intent(in) :: ids(:)
    integer :: j
    real(dp) :: spread,best
    axis=1; best=-1.0_dp
    do j=1,size(x,2)
      spread=maxval(x(ids,j))-minval(x(ids,j))
      if(spread>best) then; best=spread; axis=j; end if
    end do
  end function widest_axis

  subroutine sort_ids_axis(x,ids,axis)
    real(dp), intent(in) :: x(:,:)
    integer, intent(inout) :: ids(:)
    integer, intent(in) :: axis
    if(size(ids)>1) call quicksort_ids(x,ids,axis,1,size(ids))
  end subroutine sort_ids_axis

  recursive subroutine quicksort_ids(x,ids,axis,lo,hi)
    real(dp), intent(in) :: x(:,:)
    integer, intent(inout) :: ids(:)
    integer, intent(in) :: axis,lo,hi
    integer :: i,j,pivot,tmp
    if(lo>=hi) return
    pivot=ids((lo+hi)/2); i=lo; j=hi
    do
      do
        if(i>hi) exit
        if(.not.id_less(x,ids(i),pivot,axis)) exit
        i=i+1
      end do
      do
        if(j<lo) exit
        if(.not.id_less(x,pivot,ids(j),axis)) exit
        j=j-1
      end do
      if(i>j) exit
      tmp=ids(i); ids(i)=ids(j); ids(j)=tmp
      i=i+1; j=j-1
    end do
    if(lo<j) call quicksort_ids(x,ids,axis,lo,j)
    if(i<hi) call quicksort_ids(x,ids,axis,i,hi)
  end subroutine quicksort_ids

  pure logical function id_less(x,a,b,axis) result(less)
    real(dp), intent(in) :: x(:,:)
    integer, intent(in) :: a,b,axis
    if(x(a,axis)<x(b,axis)) then
      less=.true.
    else if(x(a,axis)>x(b,axis)) then
      less=.false.
    else
      less=a<b
    end if
  end function id_less

  subroutine query_kd(tree,q,k,idx,dist,exclude)
    type(kd_tree), intent(in) :: tree
    real(dp), intent(in) :: q(:)
    integer, intent(in) :: k,exclude
    integer, intent(out) :: idx(:)
    real(dp), intent(out) :: dist(:)
    if(size(idx)/=k .or. size(dist)/=k) error stop "query_kd: result size mismatch"
    call init_best(idx,dist)
    call search_node(tree,tree%root,q,exclude,idx,dist)
  end subroutine query_kd

  recursive subroutine search_node(tree,node,q,exclude,idx,dist)
    type(kd_tree), intent(in) :: tree
    integer, intent(in) :: node,exclude
    real(dp), intent(in) :: q(:)
    integer, intent(inout) :: idx(:)
    real(dp), intent(inout) :: dist(:)
    integer :: p,a,near,far
    real(dp) :: delta,d
    if(node==0) return
    p=tree%point(node); a=tree%axis(node); delta=q(a)-tree%x(p,a)
    if(delta<=0.0_dp) then; near=tree%left(node); far=tree%right(node)
    else; near=tree%right(node); far=tree%left(node); end if
    call search_node(tree,near,q,exclude,idx,dist)
    if(p/=exclude) then
      d=sqrt(max(0.0_dp,sum((tree%x(p,:)-q)**2)))
      call best_insert(p,d,idx,dist)
    end if
    if(abs(delta)<=dist(size(dist))) call search_node(tree,far,q,exclude,idx,dist)
  end subroutine search_node

  subroutine init_best(idx,dist)
    integer, intent(out) :: idx(:)
    real(dp), intent(out) :: dist(:)
    idx=-1; dist=huge(1.0_dp)
  end subroutine init_best

  subroutine best_insert(id,d,idx,dist)
    integer, intent(in) :: id
    real(dp), intent(in) :: d
    integer, intent(inout) :: idx(:)
    real(dp), intent(inout) :: dist(:)
    integer :: pos,j,n
    n=size(dist); pos=n+1
    do j=1,n
      if(d<dist(j)) then
        pos=j; exit
      else if(d>dist(j)) then
        cycle
      else if(idx(j)<0 .or. id<idx(j)) then
        pos=j; exit
      end if
    end do
    if(pos>n) return
    do j=n,pos+1,-1
      dist(j)=dist(j-1); idx(j)=idx(j-1)
    end do
    dist(pos)=d; idx(pos)=id
  end subroutine best_insert
end module fnn_neighbors
