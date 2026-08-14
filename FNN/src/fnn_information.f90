! FNN-fortran: modern Fortran translation of computational code from FNN 1.1.4.1.
! Modified/translated 2026 by the FNN-fortran contributors.
! SPDX-License-Identifier: GPL-2.0-or-later
! See UPSTREAM.md and upstream/FNN-1.1.4.1 for original authorship and notices.
module fnn_information
  use fnn_kinds, only : dp
  use fnn_special, only : digamma_dp, log_unit_ball_volume
  use fnn_neighbors, only : get_knn, get_knnx, mean_log_knn_distance
  use fnn_types, only : knn_result
  implicit none
  private
  public :: entropy, crossentropy, kl_divergence, kl_dist, klx_divergence, klx_dist
  public :: mutinfo, mutual_information_entropy
contains
  function entropy(x,k,algorithm) result(h)
    real(dp), target, intent(in) :: x(:,:)
    integer, intent(in) :: k
    character(len=*), intent(in), optional :: algorithm
    real(dp), allocatable :: h(:),mld(:)
    integer :: j,n,p
    n=size(x,1); p=size(x,2)
    if(k<1 .or. k>=n) error stop "entropy: require 1 <= k < n"
    mld=mean_log_knn_distance(x,k,algorithm)
    allocate(h(k))
    do j=1,k
      h(j)=digamma_dp(real(n,dp))-digamma_dp(real(j,dp))+log_unit_ball_volume(p)+real(p,dp)*mld(j)
    end do
  end function entropy

  function crossentropy(x,y,k,algorithm) result(h)
    real(dp), target, intent(in) :: x(:,:),y(:,:)
    integer, intent(in) :: k
    character(len=*), intent(in), optional :: algorithm
    real(dp), allocatable :: h(:)
    type(knn_result) :: z
    integer :: j,m,p
    if(size(x,2)/=size(y,2)) error stop "crossentropy: dimensions differ"
    m=size(y,1); p=size(x,2)
    if(k<1 .or. k>m) error stop "crossentropy: require 1 <= k <= nrow(y)"
    z=get_knnx(y,x,k,algorithm)
    allocate(h(k))
    do j=1,k
      h(j)=digamma_dp(real(m,dp))-digamma_dp(real(j,dp))+log_unit_ball_volume(p) &
        +real(p,dp)*sum(log(z%distance(:,j)))/real(size(x,1),dp)
    end do
  end function crossentropy

  function kl_divergence(x,y,k,algorithm) result(v)
    real(dp), target, intent(in) :: x(:,:),y(:,:)
    integer, intent(in) :: k
    character(len=*), intent(in), optional :: algorithm
    real(dp), allocatable :: v(:)
    type(knn_result) :: xx,xy
    integer :: j,n,m,p
    if(size(x,2)/=size(y,2)) error stop "kl_divergence: dimensions differ"
    n=size(x,1); m=size(y,1); p=size(x,2)
    if(k<1 .or. k>=n .or. k>m) error stop "kl_divergence: invalid k"
    xx=get_knn(x,k,algorithm); xy=get_knnx(y,x,k,algorithm)
    allocate(v(k))
    do j=1,k
      v(j)=log(real(m,dp)/real(n,dp))+real(p,dp)*( &
        sum(log(xy%distance(:,j)))-sum(log(xx%distance(:,j))))/real(n,dp)
    end do
  end function kl_divergence

  function kl_dist(x,y,k,algorithm) result(v)
    real(dp), target, intent(in) :: x(:,:),y(:,:)
    integer, intent(in) :: k
    character(len=*), intent(in), optional :: algorithm
    real(dp), allocatable :: v(:)
    v=kl_divergence(x,y,k,algorithm)+kl_divergence(y,x,k,algorithm)
  end function kl_dist

  function klx_divergence(x,y,k,algorithm) result(v)
    real(dp), target, intent(in) :: x(:,:),y(:,:)
    integer, intent(in) :: k
    character(len=*), intent(in), optional :: algorithm
    real(dp), allocatable :: v(:)
    v=kl_divergence(x,y,k,algorithm)
  end function klx_divergence

  function klx_dist(x,y,k,algorithm) result(v)
    real(dp), target, intent(in) :: x(:,:),y(:,:)
    integer, intent(in) :: k
    character(len=*), intent(in), optional :: algorithm
    real(dp), allocatable :: v(:)
    v=kl_dist(x,y,k,algorithm)
  end function klx_dist

  real(dp) function mutinfo(x,y,k) result(mi)
    real(dp), intent(in) :: x(:,:),y(:,:)
    integer, intent(in) :: k
    integer :: n,i,j,nx,ny
    real(dp), allocatable :: joint_dist(:)
    real(dp) :: eps,dx,dy
    if(size(x,1)/=size(y,1)) error stop "mutinfo: row counts differ"
    n=size(x,1)
    if(k<1 .or. k>=n) error stop "mutinfo: require 1 <= k < n"
    allocate(joint_dist(n-1))
    mi=0.0_dp
    do i=1,n
      nx=0; ny=0; j=0
      call fill_joint_distances(x,y,i,joint_dist)
      call sort_real(joint_dist)
      eps=joint_dist(k)
      do j=1,n
        dx=maxval(abs(x(i,:)-x(j,:)))
        dy=maxval(abs(y(i,:)-y(j,:)))
        if(dx<eps) nx=nx+1
        if(dy<eps) ny=ny+1
      end do
      mi=mi+digamma_dp(real(nx,dp))+digamma_dp(real(ny,dp))
    end do
    mi=digamma_dp(real(n,dp))+digamma_dp(real(k,dp))-mi/real(n,dp)
  end function mutinfo

  function mutual_information_entropy(x,y,k,algorithm) result(mi)
    real(dp), target, intent(in) :: x(:,:),y(:,:)
    integer, intent(in) :: k
    character(len=*), intent(in), optional :: algorithm
    real(dp), allocatable :: mi(:),hx(:),hy(:),hxy(:),xy(:,:)
    if(size(x,1)/=size(y,1)) error stop "mutual_information_entropy: row counts differ"
    allocate(xy(size(x,1),size(x,2)+size(y,2)))
    xy(:,:size(x,2))=x; xy(:,size(x,2)+1:)=y
    hx=entropy(x,k,algorithm); hy=entropy(y,k,algorithm); hxy=entropy(xy,k,algorithm)
    mi=hx+hy-hxy
  end function mutual_information_entropy

  subroutine fill_joint_distances(x,y,i,d)
    real(dp), intent(in) :: x(:,:),y(:,:)
    integer, intent(in) :: i
    real(dp), intent(out) :: d(:)
    integer :: j,t
    real(dp) :: dx,dy
    t=0
    do j=1,size(x,1)
      if(j==i) cycle
      t=t+1
      dx=maxval(abs(x(i,:)-x(j,:))); dy=maxval(abs(y(i,:)-y(j,:)))
      d(t)=max(dx,dy)
    end do
  end subroutine fill_joint_distances

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
end module fnn_information
