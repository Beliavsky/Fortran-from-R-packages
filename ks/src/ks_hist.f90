! SPDX-License-Identifier: GPL-2.0-only
module ks_hist
  use ks_kinds, only: dp
  implicit none
  private
  public :: histde_1d, histde_2d, hist_predict_1d, hist_predict_2d
contains
  subroutine histde_1d(x, binw, xmin, xmax, adj, edges, density)
    real(dp), intent(in) :: x(:), binw, xmin, xmax
    real(dp), intent(in), optional :: adj
    real(dp), allocatable, intent(out) :: edges(:), density(:)
    real(dp) :: shift, lo, hi, bw
    integer :: nb, i, j, n
    n = size(x)
    if (n <= 0 .or. binw <= 0.0_dp .or. xmax <= xmin) error stop 'histde_1d: invalid input'
    shift = 0.0_dp; if (present(adj)) shift = adj
    bw = binw
    lo = xmin + shift*bw
    hi = xmax + shift*bw
    nb = max(1, ceiling((hi-lo)/bw))
    allocate(edges(nb+1), density(nb))
    do j=1,nb+1; edges(j)=lo+real(j-1,dp)*bw; end do
    density = 0.0_dp
    do i=1,n
      j = floor((x(i)-lo)/bw) + 1
      if (j == nb+1 .and. x(i) <= edges(nb+1)) j = nb
      if (j >= 1 .and. j <= nb) density(j)=density(j)+1.0_dp
    end do
    density = density/(real(n,dp)*bw)
  end subroutine

  subroutine histde_2d(x, binw, xmin, xmax, adj, edges1, edges2, density)
    real(dp), intent(in) :: x(:,:), binw(2), xmin(2), xmax(2)
    real(dp), intent(in), optional :: adj(2)
    real(dp), allocatable, intent(out) :: edges1(:), edges2(:), density(:,:)
    real(dp) :: a(2), lo(2), hi(2)
    integer :: nb(2), i,j1,j2,n
    n=size(x,1); if (size(x,2)/=2 .or. n<=0 .or. any(binw<=0.0_dp)) error stop 'histde_2d: invalid input'
    a=0.0_dp; if (present(adj)) a=adj
    lo=xmin+a*binw; hi=xmax+a*binw
    nb=max(1,ceiling((hi-lo)/binw))
    allocate(edges1(nb(1)+1),edges2(nb(2)+1),density(nb(1),nb(2)))
    do i=1,nb(1)+1; edges1(i)=lo(1)+real(i-1,dp)*binw(1); end do
    do i=1,nb(2)+1; edges2(i)=lo(2)+real(i-1,dp)*binw(2); end do
    density=0.0_dp
    do i=1,n
      j1=floor((x(i,1)-lo(1))/binw(1))+1
      j2=floor((x(i,2)-lo(2))/binw(2))+1
      if (j1==nb(1)+1 .and. x(i,1)<=edges1(nb(1)+1)) j1=nb(1)
      if (j2==nb(2)+1 .and. x(i,2)<=edges2(nb(2)+1)) j2=nb(2)
      if (j1>=1 .and. j1<=nb(1) .and. j2>=1 .and. j2<=nb(2)) density(j1,j2)=density(j1,j2)+1.0_dp
    end do
    density=density/(real(n,dp)*product(binw))
  end subroutine

  function hist_predict_1d(x, edges, density) result(f)
    real(dp), intent(in) :: x(:), edges(:), density(:)
    real(dp) :: f(size(x)); integer :: i,j,nb
    nb=size(density); if(size(edges)/=nb+1) error stop 'hist_predict_1d: shape'
    f=0.0_dp
    do i=1,size(x)
      if (x(i)>=edges(1) .and. x(i)<=edges(nb+1)) then
        j=floor((x(i)-edges(1))/(edges(2)-edges(1)))+1
        if(j>nb) j=nb
        if(j>=1) f(i)=density(j)
      end if
    end do
  end function

  function hist_predict_2d(x, edges1, edges2, density) result(f)
    real(dp), intent(in) :: x(:,:), edges1(:), edges2(:), density(:,:)
    real(dp) :: f(size(x,1)); integer :: i,j1,j2,n1,n2
    n1=size(density,1); n2=size(density,2); f=0.0_dp
    do i=1,size(x,1)
      if(x(i,1)>=edges1(1).and.x(i,1)<=edges1(n1+1).and.x(i,2)>=edges2(1).and.x(i,2)<=edges2(n2+1)) then
        j1=floor((x(i,1)-edges1(1))/(edges1(2)-edges1(1)))+1
        j2=floor((x(i,2)-edges2(1))/(edges2(2)-edges2(1)))+1
        j1=min(j1,n1); j2=min(j2,n2)
        if(j1>=1.and.j2>=1) f(i)=density(j1,j2)
      end if
    end do
  end function
end module ks_hist
