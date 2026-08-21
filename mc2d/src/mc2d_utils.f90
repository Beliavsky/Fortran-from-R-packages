! SPDX-License-Identifier: GPL-2.0-or-later
module mc2d_utils
  use mc2d_kinds, only : dp, nan_dp
  implicit none
  private
  public :: clamp01, mean_dp, sd_dp, quantile_dp, ranks_dp, correlation_dp
  public :: sort_dp, argsort_dp, is_close
contains
  elemental real(dp) function clamp01(x) result(y)
    real(dp), intent(in) :: x
    y = min(1.0_dp, max(0.0_dp, x))
  end function clamp01

  pure logical function is_close(a,b,tol) result(ok)
    real(dp), intent(in) :: a,b
    real(dp), intent(in), optional :: tol
    real(dp) :: t
    t = sqrt(epsilon(1.0_dp)); if (present(tol)) t = tol
    ok = abs(a-b) <= t*max(1.0_dp,abs(a),abs(b))
  end function is_close

  real(dp) function mean_dp(x) result(m)
    real(dp), intent(in) :: x(:)
    if (size(x)==0) then; m=nan_dp(); else; m=sum(x)/real(size(x),dp); end if
  end function mean_dp

  real(dp) function sd_dp(x) result(s)
    real(dp), intent(in) :: x(:)
    real(dp) :: m
    if (size(x)<2) then; s=nan_dp(); return; end if
    m=mean_dp(x); s=sqrt(sum((x-m)**2)/real(size(x)-1,dp))
  end function sd_dp

  subroutine sort_dp(x)
    real(dp), intent(inout) :: x(:)
    integer :: i,j
    real(dp) :: key
    do i=2,size(x)
      key=x(i); j=i-1
      do while(j>=1)
        if(x(j)<=key) exit
        x(j+1)=x(j); j=j-1
      end do
      x(j+1)=key
    end do
  end subroutine sort_dp

  subroutine argsort_dp(x, idx)
    real(dp), intent(in) :: x(:)
    integer, intent(out) :: idx(size(x))
    integer :: i,j,key
    do i=1,size(x); idx(i)=i; end do
    do i=2,size(x)
      key=idx(i); j=i-1
      do while(j>=1)
        if(x(idx(j)) <= x(key)) exit
        idx(j+1)=idx(j); j=j-1
      end do
      idx(j+1)=key
    end do
  end subroutine argsort_dp

  real(dp) function quantile_dp(x,p) result(q)
    real(dp), intent(in) :: x(:),p
    real(dp), allocatable :: y(:)
    real(dp) :: h, frac
    integer :: j,n
    n=size(x)
    if(n==0 .or. p<0.0_dp .or. p>1.0_dp) then; q=nan_dp(); return; end if
    allocate(y(n)); y=x; call sort_dp(y)
    if(n==1) then; q=y(1); return; end if
    h = 1.0_dp + real(n-1,dp)*p
    j = floor(h); frac=h-real(j,dp)
    if(j>=n) then; q=y(n); else; q=(1.0_dp-frac)*y(j)+frac*y(j+1); end if
  end function quantile_dp

  subroutine ranks_dp(x,r)
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: r(size(x))
    integer, allocatable :: idx(:)
    integer :: i,j,k,n
    n=size(x); allocate(idx(n)); call argsort_dp(x,idx)
    i=1
    do while(i<=n)
      j=i
      do while(j<n)
        if(x(idx(j+1))/=x(idx(i))) exit
        j=j+1
      end do
      do k=i,j; r(idx(k))=0.5_dp*real(i+j,dp); end do
      i=j+1
    end do
  end subroutine ranks_dp

  real(dp) function pearson_corr(x,y) result(r)
    real(dp), intent(in) :: x(:),y(:)
    real(dp) :: mx,my,dx,dy
    if(size(x)/=size(y) .or. size(x)<2) then; r=nan_dp(); return; end if
    mx=mean_dp(x); my=mean_dp(y); dx=sum((x-mx)**2); dy=sum((y-my)**2)
    if(dx<=0.0_dp .or. dy<=0.0_dp) then; r=nan_dp(); else; r=sum((x-mx)*(y-my))/sqrt(dx*dy); end if
  end function pearson_corr

  real(dp) function kendall_corr(x,y) result(r)
    real(dp), intent(in) :: x(:),y(:)
    integer :: i,j,n
    real(dp) :: c,d,tx,ty,sx,sy
    n=size(x); if(n/=size(y) .or. n<2) then; r=nan_dp(); return; end if
    c=0; d=0; tx=0; ty=0
    do i=1,n-1; do j=i+1,n
      sx=sign(1.0_dp,x(j)-x(i)); if(x(j)==x(i)) sx=0
      sy=sign(1.0_dp,y(j)-y(i)); if(y(j)==y(i)) sy=0
      if(sx==0) tx=tx+1
      if(sy==0) ty=ty+1
      if(sx*sy>0) c=c+1
      if(sx*sy<0) d=d+1
    end do; end do
    if((c+d+tx)*(c+d+ty)<=0) then; r=nan_dp(); else; r=(c-d)/sqrt((c+d+tx)*(c+d+ty)); end if
  end function kendall_corr

  real(dp) function correlation_dp(x,y,method) result(r)
    real(dp), intent(in) :: x(:),y(:)
    character(len=*), intent(in), optional :: method
    character(len=16) :: m
    real(dp), allocatable :: rx(:),ry(:)
    m='pearson'; if(present(method)) m=adjustl(method)
    select case(trim(m))
    case('spearman')
      allocate(rx(size(x)),ry(size(y))); call ranks_dp(x,rx); call ranks_dp(y,ry); r=pearson_corr(rx,ry)
    case('kendall'); r=kendall_corr(x,y)
    case default; r=pearson_corr(x,y)
    end select
  end function correlation_dp
end module mc2d_utils
