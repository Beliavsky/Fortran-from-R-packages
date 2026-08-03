! SPDX-License-Identifier: Artistic-2.0
module ecd_timeseries
  use ecd_kinds, only : dp, ecd_ok, ecd_invalid
  use ecd_math, only : nan_dp
  implicit none
  private

  type, public :: sample_statistics
    integer :: n=0
    real(dp) :: mean=0.0_dp, variance=0.0_dp, sd=0.0_dp
    real(dp) :: skewness=0.0_dp, kurtosis=0.0_dp
    real(dp) :: minimum=0.0_dp, maximum=0.0_dp
  end type sample_statistics

  public :: difference_series, lag_series, sample_stats, lag_stats
  public :: quantilize, empirical_quantile, manage_hist_tails

contains

  subroutine difference_series(x,lag,d)
    real(dp), intent(in) :: x(:)
    integer, intent(in), optional :: lag
    real(dp), allocatable, intent(out) :: d(:)
    integer :: l
    l=1; if(present(lag))l=lag
    if(l<0 .or. l>=size(x)) then
      allocate(d(0)); return
    end if
    allocate(d(size(x)-l)); d=x(l+1:)-x(:size(x)-l)
  end subroutine difference_series

  subroutine lag_series(x,lag,current,lagged)
    real(dp), intent(in) :: x(:)
    integer, intent(in) :: lag
    real(dp), allocatable, intent(out) :: current(:),lagged(:)
    if(lag<0 .or. lag>=size(x)) then
      allocate(current(0),lagged(0)); return
    end if
    allocate(current(size(x)-lag),lagged(size(x)-lag))
    current=x(lag+1:); lagged=x(:size(x)-lag)
  end subroutine lag_series

  function sample_stats(x) result(s)
    real(dp), intent(in) :: x(:)
    type(sample_statistics) :: s
    real(dp) :: z2,z3,z4
    integer :: n
    n=size(x); s%n=n
    if(n==0) then
      s%mean=nan_dp(); s%variance=nan_dp(); s%sd=nan_dp()
      s%skewness=nan_dp(); s%kurtosis=nan_dp(); s%minimum=nan_dp(); s%maximum=nan_dp(); return
    end if
    s%mean=sum(x)/real(n,dp); s%minimum=minval(x); s%maximum=maxval(x)
    z2=sum((x-s%mean)**2); z3=sum((x-s%mean)**3); z4=sum((x-s%mean)**4)
    if(n>1) then; s%variance=z2/real(n-1,dp); else; s%variance=0.0_dp; end if
    s%sd=sqrt(max(0.0_dp,s%variance))
    if(z2>0.0_dp) then
      s%skewness=(z3/real(n,dp))/(z2/real(n,dp))**1.5_dp
      s%kurtosis=(z4/real(n,dp))/(z2/real(n,dp))**2
    else
      s%skewness=0.0_dp; s%kurtosis=0.0_dp
    end if
  end function sample_stats

  subroutine lag_stats(x,max_lag,correlation,covariance,status)
    real(dp), intent(in) :: x(:)
    integer, intent(in) :: max_lag
    real(dp), allocatable, intent(out) :: correlation(:),covariance(:)
    integer, intent(out), optional :: status
    real(dp) :: m,v
    integer :: l,n
    if(present(status))status=ecd_ok
    n=size(x)
    if(n<2 .or. max_lag<0 .or. max_lag>=n) then
      allocate(correlation(0),covariance(0)); if(present(status))status=ecd_invalid; return
    end if
    allocate(correlation(0:max_lag),covariance(0:max_lag))
    m=sum(x)/real(n,dp); v=sum((x-m)**2)/real(n,dp)
    do l=0,max_lag
      covariance(l)=dot_product(x(1:n-l)-m,x(1+l:n)-m)/real(n-l,dp)
      correlation(l)=merge(covariance(l)/v,0.0_dp,v>0.0_dp)
    end do
  end subroutine lag_stats

  function empirical_quantile(x,p) result(q)
    real(dp), intent(in) :: x(:),p
    real(dp) :: q,h,g
    real(dp), allocatable :: y(:)
    integer :: n,j
    n=size(x)
    if(n==0 .or. p<0.0_dp .or. p>1.0_dp) then; q=nan_dp(); return; end if
    allocate(y(n)); y=x; call sort_real(y)
    if(n==1) then; q=y(1); return; end if
    h=1.0_dp+real(n-1,dp)*p; j=int(floor(h)); g=h-real(j,dp)
    if(j>=n) then; q=y(n); else; q=(1.0_dp-g)*y(j)+g*y(j+1); end if
  end function empirical_quantile

  subroutine quantilize(x,n_bins,labels,breaks,status)
    real(dp), intent(in) :: x(:)
    integer, intent(in) :: n_bins
    integer, intent(out) :: labels(:)
    real(dp), allocatable, intent(out), optional :: breaks(:)
    integer, intent(out), optional :: status
    real(dp), allocatable :: b(:)
    integer :: i,j
    if(present(status))status=ecd_ok
    if(size(labels)/=size(x) .or. n_bins<1) then
      if(size(labels)>0)labels=0; if(present(status))status=ecd_invalid
      if(present(breaks))allocate(breaks(0)); return
    end if
    allocate(b(0:n_bins))
    do j=0,n_bins; b(j)=empirical_quantile(x,real(j,dp)/real(n_bins,dp)); end do
    do i=1,size(x)
      labels(i)=n_bins
      do j=1,n_bins
        if(x(i)<=b(j)) then; labels(i)=j; exit; end if
      end do
    end do
    if(present(breaks)) then; allocate(breaks(0:n_bins)); breaks=b; end if
  end subroutine quantilize

  subroutine manage_hist_tails(x,lower_p,upper_p,trimmed,lower_tail,upper_tail)
    real(dp), intent(in) :: x(:)
    real(dp), intent(in), optional :: lower_p,upper_p
    real(dp), allocatable, intent(out) :: trimmed(:),lower_tail(:),upper_tail(:)
    real(dp) :: lp,up,ql,qu
    integer :: nt,nl,nu,i,it,il,iu
    lp=0.001_dp; up=0.999_dp
    if(present(lower_p))lp=lower_p
    if(present(upper_p))up=upper_p
    ql=empirical_quantile(x,lp); qu=empirical_quantile(x,up)
    nl=count(x<ql); nu=count(x>qu); nt=size(x)-nl-nu
    allocate(trimmed(nt),lower_tail(nl),upper_tail(nu))
    it=0; il=0; iu=0
    do i=1,size(x)
      if(x(i)<ql) then; il=il+1; lower_tail(il)=x(i)
      else if(x(i)>qu) then; iu=iu+1; upper_tail(iu)=x(i)
      else; it=it+1; trimmed(it)=x(i)
      end if
    end do
  end subroutine manage_hist_tails

  subroutine sort_real(x)
    real(dp), intent(inout) :: x(:)
    integer :: i,j
    real(dp) :: key
    do i=2,size(x)
      key=x(i); j=i-1
      do while(j>=1)
        if(x(j)<=key)exit
        x(j+1)=x(j); j=j-1
      end do
      x(j+1)=key
    end do
  end subroutine sort_real

end module ecd_timeseries
