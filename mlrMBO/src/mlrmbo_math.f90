module mlrmbo_math
  use mlrmbo_kinds, only : dp, pi_dp
  use mlrmbo_rng, only : mbo_rng
  implicit none
  private
  public :: normal_pdf, normal_cdf, normal_quantile, lhs_design, mean_value, variance_value
contains
  elemental real(dp) function normal_pdf(x) result(v)
    real(dp), intent(in) :: x
    v=exp(-0.5_dp*x*x)/sqrt(2.0_dp*pi_dp)
  end function normal_pdf

  elemental real(dp) function normal_cdf(x) result(v)
    real(dp), intent(in) :: x
    v=0.5_dp*erfc(-x/sqrt(2.0_dp))
  end function normal_cdf

  real(dp) function normal_quantile(p) result(x)
    real(dp), intent(in) :: p
    real(dp) :: q,r
    real(dp), parameter :: a(6)=[-3.969683028665376e1_dp,2.209460984245205e2_dp, &
      -2.759285104469687e2_dp,1.383577518672690e2_dp,-3.066479806614716e1_dp,2.506628277459239_dp]
    real(dp), parameter :: b(5)=[-5.447609879822406e1_dp,1.615858368580409e2_dp, &
      -1.556989798598866e2_dp,6.680131188771972e1_dp,-1.328068155288572e1_dp]
    real(dp), parameter :: c(6)=[-7.784894002430293e-3_dp,-3.223964580411365e-1_dp, &
      -2.400758277161838_dp,-2.549732539343734_dp,4.374664141464968_dp,2.938163982698783_dp]
    real(dp), parameter :: d(4)=[7.784695709041462e-3_dp,3.224671290700398e-1_dp, &
      2.445134137142996_dp,3.754408661907416_dp]
    if (p <= 0.0_dp) then
      x=-huge(1.0_dp); return
    else if (p >= 1.0_dp) then
      x=huge(1.0_dp); return
    end if
    if (p < 0.02425_dp) then
      q=sqrt(-2.0_dp*log(p))
      x=(((((c(1)*q+c(2))*q+c(3))*q+c(4))*q+c(5))*q+c(6))/ &
        ((((d(1)*q+d(2))*q+d(3))*q+d(4))*q+1.0_dp)
    else if (p > 0.97575_dp) then
      q=sqrt(-2.0_dp*log(1.0_dp-p))
      x=-(((((c(1)*q+c(2))*q+c(3))*q+c(4))*q+c(5))*q+c(6))/ &
        ((((d(1)*q+d(2))*q+d(3))*q+d(4))*q+1.0_dp)
    else
      q=p-0.5_dp; r=q*q
      x=(((((a(1)*r+a(2))*r+a(3))*r+a(4))*r+a(5))*r+a(6))*q/ &
        (((((b(1)*r+b(2))*r+b(3))*r+b(4))*r+b(5))*r+1.0_dp)
    end if
  end function normal_quantile

  subroutine lhs_design(rng,n,lo,hi,x)
    type(mbo_rng), intent(inout) :: rng
    integer, intent(in) :: n
    real(dp), intent(in) :: lo(:),hi(:)
    real(dp), allocatable, intent(out) :: x(:,:)
    integer :: d,i,j,k,it
    integer, allocatable :: perm(:)
    real(dp) :: tmp
    d=size(lo); if(size(hi)/=d .or. n<1) error stop 'lhs_design: invalid dimensions'
    allocate(x(n,d),perm(n))
    do j=1,d
      perm=[(i,i=1,n)]
      do i=n,2,-1
        k=rng%randint(1,i)
        if(k/=i) then
          it=perm(i); perm(i)=perm(k); perm(k)=it
        end if
      end do
      do i=1,n
        tmp=(real(perm(i)-1,dp)+rng%uniform())/real(n,dp)
        x(i,j)=lo(j)+(hi(j)-lo(j))*tmp
      end do
    end do
  end subroutine lhs_design

  pure real(dp) function mean_value(x) result(v)
    real(dp), intent(in) :: x(:)
    if(size(x)==0) then; v=0.0_dp; else; v=sum(x)/real(size(x),dp); end if
  end function mean_value

  pure real(dp) function variance_value(x) result(v)
    real(dp), intent(in) :: x(:)
    real(dp) :: m
    if(size(x)<=1) then; v=0.0_dp; return; end if
    m=mean_value(x); v=sum((x-m)**2)/real(size(x)-1,dp)
  end function variance_value
end module mlrmbo_math
