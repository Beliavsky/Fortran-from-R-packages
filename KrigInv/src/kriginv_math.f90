module kriginv_math
  use kriginv_kinds, only : dp
  implicit none
  private
  real(dp), parameter :: pi = acos(-1.0_dp)
  public :: normal_pdf, normal_cdf, normal_quantile, bvn_cdf
contains
  elemental real(dp) function normal_pdf(x) result(v)
    real(dp), intent(in) :: x
    v = exp(-0.5_dp*x*x)/sqrt(2.0_dp*pi)
  end function normal_pdf

  elemental real(dp) function normal_cdf(x) result(v)
    real(dp), intent(in) :: x
    v = 0.5_dp*erfc(-x/sqrt(2.0_dp))
  end function normal_cdf

  elemental real(dp) function normal_quantile(p) result(x)
    real(dp), intent(in) :: p
    real(dp), parameter :: a(6) = [ &
      -3.969683028665376d1, 2.209460984245205d2, -2.759285104469687d2, &
       1.383577518672690d2, -3.066479806614716d1, 2.506628277459239d0 ]
    real(dp), parameter :: b(5) = [ &
      -5.447609879822406d1, 1.615858368580409d2, -1.556989798598866d2, &
       6.680131188771972d1, -1.328068155288572d1 ]
    real(dp), parameter :: c(6) = [ &
      -7.784894002430293d-3, -3.223964580411365d-1, -2.400758277161838d0, &
      -2.549732539343734d0, 4.374664141464968d0, 2.938163982698783d0 ]
    real(dp), parameter :: d(4) = [ &
       7.784695709041462d-3, 3.224671290700398d-1, 2.445134137142996d0, 3.754408661907416d0 ]
    real(dp), parameter :: plow = 0.02425_dp, phigh = 1.0_dp-plow
    real(dp) :: q, r, e, u
    if (p <= 0.0_dp) then
      x = -huge(1.0_dp)
      return
    else if (p >= 1.0_dp) then
      x = huge(1.0_dp)
      return
    end if
    if (p < plow) then
      q = sqrt(-2.0_dp*log(p))
      x = (((((c(1)*q+c(2))*q+c(3))*q+c(4))*q+c(5))*q+c(6))/ &
          ((((d(1)*q+d(2))*q+d(3))*q+d(4))*q+1.0_dp)
    else if (p > phigh) then
      q = sqrt(-2.0_dp*log(1.0_dp-p))
      x = -(((((c(1)*q+c(2))*q+c(3))*q+c(4))*q+c(5))*q+c(6))/ &
           ((((d(1)*q+d(2))*q+d(3))*q+d(4))*q+1.0_dp)
    else
      q = p-0.5_dp
      r = q*q
      x = (((((a(1)*r+a(2))*r+a(3))*r+a(4))*r+a(5))*r+a(6))*q/ &
          (((((b(1)*r+b(2))*r+b(3))*r+b(4))*r+b(5))*r+1.0_dp)
    end if
    e = normal_cdf(x)-p
    u = e/max(normal_pdf(x),tiny(1.0_dp))
    x = x-u/(1.0_dp+0.5_dp*x*u)
  end function normal_quantile

  real(dp) function bvn_cdf(a,b,rho) result(p)
    real(dp), intent(in) :: a,b,rho
    real(dp), parameter :: xn(32) = [ &
      -0.9972638618494816_dp,-0.9856115115452684_dp,-0.9647622555875064_dp,-0.9349060759377397_dp, &
      -0.8963211557660521_dp,-0.8493676137325700_dp,-0.7944837959679424_dp,-0.7321821187402897_dp, &
      -0.6630442669302152_dp,-0.5877157572407623_dp,-0.5068999089322294_dp,-0.4213512761306353_dp, &
      -0.3318686022821277_dp,-0.2392873622521371_dp,-0.1444719615827965_dp,-0.0483076656877383_dp, &
       0.0483076656877383_dp, 0.1444719615827965_dp, 0.2392873622521371_dp, 0.3318686022821277_dp, &
       0.4213512761306353_dp, 0.5068999089322294_dp, 0.5877157572407623_dp, 0.6630442669302152_dp, &
       0.7321821187402897_dp, 0.7944837959679424_dp, 0.8493676137325700_dp, 0.8963211557660521_dp, &
       0.9349060759377397_dp, 0.9647622555875064_dp, 0.9856115115452684_dp, 0.9972638618494816_dp ]
    real(dp), parameter :: wn(32) = [ &
      0.0070186100094701_dp,0.0162743947309057_dp,0.0253920653092621_dp,0.0342738629130214_dp, &
      0.0428358980222267_dp,0.0509980592623762_dp,0.0586840934785355_dp,0.0658222227763618_dp, &
      0.0723457941088485_dp,0.0781938957870703_dp,0.0833119242269468_dp,0.0876520930044038_dp, &
      0.0911738786957639_dp,0.0938443990808046_dp,0.0956387200792749_dp,0.0965400885147278_dp, &
      0.0965400885147278_dp,0.0956387200792749_dp,0.0938443990808046_dp,0.0911738786957639_dp, &
      0.0876520930044038_dp,0.0833119242269468_dp,0.0781938957870703_dp,0.0723457941088485_dp, &
      0.0658222227763618_dp,0.0586840934785355_dp,0.0509980592623762_dp,0.0428358980222267_dp, &
      0.0342738629130214_dp,0.0253920653092621_dp,0.0162743947309057_dp,0.0070186100094701_dp ]
    real(dp) :: r,lo,hi,mid,half,t,den,s
    integer :: i,k,nseg
    r=max(-1.0_dp,min(1.0_dp,rho))
    if(r>1.0_dp-1.0e-12_dp) then
      p=normal_cdf(min(a,b)); return
    else if(r<-1.0_dp+1.0e-12_dp) then
      p=max(0.0_dp,normal_cdf(a)-normal_cdf(-b)); return
    else if(abs(r)<1.0e-15_dp) then
      p=normal_cdf(a)*normal_cdf(b); return
    end if
    ! Plackett's identity integrates d Phi_2 / d rho from zero to rho.
    ! Splitting the correlation interval keeps the fixed Gauss-Legendre
    ! rule accurate even for |rho| close to one.
    s=0.0_dp; nseg=max(1,ceiling(8.0_dp*abs(r)))
    do k=1,nseg
      lo=r*real(k-1,dp)/real(nseg,dp)
      hi=r*real(k,dp)/real(nseg,dp)
      mid=0.5_dp*(lo+hi); half=0.5_dp*(hi-lo)
      do i=1,32
        t=mid+half*xn(i); den=max(1.0e-30_dp,1.0_dp-t*t)
        s=s+half*wn(i)*exp(-(a*a-2.0_dp*a*b*t+b*b)/(2.0_dp*den))/sqrt(den)
      end do
    end do
    p=normal_cdf(a)*normal_cdf(b)+s/(2.0_dp*pi)
    p=max(0.0_dp,min(min(normal_cdf(a),normal_cdf(b)),p))
  end function bvn_cdf
end module kriginv_math
