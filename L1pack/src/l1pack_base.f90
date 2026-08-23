module l1pack_base
  use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
  use fastmatrix_base, only: dp, pi, normal_rand, chol_lower, inverse_matrix
  implicit none
  private

  real(dp), parameter :: sqrt2 = sqrt(2.0_dp)
  real(dp), parameter :: inv_sqrt2 = 1.0_dp / sqrt2

  public :: dp, pi, sqrt2, inv_sqrt2, normal_rand, chol_lower, inverse_matrix
  public :: normal_cdf_l1, normal_quantile_l1, gamma_rand_l1
  public :: sample_quantile_l1, median_l1, mahalanobis_one, logdet_chol
  public :: weighted_center_scatter, weighted_center, as78_median_center
  public :: bessel_k_ratio_prev, finite_real

contains

  pure logical function finite_real(x)
    real(dp), intent(in) :: x
    finite_real = abs(x) <= huge(x)
  end function finite_real

  pure real(dp) function normal_cdf_l1(x) result(p)
    real(dp), intent(in) :: x
    p = 0.5_dp * erfc(-x / sqrt2)
  end function normal_cdf_l1

  pure real(dp) function normal_quantile_l1(p) result(x)
    real(dp), intent(in) :: p
    real(dp), parameter :: a(6) = [ &
      -3.969683028665376e1_dp, 2.209460984245205e2_dp, &
      -2.759285104469687e2_dp, 1.383577518672690e2_dp, &
      -3.066479806614716e1_dp, 2.506628277459239_dp ]
    real(dp), parameter :: b(5) = [ &
      -5.447609879822406e1_dp, 1.615858368580409e2_dp, &
      -1.556989798598866e2_dp, 6.680131188771972e1_dp, &
      -1.328068155288572e1_dp ]
    real(dp), parameter :: c(6) = [ &
      -7.784894002430293e-3_dp, -3.223964580411365e-1_dp, &
      -2.400758277161838_dp, -2.549732539343734_dp, &
       4.374664141464968_dp, 2.938163982698783_dp ]
    real(dp), parameter :: d(4) = [ &
       7.784695709041462e-3_dp, 3.224671290700398e-1_dp, &
       2.445134137142996_dp, 3.754408661907416_dp ]
    real(dp), parameter :: plow = 0.02425_dp, phigh = 1.0_dp - plow
    real(dp) :: q, r, e, u

    if (p <= 0.0_dp) then
      x = -huge(1.0_dp)
      return
    else if (p >= 1.0_dp) then
      x = huge(1.0_dp)
      return
    end if

    if (p < plow) then
      q = sqrt(-2.0_dp * log(p))
      x = (((((c(1)*q+c(2))*q+c(3))*q+c(4))*q+c(5))*q+c(6)) / &
          ((((d(1)*q+d(2))*q+d(3))*q+d(4))*q+1.0_dp)
    else if (p <= phigh) then
      q = p - 0.5_dp
      r = q*q
      x = (((((a(1)*r+a(2))*r+a(3))*r+a(4))*r+a(5))*r+a(6))*q / &
          (((((b(1)*r+b(2))*r+b(3))*r+b(4))*r+b(5))*r+1.0_dp)
    else
      q = sqrt(-2.0_dp * log(1.0_dp-p))
      x = -(((((c(1)*q+c(2))*q+c(3))*q+c(4))*q+c(5))*q+c(6)) / &
           ((((d(1)*q+d(2))*q+d(3))*q+d(4))*q+1.0_dp)
    end if

    ! One Halley refinement substantially improves tail accuracy.
    e = normal_cdf_l1(x) - p
    u = e * sqrt(2.0_dp*pi) * exp(0.5_dp*x*x)
    x = x - u / (1.0_dp + 0.5_dp*x*u)
  end function normal_quantile_l1

  recursive real(dp) function gamma_rand_l1(shape, scale) result(x)
    real(dp), intent(in) :: shape, scale
    real(dp) :: d, c, z, u, v, g

    if (shape <= 0.0_dp .or. scale <= 0.0_dp) then
      x = ieee_value(0.0_dp, ieee_quiet_nan)
      return
    end if
    if (shape < 1.0_dp) then
      call random_number(u)
      g = gamma_rand_l1(shape + 1.0_dp, 1.0_dp)
      x = scale * g * u**(1.0_dp/shape)
      return
    end if

    d = shape - 1.0_dp/3.0_dp
    c = 1.0_dp / sqrt(9.0_dp*d)
    do
      z = normal_rand()
      v = 1.0_dp + c*z
      if (v <= 0.0_dp) cycle
      v = v*v*v
      call random_number(u)
      if (u < 1.0_dp - 0.0331_dp*z**4) exit
      if (log(u) < 0.5_dp*z*z + d*(1.0_dp-v+log(v))) exit
    end do
    x = scale*d*v
  end function gamma_rand_l1

  pure real(dp) function median_l1(x) result(m)
    real(dp), intent(in) :: x(:)
    real(dp), allocatable :: z(:)
    real(dp) :: tmp
    integer :: i, j, n
    n = size(x)
    allocate(z(n))
    z = x
    do i = 2, n
      tmp = z(i)
      j = i - 1
      do while (j >= 1)
        if (z(j) <= tmp) exit
        z(j+1) = z(j)
        j = j - 1
      end do
      z(j+1) = tmp
    end do
    if (mod(n,2) == 0) then
      m = 0.5_dp*(z(n/2)+z(n/2+1))
    else
      m = z((n+1)/2)
    end if
  end function median_l1

  pure real(dp) function sample_quantile_l1(x, k) result(q)
    real(dp), intent(in) :: x(:)
    integer, intent(in) :: k
    real(dp), allocatable :: z(:)
    real(dp) :: tmp
    integer :: i, j, kk, n
    n = size(x)
    kk = min(max(k,1),n)
    allocate(z(n))
    z = x
    do i = 2, n
      tmp = z(i)
      j = i - 1
      do while (j >= 1)
        if (z(j) <= tmp) exit
        z(j+1) = z(j)
        j = j - 1
      end do
      z(j+1) = tmp
    end do
    q = z(kk)
  end function sample_quantile_l1

  subroutine weighted_center(x, w, center)
    real(dp), intent(in) :: x(:,:), w(:)
    real(dp), intent(out) :: center(:)
    real(dp) :: sw
    sw = sum(w)
    if (sw <= tiny(1.0_dp)) then
      center = sum(x, dim=1) / real(size(x,1),dp)
    else
      center = matmul(transpose(x),w) / sw
    end if
  end subroutine weighted_center

  subroutine weighted_center_scatter(x, w, center, scatter)
    real(dp), intent(in) :: x(:,:), w(:)
    real(dp), intent(out) :: center(:), scatter(:,:)
    real(dp) :: z(size(x,2)), sw
    integer :: i, j, n
    n = size(x,1)
    call weighted_center(x,w,center)
    scatter = 0.0_dp
    do i = 1, n
      z = x(i,:) - center
      do j = 1, size(z)
        scatter(j,:) = scatter(j,:) + w(i)*z(j)*z
      end do
    end do
    ! This matches fastmatrix::FM_center_and_Scatter: divide by n,
    ! not by the sum of the EM weights.
    scatter = scatter / real(n,dp)
    sw = sum(abs(scatter))
    if (sw <= tiny(1.0_dp)) then
      scatter = 0.0_dp
      do j=1,size(scatter,1)
        scatter(j,j)=epsilon(1.0_dp)
      end do
    end if
  end subroutine weighted_center_scatter

  function mahalanobis_one(x, center, scatter) result(d2)
    real(dp), intent(in) :: x(:), center(:), scatter(:,:)
    real(dp) :: d2
    real(dp) :: inv(size(scatter,1),size(scatter,2)), z(size(x))
    integer :: info
    call inverse_matrix(scatter, inv, info)
    if (info /= 0) then
      d2 = huge(1.0_dp)
      return
    end if
    z = x-center
    d2 = max(0.0_dp, dot_product(z,matmul(inv,z)))
  end function mahalanobis_one

  function logdet_chol(scatter, info) result(v)
    real(dp), intent(in) :: scatter(:,:)
    integer, intent(out), optional :: info
    real(dp) :: v, l(size(scatter,1),size(scatter,2))
    integer :: i, ier
    call chol_lower(scatter,l,ier)
    if (present(info)) info=ier
    if (ier /= 0) then
      v = huge(1.0_dp)
      return
    end if
    v = 0.0_dp
    do i=1,size(scatter,1)
      v = v + log(l(i,i))
    end do
  end function logdet_chol

  subroutine as78_median_center(x, median, iterations)
    real(dp), intent(in) :: x(:,:)
    real(dp), intent(out) :: median(:)
    integer, intent(out) :: iterations
    real(dp), parameter :: lepsd=1.0e-4_dp, lepsr=1.0e-5_dp, lepsi=1.0e-6_dp
    integer, parameter :: icount=200, lcount=100
    integer :: n,p,i,j,l,lc,ii,ll
    real(dp) :: accum,comp,corner,d,dd,diam,delta,epsd,epsi,epsr
    real(dp) :: lambda,slam,u1,u2
    real(dp), allocatable :: c(:),z(:)

    n=size(x,1); p=size(x,2)
    allocate(c(p),z(p))
    ll=0; ii=1
    if(n==1) then
      median=x(1,:); iterations=0; return
    end if
    diam=0.0_dp
    do i=2,n
      do j=1,i-1
        accum=sum((x(i,:)-x(j,:))**2)
        diam=max(accum,diam)
      end do
    end do
    diam=sqrt(diam)
    if(diam<=tiny(1.0_dp)) then
      median=x(1,:); iterations=0; return
    end if
    epsr=lepsr*diam; epsi=lepsi*diam; epsd=lepsd*diam
    median=sum(x,dim=1)/real(n,dp)

    do l=1,icount
      corner=0.0_dp; c=0.0_dp; ll=l
      do i=1,n
        dd=sqrt(sum((x(i,:)-median)**2))
        if(dd<=epsd) then
          corner=corner+1.0_dp; ii=i
        else
          c=c+(x(i,:)-median)/dd
        end if
      end do
      d=sqrt(sum(c*c)); dd=d
      if(corner>0.5_dp) then
        if(d<=corner) then
          median=x(ii,:); iterations=-ll; return
        end if
        d=d-corner
      end if
      if(d<=epsr) then
        iterations=ll; return
      end if
      if(dd<=tiny(1.0_dp)) then
        iterations=ll; return
      end if
      c=c/dd
      u1=0.0_dp; u2=diam; lambda=0.0_dp
      do lc=1,lcount
        comp=0.0_dp
        lambda=0.5_dp*(u1+u2)
        slam=lambda*lambda
        z=median+lambda*c
        do i=1,n
          delta=sum((x(i,:)-z)**2)
          dd=sqrt(delta)
          if(dd<epsd) exit
          d=slam-sum((x(i,:)-median)**2)
          comp=comp-(d+delta)/dd
        end do
        if(dd<epsd) exit
        if(comp>0.0_dp) then
          u1=lambda
        else
          u2=lambda
        end if
        if((u2-u1)<=epsi) exit
      end do
      median=median+c*lambda
    end do
    iterations=ll
  end subroutine as78_median_center

  pure real(dp) function i0_approx(x) result(v)
    real(dp), intent(in) :: x
    real(dp) :: ax,y
    ax=abs(x)
    if(ax<3.75_dp) then
      y=(x/3.75_dp)**2
      v=1.0_dp+y*(3.5156229_dp+y*(3.0899424_dp+y*(1.2067492_dp+ &
        y*(0.2659732_dp+y*(0.0360768_dp+y*0.0045813_dp)))))
    else
      y=3.75_dp/ax
      v=(exp(ax)/sqrt(ax))*(0.39894228_dp+y*(0.01328592_dp+y*(0.00225319_dp+ &
        y*(-0.00157565_dp+y*(0.00916281_dp+y*(-0.02057706_dp+y*(0.02635537_dp+ &
        y*(-0.01647633_dp+y*0.00392377_dp))))))))
    end if
  end function i0_approx

  pure real(dp) function i1_approx(x) result(v)
    real(dp), intent(in) :: x
    real(dp) :: ax,y,ans
    ax=abs(x)
    if(ax<3.75_dp) then
      y=(x/3.75_dp)**2
      ans=ax*(0.5_dp+y*(0.87890594_dp+y*(0.51498869_dp+y*(0.15084934_dp+ &
        y*(0.02658733_dp+y*(0.00301532_dp+y*0.00032411_dp))))))
    else
      y=3.75_dp/ax
      ans=(exp(ax)/sqrt(ax))*(0.39894228_dp+y*(-0.03988024_dp+y*(-0.00362018_dp+ &
        y*(0.00163801_dp+y*(-0.01031555_dp+y*(0.02282967_dp+y*(-0.02895312_dp+ &
        y*(0.01787654_dp-y*0.00420059_dp))))))))
    end if
    if(x<0.0_dp) ans=-ans
    v=ans
  end function i1_approx

  pure subroutine k01_scaled(x,k0s,k1s)
    real(dp), intent(in) :: x
    real(dp), intent(out) :: k0s,k1s
    real(dp) :: y,k0,k1
    if(x<=2.0_dp) then
      y=x*x/4.0_dp
      k0=-log(x/2.0_dp)*i0_approx(x)+(-0.57721566_dp+y*(0.42278420_dp+ &
        y*(0.23069756_dp+y*(0.03488590_dp+y*(0.00262698_dp+ &
        y*(0.00010750_dp+y*0.00000740_dp))))))
      k1=log(x/2.0_dp)*i1_approx(x)+(1.0_dp/x)*(1.0_dp+y*(0.15443144_dp+ &
        y*(-0.67278579_dp+y*(-0.18156897_dp+y*(-0.01919402_dp+ &
        y*(-0.00110404_dp+y*(-0.00004686_dp)))))))
      k0s=k0*exp(x); k1s=k1*exp(x)
    else
      y=2.0_dp/x
      k0s=(1.0_dp/sqrt(x))*(1.25331414_dp+y*(-0.07832358_dp+y*(0.02189568_dp+ &
        y*(-0.01062446_dp+y*(0.00587872_dp+y*(-0.00251540_dp+y*0.00053208_dp))))))
      k1s=(1.0_dp/sqrt(x))*(1.25331414_dp+y*(0.23498619_dp+y*(-0.03655620_dp+ &
        y*(0.01504268_dp+y*(-0.00780353_dp+y*(0.00325614_dp-y*0.00068245_dp))))))
    end if
  end subroutine k01_scaled

  pure real(dp) function half_integer_k_scaled(n,x) result(v)
    integer, intent(in) :: n
    real(dp), intent(in) :: x
    real(dp) :: term,s
    integer :: k
    ! K_{n+1/2}(x) * exp(x); common sqrt(pi/(2x)) retained.
    s=1.0_dp; term=1.0_dp
    do k=1,n
      term=term*real((n+k)*(n-k+1),dp)/(real(k,dp)*2.0_dp*x)
      s=s+term
    end do
    v=sqrt(pi/(2.0_dp*x))*s
  end function half_integer_k_scaled

  pure real(dp) function bessel_k_ratio_prev(nu,x) result(r)
    real(dp), intent(in) :: nu,x
    integer :: n,k
    real(dp) :: km1,k0,k1
    if(x<=1.0e-12_dp) then
      if(abs(nu-0.5_dp)<1.0e-12_dp) then
        r=1.0_dp
      else if(nu>1.0_dp) then
        r=x/(2.0_dp*(nu-1.0_dp))
      else
        r=max(x*abs(log(max(x,1.0e-300_dp))),1.0e-300_dp)
      end if
      return
    end if
    if(abs(nu-nint(nu))<1.0e-12_dp) then
      n=nint(nu)
      call k01_scaled(x,k0,k1)
      if(n==0) then
        r=k1/k0
      else if(n==1) then
        r=k0/k1
      else
        km1=k0
        do k=1,n-1
          k0=k1
          k1=km1+2.0_dp*real(k,dp)*k0/x
          km1=k0
        end do
        r=k0/k1
      end if
    else
      n=nint(nu-0.5_dp)
      if(n<=0) then
        r=1.0_dp
      else
        r=half_integer_k_scaled(n-1,x)/half_integer_k_scaled(n,x)
      end if
    end if
  end function bessel_k_ratio_prev

end module l1pack_base
