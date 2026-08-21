! Derived from locfit src/weight.c, GPL-2-or-later.
module locfit_kernels
  use locfit_kinds, only : dp
  use locfit_constants
  implicit none
  private
  public :: kernel_weight, kernel_is_compact, kernel_log_derivative
  public :: kernel_derivative_over_u, kernel_second_radial
  public :: kernel_integral_moment, kernel_convolution, kernel_convolution1
  public :: kernel_convolution4, kernel_convolution5, kernel_convolution6
  public :: kernel_wikk, kernel_taylor, weight_derivative_factor
  public :: rho_distance, observation_weight

contains

  pure real(dp) function kernel_weight(u, ker) result(w)
    real(dp), intent(in) :: u
    integer, intent(in) :: ker
    real(dp) :: a, v
    a = abs(u)
    select case (ker)
    case (wrect)
      w = merge(1.0_dp, 0.0_dp, a <= 1.0_dp)
    case (wepan)
      w = merge(1.0_dp-a*a, 0.0_dp, a <= 1.0_dp)
    case (wbisq)
      if (a > 1.0_dp) then
        w = 0.0_dp
      else
        v = 1.0_dp-a*a
        w = v*v
      end if
    case (wtcub)
      if (a > 1.0_dp) then
        w = 0.0_dp
      else
        v = 1.0_dp-a*a*a
        w = v*v*v
      end if
    case (wtrwt)
      if (a > 1.0_dp) then
        w = 0.0_dp
      else
        v = 1.0_dp-a*a
        w = v*v*v
      end if
    case (wququ)
      if (a > 1.0_dp) then
        w = 0.0_dp
      else
        v = 1.0_dp-a*a
        w = v*v*v*v
      end if
    case (wtria)
      w = merge(1.0_dp-a, 0.0_dp, a <= 1.0_dp)
    case (w6cub)
      if (a > 1.0_dp) then
        w = 0.0_dp
      else
        v = 1.0_dp-a*a*a
        v = v*v*v
        w = v*v
      end if
    case (wgaus)
      w = exp(-0.5_dp*(gfact*a)**2)
    case (wexpl)
      w = exp(-efact*a)
    case (wmacl)
      w = 1.0_dp/(a+1.0e-100_dp)**2
    case (wparm)
      w = 1.0_dp
    case default
      w = 0.0_dp
    end select
  end function kernel_weight

  pure logical function kernel_is_compact(ker) result(ans)
    integer, intent(in) :: ker
    ans = all(ker /= [wexpl, wgaus, wmacl, wparm])
  end function kernel_is_compact

  pure real(dp) function sgn(x) result(y)
    real(dp), intent(in) :: x
    if (x > 0.0_dp) then
      y = 1.0_dp
    else if (x < 0.0_dp) then
      y = -1.0_dp
    else
      y = 0.0_dp
    end if
  end function sgn

  pure real(dp) function kernel_log_derivative(u, ker) result(v)
    real(dp), intent(in) :: u
    integer, intent(in) :: ker
    real(dp), parameter :: eps = 1.0e-10_dp
    if (ker == wgaus) then
      v = -gfact*gfact*u
      return
    end if
    if (ker == wparm) then
      v = 0.0_dp
      return
    end if
    if (abs(u) >= 1.0_dp) then
      v = 0.0_dp
      return
    end if
    select case (ker)
    case (wrect)
      v = 0.0_dp
    case (wtria)
      v = -sgn(u)/(1.0_dp-abs(u)+eps)
    case (wepan)
      v = -2.0_dp*u/(1.0_dp-u*u+eps)
    case (wbisq)
      v = -4.0_dp*u/(1.0_dp-u*u+eps)
    case (wtrwt)
      v = -6.0_dp*u/(1.0_dp-u*u+eps)
    case (wtcub)
      v = -9.0_dp*sgn(u)*u*u/(1.0_dp-u*u*abs(u)+eps)
    case (wexpl)
      v = merge(-efact, efact, u > 0.0_dp)
    case default
      v = 0.0_dp
    end select
  end function kernel_log_derivative

  pure real(dp) function kernel_derivative_over_u(u, ker) result(v)
    real(dp), intent(in) :: u
    integer, intent(in) :: ker
    real(dp) :: a
    if (ker == wgaus) then
      v = -(gfact*gfact)*exp(-0.5_dp*(gfact*u)**2)
      return
    end if
    if (ker == wparm) then
      v = 0.0_dp
      return
    end if
    if (abs(u) > 1.0_dp) then
      v = 0.0_dp
      return
    end if
    select case (ker)
    case (wepan)
      v = -2.0_dp
    case (wbisq)
      v = -4.0_dp*(1.0_dp-u*u)
    case (wtcub)
      a = 1.0_dp-u*u*u
      v = -9.0_dp*a*a*u
    case (wtrwt)
      a = 1.0_dp-u*u
      v = -6.0_dp*a*a
    case default
      v = 0.0_dp
    end select
  end function kernel_derivative_over_u

  pure real(dp) function kernel_second_radial(u, ker) result(v)
    real(dp), intent(in) :: u
    integer, intent(in) :: ker
    real(dp) :: a
    if (ker == wgaus) then
      v = (u*gfact*gfact)**2*exp(-0.5_dp*(u*gfact)**2)
      return
    end if
    if (ker == wparm .or. u > 1.0_dp) then
      v = 0.0_dp
      return
    end if
    select case (ker)
    case (wbisq)
      v = 12.0_dp*u*u
    case (wtcub)
      a = 1.0_dp-u*u*u
      v = -9.0_dp*u*a*a + 54.0_dp*u**4*a
    case (wtrwt)
      v = 24.0_dp*u*u*(1.0_dp-u*u)
    case default
      v = 0.0_dp
    end select
  end function kernel_second_radial

  pure integer function factorial_int(n) result(f)
    integer, intent(in) :: n
    integer :: i
    f = 1
    do i = 2, n
      f = f*i
    end do
  end function factorial_int

  pure real(dp) function kernel_integral_moment(d, powers, ker) result(ans)
    integer, intent(in) :: d, ker
    integer, intent(in), optional :: powers(:)
    integer :: dj, k, nj
    real(dp) :: base, z
    nj = 0
    if (present(powers)) nj = size(powers)
    dj = d
    if (present(powers)) dj = dj + sum(powers)
    select case (ker)
    case (wrect)
      base = 1.0_dp/real(dj,dp)
    case (wepan)
      base = 2.0_dp/(real(dj,dp)*real(dj+2,dp))
    case (wbisq)
      base = 8.0_dp/(real(dj,dp)*real(dj+2,dp)*real(dj+4,dp))
    case (wtcub)
      base = 162.0_dp/(real(dj,dp)*real(dj+3,dp)*real(dj+6,dp)*real(dj+9,dp))
    case (wtrwt)
      base = 48.0_dp/(real(dj,dp)*real(dj+2,dp)*real(dj+4,dp)*real(dj+6,dp))
    case (wtria)
      base = 1.0_dp/(real(dj,dp)*real(dj+1,dp))
    case (wququ)
      base = 384.0_dp/(real(dj,dp)*real(dj+2,dp)*real(dj+4,dp)*real(dj+6,dp)*real(dj+8,dp))
    case (w6cub)
      base = 524880.0_dp/(real(dj,dp)*real(dj+3,dp)*real(dj+6,dp)*real(dj+9,dp)* &
        real(dj+12,dp)*real(dj+15,dp)*real(dj+18,dp))
    case (wgaus)
      ans = exp(real(d,dp)*log(s2pi/gfact))
      if (present(powers)) then
        do k = 1, nj
          select case (powers(k))
          case (2)
            ans = ans/(gfact*gfact)
          case (4)
            ans = ans*3.0_dp/(gfact**4)
          case default
            if (mod(powers(k),2) /= 0) then
              ans = 0.0_dp
              return
            end if
          end select
        end do
      end if
      return
    case (wexpl)
      base = real(factorial_int(dj-1),dp)/(efact**dj)
    case default
      ans = 0.0_dp
      return
    end select
    if (d == 1 .and. nj == 0) then
      ans = 2.0_dp*base
      return
    end if
    z = real(d-nj,dp)*logpi/2.0_dp - log_gamma(real(dj,dp)/2.0_dp)
    if (present(powers)) then
      do k = 1, nj
        z = z + log_gamma(real(powers(k)+1,dp)/2.0_dp)
      end do
    end if
    ans = 2.0_dp*base*exp(z)
  end function kernel_integral_moment

  pure subroutine kernel_taylor(x,ker,f,ncoef,status)
    ! Coefficients of W(x+z) as a polynomial in z; translation of wtaylor().
    real(dp),intent(in)::x
    integer,intent(in)::ker
    real(dp),intent(out)::f(:)
    integer,intent(out)::ncoef,status
    real(dp)::v
    f=0.0_dp;ncoef=0;status=lf_ok
    select case(ker)
    case(wrect)
      ncoef=1;f(1)=1.0_dp
    case(wepan)
      ncoef=3;f(1)=1.0_dp-x*x;f(2)=-2.0_dp*x;f(3)=-1.0_dp
    case(wbisq)
      ncoef=5;v=1.0_dp-x*x
      f(1)=v*v;f(2)=-4.0_dp*x*v;f(3)=4.0_dp-6.0_dp*v;f(4)=4.0_dp*x;f(5)=1.0_dp
    case(wtcub)
      ncoef=10
      if(x==1.0_dp)then
        f(1:10)=[0.0_dp,0.0_dp,0.0_dp,-27.0_dp,-81.0_dp,-108.0_dp,-81.0_dp,-36.0_dp,-9.0_dp,-1.0_dp]
      else if(x==0.0_dp)then
        f(1:10)=[1.0_dp,0.0_dp,0.0_dp,-3.0_dp,0.0_dp,0.0_dp,3.0_dp,0.0_dp,0.0_dp,-1.0_dp]
      else
        v=1.0_dp-x*x*x
        f(1)=v**3;f(2)=-9.0_dp*v*v*x*x;f(3)=x*v*(27.0_dp-36.0_dp*v)
        f(4)=-27.0_dp+v*(108.0_dp-84.0_dp*v);f(5)=-3.0_dp*x*x*(27.0_dp-42.0_dp*v)
        f(6)=x*(-108.0_dp+126.0_dp*v);f(7)=-81.0_dp+84.0_dp*v
        f(8)=-36.0_dp*x*x;f(9)=-9.0_dp*x;f(10)=-1.0_dp
      end if
    case(wtrwt)
      ncoef=7;v=1.0_dp-x*x
      f(1)=v**3;f(2)=-6.0_dp*x*v*v;f(3)=v*(12.0_dp-15.0_dp*v)
      f(4)=x*(20.0_dp*v-8.0_dp);f(5)=15.0_dp*v-12.0_dp;f(6)=-6.0_dp;f(7)=-1.0_dp
    case(wtria)
      ncoef=2;f(1)=1.0_dp-x;f(2)=-1.0_dp
    case(wququ)
      ncoef=9;v=1.0_dp-x*x
      f(1)=v**4;f(2)=-8.0_dp*x*v**3;f(3)=v*v*(24.0_dp-28.0_dp*v)
      f(4)=v*x*(56.0_dp*v-32.0_dp);f(5)=(70.0_dp*v-80.0_dp)*v+16.0_dp
      f(6)=x*(32.0_dp-56.0_dp*v);f(7)=24.0_dp-28.0_dp*v;f(8)=8.0_dp*x;f(9)=1.0_dp
    case(w6cub)
      ncoef=19;v=1.0_dp-x*x*x
      f(1)=v**6;f(2)=-18.0_dp*x*x*v**5;f(3)=x*v**4*(135.0_dp-153.0_dp*v)
      f(4)=v**3*(-540.0_dp+v*(1350.0_dp-816.0_dp*v))
      f(5)=x*x*v*v*(1215.0_dp-v*(4050.0_dp-v*3060.0_dp))
      f(6)=x*v*(-1458.0_dp+v*(9234.0_dp+v*(-16254.0_dp+v*8568.0_dp)))
      f(7)=729.0_dp-v*(10206.0_dp-v*(35154.0_dp-v*(44226.0_dp-v*18564.0_dp)))
      f(8)=x*x*(4374.0_dp-v*(30132.0_dp-v*(56862.0_dp-v*31824.0_dp)))
      f(9)=x*(12393.0_dp-v*(61479.0_dp-v*(92664.0_dp-v*43758.0_dp)))
      f(10)=21870.0_dp-v*(89100.0_dp-v*(115830.0_dp-v*48620.0_dp))
      f(11)=x*x*(26730.0_dp-v*(69498.0_dp-v*43758.0_dp))
      f(12)=x*(23814.0_dp-v*(55458.0_dp-v*31824.0_dp))
      f(13)=15849.0_dp-v*(34398.0_dp-v*18564.0_dp)
      f(14)=x*x*(7938.0_dp-8568.0_dp*v);f(15)=x*(2970.0_dp-3060.0_dp*v)
      f(16)=810.0_dp-816.0_dp*v;f(17)=153.0_dp*x*x;f(18)=18.0_dp*x;f(19)=1.0_dp
    case default
      status=lf_err
    end select
  end subroutine kernel_taylor

  pure real(dp) function weight_derivative_factor(u,scale,ker,kt,h,style,di) result(v)
    ! (d/d target W)/W for one predictor direction; translation of weightd().
    real(dp),intent(in)::u,scale,h,di
    integer,intent(in)::ker,kt,style
    real(dp)::arg
    if(scale<=0.0_dp .or. h<=0.0_dp)then;v=0.0_dp;return;end if
    if(style==stangl)then
      if(kt==kprod)then
        arg=2.0_dp*sin(u/(2.0_dp*scale))
        v=-kernel_log_derivative(arg,ker)*cos(u/(2.0_dp*scale))/(h*scale)
      else
        if(di==0.0_dp)then;v=0.0_dp;return;end if
        v=-kernel_log_derivative(di/h,ker)*sin(u/scale)/(h*scale*di)
      end if
    else if(style==stcpar)then
      v=0.0_dp
    else if(kt==kprod)then
      v=-kernel_log_derivative(u/(h*scale),ker)/(h*scale)
    else
      if(di==0.0_dp)then;v=0.0_dp;return;end if
      v=-kernel_log_derivative(di/h,ker)*u/(h*di*scale*scale)
    end if
  end function weight_derivative_factor

  pure real(dp) function kernel_convolution(vin, ker) result(ans)
    real(dp), intent(in) :: vin
    integer, intent(in) :: ker
    real(dp) :: v, v2
    v = abs(vin)
    select case (ker)
    case (wgaus)
      ans = sqrpi/gfact*exp(-(gfact*v)**2/4.0_dp)
    case (wrect)
      ans = merge(2.0_dp-v, 0.0_dp, v <= 2.0_dp)
    case (wepan)
      if (v > 2.0_dp) then
        ans = 0.0_dp
      else
        ans = (2.0_dp-v)*(16.0_dp+v*(8.0_dp-v*(16.0_dp-v*(2.0_dp+v))))/30.0_dp
      end if
    case (wbisq)
      if (v > 2.0_dp) then
        ans = 0.0_dp
      else
        v2 = 2.0_dp-v
        ans = v2**5*(16.0_dp+v*(40.0_dp+v*(36.0_dp+v*(10.0_dp+v))))/630.0_dp
      end if
    case default
      ans = 0.0_dp
    end select
  end function kernel_convolution

  pure real(dp) function kernel_convolution1(vin, ker) result(ans)
    real(dp), intent(in) :: vin
    integer, intent(in) :: ker
    real(dp) :: v, v2
    v = abs(vin)
    select case (ker)
    case (wgaus)
      ans = -0.5_dp*sqrpi*gfact*exp(-(gfact*v)**2/4.0_dp)
    case (wrect)
      ans = merge(1.0_dp, 0.0_dp, v <= 2.0_dp)
    case (wepan)
      ans = merge((-16.0_dp+v*(12.0_dp-v*v))/6.0_dp, 0.0_dp, v <= 2.0_dp)
    case (wbisq)
      if (v > 2.0_dp) then
        ans = 0.0_dp
      else
        v2 = 2.0_dp-v
        ans = -v2**4*(32.0_dp+v*(64.0_dp+v*(24.0_dp+3.0_dp*v)))/210.0_dp
      end if
    case default
      ans = 0.0_dp
    end select
  end function kernel_convolution1

  pure real(dp) function kernel_convolution4(v, ker) result(ans)
    real(dp), intent(in) :: v
    integer, intent(in) :: ker
    real(dp) :: gv
    if (ker /= wgaus) then
      ans = 0.0_dp
      return
    end if
    gv = gfact*v
    ans = exp(-gv*gv/4.0_dp)*gfact**3*(12.0_dp-gv*gv*(12.0_dp-gv*gv))*sqrpi/16.0_dp
  end function kernel_convolution4

  pure real(dp) function kernel_convolution5(v, ker) result(ans)
    real(dp), intent(in) :: v
    integer, intent(in) :: ker
    real(dp) :: gv
    if (ker /= wgaus) then
      ans = 0.0_dp
      return
    end if
    gv = gfact*v
    ans = -exp(-gv*gv/4.0_dp)*gfact**4*gv*(60.0_dp-gv*gv*(20.0_dp-gv*gv))*sqrpi/32.0_dp
  end function kernel_convolution5

  pure real(dp) function kernel_convolution6(v, ker) result(ans)
    real(dp), intent(in) :: v
    integer, intent(in) :: ker
    real(dp) :: gv, z
    if (ker /= wgaus) then
      ans = 0.0_dp
      return
    end if
    gv = (gfact*v)**2
    z = exp(-gv/4.0_dp)*(-120.0_dp+gv*(180.0_dp-gv*(30.0_dp-gv)))*0.02769459142_dp
    ans = z*gfact**5
  end function kernel_convolution6

  pure real(dp) function kernel_wikk(ker, degree) result(ans)
    integer, intent(in) :: ker, degree
    ans = 0.0_dp
    select case (degree)
    case (0,1)
      select case (ker)
      case (wrect); ans = 4.5_dp
      case (wepan); ans = 15.0_dp
      case (wbisq); ans = 35.0_dp
      case (wgaus); ans = 0.2820947918_dp*gfact**5
      case (wtcub); ans = 34.15211105_dp
      case (wtrwt); ans = 66.08391608_dp
      end select
    case (2,3)
      select case (ker)
      case (wrect); ans = 11025.0_dp
      case (wepan); ans = 39690.0_dp
      case (wbisq); ans = 110346.9231_dp
      case (wgaus); ans = 14527.43412_dp
      case (wtcub); ans = 126500.5904_dp
      case (wtrwt); ans = 254371.7647_dp
      end select
    end select
  end function kernel_wikk

  pure real(dp) function rho_distance(delta, scale, kt, style) result(rho)
    real(dp), intent(in) :: delta(:), scale(:)
    integer, intent(in) :: kt
    integer, intent(in), optional :: style(:)
    real(dp) :: r(size(delta))
    integer :: i, sty
    do i = 1, size(delta)
      sty = 0
      if (present(style)) sty = style(i)
      select case (sty)
      case (stangl)
        r(i) = 2.0_dp*sin(delta(i)/(2.0_dp*scale(i)))
      case (stcpar)
        r(i) = 0.0_dp
      case default
        r(i) = delta(i)/scale(i)
      end select
    end do
    if (size(delta) == 1) then
      rho = abs(r(1))
    else if (kt == kprod) then
      rho = maxval(abs(r))
    else
      rho = sqrt(sum(r*r))
    end if
  end function rho_distance

  pure real(dp) function observation_weight(x, target, scale, h, ker, kt, style) result(w)
    real(dp), intent(in) :: x(:), target(:), scale(:), h
    integer, intent(in) :: ker, kt
    integer, intent(in), optional :: style(:)
    real(dp) :: u(size(x)), di
    integer :: i, sty
    u = x-target
    do i = 1, size(x)
      sty = 0
      if (present(style)) sty = style(i)
      if (sty == stleft .and. u(i) > 0.0_dp) then
        w = 0.0_dp; return
      end if
      if (sty == strigh .and. u(i) < 0.0_dp) then
        w = 0.0_dp; return
      end if
    end do
    if (kt == kprod) then
      w = 1.0_dp
      do i = 1, size(x)
        sty = 0
        if (present(style)) sty = style(i)
        select case (sty)
        case (stangl)
          w = w*kernel_weight(2.0_dp*abs(sin(u(i)/(2.0_dp*scale(i))))/h, ker)
        case (stcpar)
          cycle
        case default
          if (h == 0.0_dp) then
            if (u(i) /= 0.0_dp) then
              w = 0.0_dp; return
            end if
          else
            w = w*kernel_weight(abs(u(i))/(h*scale(i)), ker)
          end if
        end select
        if (w == 0.0_dp) return
      end do
    else
      if (present(style)) then
        di = rho_distance(u, scale, ksph, style)
      else
        di = rho_distance(u, scale, ksph)
      end if
      if (h == 0.0_dp) then
        w = merge(1.0_dp, 0.0_dp, di == 0.0_dp)
      else
        w = kernel_weight(di/h, ker)
      end if
    end if
  end function observation_weight

end module locfit_kernels
