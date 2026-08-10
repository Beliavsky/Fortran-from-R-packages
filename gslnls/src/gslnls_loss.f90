! SPDX-License-Identifier: LGPL-3.0-only
module gslnls_loss
  use gslnls_kinds, only : dp
  use gslnls_types, only : nls_loss, LOSS_DEFAULT, LOSS_HUBER, LOSS_BARRON
  use gslnls_types, only : LOSS_BISQUARE, LOSS_WELSH, LOSS_OPTIMAL
  use gslnls_types, only : LOSS_HAMPEL, LOSS_GGW, LOSS_LQQ
  implicit none
  private
  public :: loss_psi, loss_dpsi, robust_weights

contains

  pure real(dp) function loss_psi(x, loss) result(v)
    real(dp), intent(in) :: x
    type(nls_loss), intent(in) :: loss
    real(dp) :: alpha, c, c2, xc2, ax, ac, a, b, r, sx, u
    real(dp) :: r1, r2, r3, r4, a2, min_exp, k01, s5, s6
    select case (loss%kind)
    case (LOSS_DEFAULT)
      v = x
    case (LOSS_HUBER)
      c = loss%cc(1)
      if (x <= -c) then; v = -c
      else if (x < c) then; v = x
      else; v = c
      end if
    case (LOSS_BARRON)
      alpha = min(loss%cc(1), 2.0_dp); c2 = loss%cc(2)**2; xc2 = x*x/c2
      if (abs(alpha-2.0_dp) < sqrt(epsilon(1.0_dp))) then
        v = x/c2
      else if (abs(alpha) < sqrt(epsilon(1.0_dp))) then
        v = 2.0_dp*x/(x*x+2.0_dp*c2)
      else if (alpha > -1.0e8_dp) then
        v = x/c2 * (xc2/abs(alpha-2.0_dp)+1.0_dp)**(0.5_dp*alpha-1.0_dp)
      else
        v = x/c2 * exp(-0.5_dp*xc2)
      end if
    case (LOSS_BISQUARE)
      c = loss%cc(1)
      if (abs(x) > c) then
        v = 0.0_dp
      else
        a = x/c; u = 1.0_dp-a*a; v = x*u*u
      end if
    case (LOSS_WELSH)
      c = loss%cc(1); a = x/c
      if (abs(a) > 37.7_dp) then; v = 0.0_dp
      else; v = x*exp(-0.5_dp*a*a)
      end if
    case (LOSS_OPTIMAL)
      c = loss%cc(1); ac = x/c; ax = abs(ac)
      r1=-1.944_dp; r2=1.728_dp; r3=-0.312_dp; r4=0.016_dp
      if (ax > 3.0_dp) then
        v = 0.0_dp
      else if (ax > 2.0_dp) then
        a2 = ac*ac
        v = c*((((r4*a2+r3)*a2+r2)*a2+r1)*ac)
        if (ac > 0.0_dp) then; v=max(0.0_dp,v)
        else; v=-abs(v)
        end if
      else
        v = x
      end if
    case (LOSS_HAMPEL)
      c=loss%cc(1); a=1.5_dp*c; b=3.5_dp*c; r=8.0_dp*c
      if (x < 0.0_dp) then; sx=-1.0_dp; u=-x
      else; sx=1.0_dp; u=x
      end if
      if (u <= a) then; v=x
      else if (u <= b) then; v=sx*a
      else if (u <= r) then; v=sx*a*(r-u)/(r-b)
      else; v=0.0_dp
      end if
    case (LOSS_GGW)
      ax=abs(x); a=loss%cc(1); b=loss%cc(2); c=loss%cc(3)
      if (ax < c) then
        v=x
      else
        min_exp=-708.4_dp; u=-(ax-c)**b/(2.0_dp*a)
        if (u < min_exp) then; v=0.0_dp
        else; v=x*exp(u)
        end if
      end if
    case (LOSS_LQQ)
      ax=abs(x); b=loss%cc(1); c=loss%cc(2); r=loss%cc(3)
      if (ax <= c) then
        v=x
      else
        k01=b+c
        if (ax <= k01) then
          sx=merge(1.0_dp,-1.0_dp,x>0.0_dp)
          v=sx*(ax-r*(ax-c)**2/(2.0_dp*b))
        else
          s5=r-1.0_dp; s6=-2.0_dp*k01+b*r
          if (ax < k01-s6/s5) then
            sx=merge(1.0_dp,-1.0_dp,x>0.0_dp)
            v=sx*(-s6/2.0_dp-s5*s5/s6*((ax-k01)**2/2.0_dp+s6/s5*(ax-k01)))
          else
            v=0.0_dp
          end if
        end if
      end if
    case default
      v=x
    end select
  end function loss_psi

  pure real(dp) function loss_dpsi(x, loss) result(v)
    real(dp), intent(in) :: x
    type(nls_loss), intent(in) :: loss
    real(dp) :: alpha,c,c2,x2,denom,ax,ac,a,b,r,r1,r2,r3,r4,s5,aa,k01,ea
    select case(loss%kind)
    case(LOSS_DEFAULT)
      v=1.0_dp
    case(LOSS_HUBER)
      c=loss%cc(1); v=merge(0.0_dp,1.0_dp,abs(x)>=c)
    case(LOSS_BARRON)
      alpha=min(loss%cc(1),2.0_dp); c2=loss%cc(2)**2; x2=x*x
      if(abs(alpha-2.0_dp)<sqrt(epsilon(1.0_dp))) then
        v=1.0_dp/c2
      else if(abs(alpha)<sqrt(epsilon(1.0_dp))) then
        v=-2.0_dp*(x2-2.0_dp*c2)/(2.0_dp*c2+x2)**2
      else if(alpha>-1.0e8_dp) then
        denom=x2-(alpha-2.0_dp)*c2
        v=(alpha-2.0_dp)*((alpha-2.0_dp)*c2-(alpha-1.0_dp)*x2) &
          *(1.0_dp-x2/((alpha-2.0_dp)*c2))**(0.5_dp*alpha)/(denom*denom)
      else
        v=exp(-x2/(2.0_dp*c2))*(c2-x2)/(c2*c2)
      end if
    case(LOSS_BISQUARE)
      c=loss%cc(1)
      if(abs(x)>c) then; v=0.0_dp
      else; ac=x/c; x2=ac*ac; v=(1.0_dp-x2)*(1.0_dp-5.0_dp*x2)
      end if
    case(LOSS_WELSH)
      c=loss%cc(1); ac=x/c
      if(abs(ac)>37.7_dp) then; v=0.0_dp
      else; v=exp(-ac*ac/2.0_dp)*(1.0_dp-ac*ac)
      end if
    case(LOSS_OPTIMAL)
      c=loss%cc(1); ac=x/c; ax=abs(ac)
      if(ax>3.0_dp) then; v=0.0_dp
      else if(ax>2.0_dp) then
        r1=-1.944_dp; r2=1.728_dp; r3=-0.312_dp; r4=0.016_dp
        ax=ax*ax; v=r1+ax*(3.0_dp*r2+ax*(5.0_dp*r3+ax*7.0_dp*r4))
      else; v=1.0_dp
      end if
    case(LOSS_HAMPEL)
      c=loss%cc(1); a=1.5_dp*c; b=3.5_dp*c; r=8.0_dp*c; ax=abs(x)
      if(ax<=a) then; v=1.0_dp
      else if(ax<=b) then; v=0.0_dp
      else if(ax<=r) then; v=a/(b-r)
      else; v=0.0_dp
      end if
    case(LOSS_GGW)
      ax=abs(x); a=2.0_dp*loss%cc(1); b=loss%cc(2); c=loss%cc(3)
      if(ax<c) then; v=1.0_dp
      else
        ea=-(ax-c)**b/a
        if(ea < -708.4_dp) then; v=0.0_dp
        else; v=exp(ea)*(1.0_dp-b/a*ax*(ax-c)**(b-1.0_dp))
        end if
      end if
    case(LOSS_LQQ)
      ax=abs(x); b=loss%cc(1); c=loss%cc(2); r=loss%cc(3)
      if(ax<=c) then; v=1.0_dp
      else
        k01=b+c
        if(ax<=k01) then; v=1.0_dp-r/b*(ax-c)
        else
          s5=1.0_dp-r; aa=(b*r-2.0_dp*k01)/s5
          if(ax<k01+aa) then; v=-s5*((ax-k01)/aa-1.0_dp)
          else; v=0.0_dp
          end if
        end if
      end if
    case default
      v=1.0_dp
    end select
  end function loss_dpsi

  subroutine robust_weights(resid, loss, sigma, weights, psi, dpsi)
    real(dp), intent(in) :: resid(:)
    type(nls_loss), intent(in) :: loss
    real(dp), intent(out) :: sigma, weights(:), psi(:), dpsi(:)
    real(dp), allocatable :: ar(:)
    real(dp) :: z, sw
    integer :: i, n
    n=size(resid); allocate(ar(n)); ar=abs(resid)
    call sort_local(ar)
    if(mod(n,2)==1) then; sigma=1.482602218505602_dp*ar((n+1)/2)
    else; sigma=1.482602218505602_dp*0.5_dp*(ar(n/2)+ar(n/2+1))
    end if
    sigma=max(sigma,sqrt(tiny(1.0_dp)))
    sw=0.0_dp
    do i=1,n
      z=resid(i)/sigma
      psi(i)=loss_psi(z,loss); dpsi(i)=loss_dpsi(z,loss)
      if(abs(z)>sqrt(epsilon(1.0_dp))) then
        weights(i)=max(psi(i)/z,epsilon(1.0_dp))
      else
        weights(i)=max(loss_dpsi(0.0_dp,loss),epsilon(1.0_dp))
      end if
      sw=sw+weights(i)
    end do
    if(sw>0.0_dp) weights=weights*real(n,dp)/sw
  contains
    subroutine sort_local(a)
      real(dp), intent(inout) :: a(:)
      real(dp) :: t
      integer :: ii,jj
      do ii=2,size(a)
        t=a(ii); jj=ii-1
        do while(jj>=1)
          if(a(jj)<=t) exit
          a(jj+1)=a(jj); jj=jj-1
        end do
        a(jj+1)=t
      end do
    end subroutine sort_local
  end subroutine robust_weights

end module gslnls_loss
