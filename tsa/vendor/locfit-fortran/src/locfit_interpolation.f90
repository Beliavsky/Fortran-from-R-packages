! Derived from locfit src/ev_interp.c, GPL-2-or-later.
module locfit_interpolation
  use locfit_kinds, only : dp
  implicit none
  private
  public :: linear_interp, hermite2, cubic_interp, cubic_interp_derivative
  public :: rectcell_interp
contains
  pure real(dp) function linear_interp(h,d,f0,f1) result(v)
    real(dp),intent(in)::h,d,f0,f1
    if(d==0.0_dp)then
      v=f0
    else
      v=((d-h)*f0+h*f1)/d
    end if
  end function linear_interp

  pure subroutine hermite2(x,z,phi)
    real(dp),intent(in)::x,z
    real(dp),intent(out)::phi(4)
    real(dp)::h
    if(z==0.0_dp)then
      phi=[1.0_dp,0.0_dp,0.0_dp,0.0_dp];return
    end if
    h=x/z
    if(h<0.0_dp)then
      phi=[1.0_dp,0.0_dp,h,0.0_dp]
    else if(h>1.0_dp)then
      phi=[0.0_dp,1.0_dp,0.0_dp,h-1.0_dp]
    else
      phi(2)=h*h*(3.0_dp-2.0_dp*h)
      phi(1)=1.0_dp-phi(2)
      phi(3)=h*(1.0_dp-h)**2
      phi(4)=h*h*(h-1.0_dp)
    end if
  end subroutine hermite2

  pure real(dp) function cubic_interp(h,f0,f1,d0,d1) result(v)
    real(dp),intent(in)::h,f0,f1,d0,d1
    real(dp)::phi(4)
    call hermite2(h,1.0_dp,phi)
    v=phi(1)*f0+phi(2)*f1+phi(3)*d0+phi(4)*d1
  end function cubic_interp

  pure real(dp) function cubic_interp_derivative(h,f0,f1,d0,d1) result(v)
    real(dp),intent(in)::h,f0,f1,d0,d1
    real(dp)::phi(4)
    phi(2)=6.0_dp*h*(1.0_dp-h);phi(1)=-phi(2)
    phi(3)=(1.0_dp-h)*(1.0_dp-3.0_dp*h)
    phi(4)=h*(3.0_dp*h-2.0_dp)
    v=phi(1)*f0+phi(2)*f1+phi(3)*d0+phi(4)*d1
  end function cubic_interp_derivative

  pure real(dp) function rectcell_interp(x,values,ll,ur,ncoef) result(v)
    real(dp),intent(in)::x(:),values(:,:),ll(:),ur(:)
    integer,intent(in)::ncoef
    real(dp)::vv(size(values,1),size(values,2)),phi(4)
    integer::d,i,j,k,tk,nvert
    d=size(x);nvert=2**d
    if(size(values,1)<nvert)then;v=0.0_dp;return;end if
    vv=values
    if(ncoef==1)then
      do i=d,1,-1
        tk=2**(i-1)
        do j=1,tk
          vv(j,1)=linear_interp(x(i)-ll(i),ur(i)-ll(i),vv(j,1),vv(j+tk,1))
        end do
      end do
      v=vv(1,1);return
    end if
    if(ncoef==d+1)then
      do i=d,1,-1
        call hermite2(x(i)-ll(i),ur(i)-ll(i),phi)
        tk=2**(i-1);phi(3:4)=phi(3:4)*(ur(i)-ll(i))
        do j=1,tk
          vv(j,1)=phi(1)*vv(j,1)+phi(2)*vv(j+tk,1)+phi(3)*vv(j,i+1)+phi(4)*vv(j+tk,i+1)
          do k=2,i+1
            vv(j,k)=phi(1)*vv(j,k)+phi(2)*vv(j+tk,k)
          end do
        end do
      end do
      v=vv(1,1);return
    end if
    do i=d,1,-1
      call hermite2(x(i)-ll(i),ur(i)-ll(i),phi)
      tk=2**(i-1);phi(3:4)=phi(3:4)*(ur(i)-ll(i))
      do j=1,tk
        do k=1,tk
          vv(j,k)=phi(1)*vv(j,k)+phi(2)*vv(j+tk,k)+phi(3)*vv(j,k+tk)+phi(4)*vv(j+tk,k+tk)
        end do
      end do
    end do
    v=vv(1,1)
  end function rectcell_interp
end module locfit_interpolation
