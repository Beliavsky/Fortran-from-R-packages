program test_stiff
  use desolve_kinds, only : dp
  use desolve_types, only : ode_result, complex_ode_result
  use desolve_stiff, only : radau, daspk, zvode
  implicit none
  type(ode_result) :: sr, sd
  type(complex_ode_result) :: sz
  real(dp) :: y0(1), yp0(1), tt(3)
  complex(dp) :: z0(1)
  y0=1.0_dp;yp0=-1.0_dp;tt=[0.0_dp,0.5_dp,1.0_dp]
  sr=radau(rhs,y0,tt,rtol=1e-9_dp,atol=1e-11_dp)
  if(.not.sr%ok())then;print *,sr%status,sr%message;error stop 'radau status';end if
  if(abs(sr%y(1,3)-exp(-1.0_dp))>2e-7_dp)then;print *,sr%y(1,3);error stop 'radau accuracy';end if
  sd=daspk(resid,y0,yp0,tt,rtol=1e-9_dp,atol=1e-11_dp)
  if(.not.sd%ok())then;print *,sd%status,sd%message;error stop 'daspk status';end if
  if(abs(sd%y(1,3)-exp(-1.0_dp))>2e-6_dp)then;print *,sd%y(1,3);error stop 'daspk accuracy';end if
  z0(1)=(1.0_dp,0.0_dp)
  sz=zvode(zrhs,z0,tt,rtol=1e-9_dp,atol=1e-11_dp,mf=10)
  if(.not.sz%ok())then;print *,sz%status,sz%message;error stop 'zvode status';end if
  if(abs(sz%y(1,3)-exp(cmplx(0.0_dp,1.0_dp,dp)))>2e-7_dp)then
    print *,sz%y(1,3);error stop 'zvode accuracy'
  end if
  print *, 'test_stiff: PASS'
contains
  subroutine rhs(t,y,dy)
    real(dp),intent(in)::t,y(:);real(dp),intent(out)::dy(:)
    dy=-y;if(t < -huge(1.0_dp))stop
  end subroutine rhs
  subroutine resid(t,y,yp,r)
    real(dp),intent(in)::t,y(:),yp(:);real(dp),intent(out)::r(:)
    r=yp+y;if(t < -huge(1.0_dp))stop
  end subroutine resid
  subroutine zrhs(t,y,dy)
    real(dp),intent(in)::t;complex(dp),intent(in)::y(:);complex(dp),intent(out)::dy(:)
    dy=cmplx(0.0_dp,1.0_dp,dp)*y;if(t < -huge(1.0_dp))stop
  end subroutine zrhs
end program test_stiff
