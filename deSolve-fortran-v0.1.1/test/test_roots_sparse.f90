program test_roots_sparse
  use desolve_kinds, only : dp
  use desolve_types, only : ode_result, lsodar_result
  use desolve_roots_sparse, only : lsodes, lsodar
  implicit none
  type(ode_result) :: s
  type(lsodar_result) :: r
  real(dp)::y0(1),tt(3)
  y0=1.0_dp;tt=[0.0_dp,0.5_dp,1.0_dp]
  s=lsodes(rhs,y0,tt,rtol=1e-9_dp,atol=1e-11_dp)
  if(.not.s%ok())then;print *,s%status,s%message;error stop 'lsodes status';end if
  if(abs(s%y(1,3)-exp(-1.0_dp))>2e-7_dp)error stop 'lsodes accuracy'
  r=lsodar(rhs,root,y0,tt,1,rtol=1e-10_dp,atol=1e-12_dp)
  if(.not.r%solution%ok())then;print *,r%solution%status,r%solution%message;error stop 'lsodar status';end if
  if(r%nroots/=1)then;print *,'nroots',r%nroots;error stop 'lsodar root count';end if
  if(abs(r%root_time(1)-log(2.0_dp))>2e-8_dp)then;print *,r%root_time;error stop 'lsodar root time';end if
  print *,'test_roots_sparse: PASS'
contains
  subroutine rhs(t,y,dy)
    real(dp),intent(in)::t,y(:);real(dp),intent(out)::dy(:);dy=-y;if(t< -huge(1.0_dp))stop
  end subroutine rhs
  subroutine root(t,y,g)
    real(dp),intent(in)::t,y(:);real(dp),intent(out)::g(:);g(1)=y(1)-0.5_dp;if(t< -huge(1.0_dp))stop
  end subroutine root
end program test_roots_sparse
