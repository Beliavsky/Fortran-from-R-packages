! SPDX-License-Identifier: GPL-2.0-only
program test_periods
    use elliptic, only : dp, elliptic_parameters, parameters_from_g, g_from_periods, &
        primitive_periods, is_primitive, fundamental_parallelogram, wp
    implicit none
    integer::fails,mn(2)
    complex(dp)::g(2),gg(2),p(2),pp(2),z,zr
    type(elliptic_parameters)::par
    fails=0
    g=[(2.0_dp,0.3_dp),(1.0_dp,-0.99_dp)];par=parameters_from_g(g)
    gg=g_from_periods(par%omega)
    call chk(gg(1),g(1),2e-10_dp);call chk(gg(2),g(2),2e-10_dp)
    p=[(3.0_dp,1.0_dp),(1.0_dp,4.0_dp)];pp=primitive_periods(p)
    if(.not.is_primitive(pp))then;fails=fails+1;print *,'FAIL primitive';end if
    z=(13.2_dp,-9.1_dp);zr=fundamental_parallelogram(z,2.0_dp*par%omega,mn)
    call chk(wp(z,params=par),wp(zr,params=par,use_fpp=.false.),2e-10_dp)
    if(fails>0)error stop 1
    print *, 'test_periods: PASS'
contains
    subroutine chk(a,b,tol)
        complex(dp),intent(in)::a,b;real(dp),intent(in)::tol
        if(abs(a-b)>tol)then;fails=fails+1;print *,'FAIL',a,b,abs(a-b),tol;end if
    end subroutine
end program
