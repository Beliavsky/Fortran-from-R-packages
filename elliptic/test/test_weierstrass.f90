! SPDX-License-Identifier: GPL-2.0-only
program test_weierstrass
    use elliptic, only : dp, elliptic_parameters, parameters_from_g, wp, &
        wp_prime, weierstrass_sigma, e18_10_9
    implicit none
    integer :: fails
    complex(dp) :: g(2),z,pv,pd,eq
    complex(dp) :: id(3)
    type(elliptic_parameters)::p
    fails=0
    g=[(10.0_dp,0.0_dp),(2.0_dp,0.0_dp)]
    call chk(wp((0.07_dp,0.1_dp),g=g),(-22.97450010_dp,-63.0532328_dp),1e-7_dp)
    g=[(2.0_dp,0.3_dp),(1.0_dp,-0.99_dp)];p=parameters_from_g(g)
    call chk(weierstrass_sigma((1.0_dp,0.4_dp),params=p), &
        (1.006555817_dp,0.3865197102_dp),2e-9_dp)
    g=[(1.0_dp,-0.4_dp),(2.1_dp,-0.7_dp)];p=parameters_from_g(g)
    call chk(weierstrass_sigma((10.0_dp,-8.0_dp),params=p), &
        (-1.033893831e18_dp,6.898810975e17_dp),2e11_dp)
    g=[(2.0_dp,0.0_dp),(3.0_dp,0.0_dp)];p=parameters_from_g(g)
    call chk(weierstrass_sigma((4.0_dp,0.0_dp),params=p),(-80.74922381_dp,0.0_dp),5e-7_dp)
    g=[(1.0_dp,1.0_dp),(2.0_dp,-0.33_dp)];p=parameters_from_g(g)
    call chk(wp((1.0_dp,0.3_dp),params=p),(0.8231651984_dp,-0.3567903513_dp),2e-10_dp)
    g=[(0.3123_dp,10.0_dp),(0.1_dp,-0.2222_dp)];p=parameters_from_g(g)
    call chk(wp((-4.0_dp,-4.0_dp),params=p),(-1.118985985_dp,-1.038221043_dp),5e-9_dp)
    g=[(1.44_dp,0.1_dp),(-0.3_dp,0.99_dp)];p=parameters_from_g(g);z=(1.7_dp,-0.8_dp)
    pv=wp(z,params=p);pd=wp_prime(z,params=p);eq=pd*pd-(4.0_dp*pv**3-g(1)*pv-g(2))
    call chk(eq,(0.0_dp,0.0_dp),2e-10_dp)
    id=e18_10_9(p)
    if(maxval(abs(id))>5e-10_dp)then;fails=fails+1;print *,'FAIL e18.10.9',id;end if
    if(fails>0)error stop 1
    print *, 'test_weierstrass: PASS'
contains
    subroutine chk(a,b,tol)
        complex(dp),intent(in)::a,b;real(dp),intent(in)::tol
        if(abs(a-b)>tol)then;fails=fails+1;print *,'FAIL',a,b,abs(a-b),tol;end if
    end subroutine
end program
