! SPDX-License-Identifier: GPL-2.0-only
program test_theta_jacobi
    use elliptic, only : dp, sn, cn, dn, dc, cs, sc, theta_big, &
        e16_28_1, e16_28_2, e16_28_3, e16_28_4, e16_28_5
    implicit none
    integer :: fails
    complex(dp) :: m, u, lhs
    fails=0
    call chk(dn((0.2_dp,0.0_dp),(0.19_dp,0.0_dp)),(0.996253_dp,0.0_dp),2e-6_dp)
    call chk(dn((0.2_dp,0.0_dp),(0.81_dp,0.0_dp)),(0.98406_dp,0.0_dp),2e-5_dp)
    call chk(cn((0.2_dp,0.0_dp),(0.81_dp,0.0_dp)),(0.980278_dp,0.0_dp),2e-6_dp)
    call chk(dc((0.672_dp,0.0_dp),(0.36_dp,0.0_dp)),(1.174_dp,0.0_dp),2e-4_dp)
    call chk(theta_big((0.6_dp,0.0_dp),(0.36_dp,0.0_dp)),(0.97357_dp,0.0_dp),2e-5_dp)
    call chk(cs((0.5360162_dp,0.0_dp),(0.09_dp,0.0_dp)),(1.6918083_dp,0.0_dp),2e-7_dp)
    call chk(sn((0.61802_dp,0.0_dp),(0.5_dp,0.0_dp)),(0.56458_dp,0.0_dp),2e-5_dp)
    call chk(sc((0.61802_dp,0.0_dp),(0.5_dp,0.0_dp)),(0.68402_dp,0.0_dp),2e-5_dp)
    m=(0.1_dp,0.1123312_dp);u=(1.7_dp,0.4_dp)
    lhs=sn((0.0_dp,1.0_dp)*u,m)-(0.0_dp,1.0_dp)*sc(u,1.0_dp-m)
    call chk(lhs,(0.0_dp,0.0_dp),2e-10_dp)
    call chk(cn((0.0_dp,1.0_dp)*u,m)-1.0_dp/cn(u,1.0_dp-m),(0.0_dp,0.0_dp),2e-10_dp)
    call chk(e16_28_1((1.3_dp,0.2_dp),(0.234_dp,0.1_dp)),(0.0_dp,0.0_dp),2e-13_dp)
    call chk(e16_28_2((1.3_dp,0.2_dp),(0.234_dp,0.1_dp)),(0.0_dp,0.0_dp),2e-13_dp)
    call chk(e16_28_3((1.3_dp,0.2_dp),(0.234_dp,0.1_dp)),(0.0_dp,0.0_dp),2e-13_dp)
    call chk(e16_28_4((1.3_dp,0.2_dp),(0.234_dp,0.1_dp)),(0.0_dp,0.0_dp),2e-13_dp)
    call chk(e16_28_5((0.234_dp,0.1_dp)),(0.0_dp,0.0_dp),2e-13_dp)
    if(fails>0)error stop 1
    print *, 'test_theta_jacobi: PASS'
contains
    subroutine chk(a,b,tol)
        complex(dp),intent(in)::a,b
        real(dp),intent(in)::tol
        if(abs(a-b)>tol)then
            fails=fails+1
            print *, 'FAIL',a,b,abs(a-b),tol
        end if
    end subroutine
end program
