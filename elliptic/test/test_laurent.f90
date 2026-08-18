! SPDX-License-Identifier: GPL-2.0-only
program test_laurent
    use elliptic, only : dp, ck_coefficients, sigma_laurent, wp_laurent, &
        wp_prime_laurent, zeta_laurent
    implicit none
    complex(dp),allocatable::c(:)
    complex(dp)::g(2),z,p,pd,zz
    integer::fails
    fails=0
    c=ck_coefficients([(0.0_dp,0.0_dp),(1.0_dp,0.0_dp)],12)
    call chk(c(3),cmplx(1.0_dp/28.0_dp,0.0_dp,dp),1e-16_dp)
    call chk(c(6),cmplx(1.0_dp/10192.0_dp,0.0_dp,dp),1e-16_dp)
    call chk(c(9),cmplx(1.0_dp/5422144.0_dp,0.0_dp,dp),1e-18_dp)
    call chk(sigma_laurent((0.4_dp,1.3_dp),[(8.0_dp,0.0_dp),(4.0_dp,0.0_dp)],8), &
        (0.278080_dp,1.272785_dp),7e-8_dp)
    call chk(sigma_laurent((0.8_dp,0.4_dp),[(7.0_dp,0.0_dp),(6.0_dp,0.0_dp)],8), &
        (0.81465765_dp,0.38819473_dp),1e-8_dp)
    g=[(1.2_dp,0.1_dp),(0.3_dp,-0.2_dp)];z=(0.08_dp,0.03_dp)
    p=wp_laurent(z,g);pd=wp_prime_laurent(z,g);zz=zeta_laurent(z,g)
    call chk(pd*pd-(4.0_dp*p**3-g(1)*p-g(2)),(0.0_dp,0.0_dp),1e-7_dp)
    call chk(-pd/(2.0_dp*p*p),z,2e-4_dp)
    if(abs(zz-1.0_dp/z)>0.1_dp)then;fails=fails+1;print *,'FAIL zeta Laurent';end if
    if(fails>0)error stop 1
    print *, 'test_laurent: PASS'
contains
    subroutine chk(a,b,tol)
        complex(dp),intent(in)::a,b;real(dp),intent(in)::tol
        if(abs(a-b)>tol)then;fails=fails+1;print *,'FAIL',a,b,abs(a-b),tol;end if
    end subroutine
end program
