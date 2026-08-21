! SPDX-License-Identifier: GPL-2.0-only
program test_modular_arithmetic
    use elliptic, only : dp, i8, pi, dedekind_eta, eta_series, j_invariant, &
        modular_lambda, divisor_sigma, totient, mobius, liouville, g2_from_periods, &
        g2_divisor, g2_lambert, g3_from_periods, g3_divisor, g3_lambert
    implicit none
    integer::fails
    complex(dp)::tau,b(2),l0
    fails=0;tau=(0.0_dp,1.0_dp)
    call chk(dedekind_eta(tau),cmplx(gamma(0.25_dp)/(2.0_dp*pi**0.75_dp),0.0_dp,dp),2e-15_dp)
    call chk(dedekind_eta((0.2_dp,1.2_dp)),eta_series((0.2_dp,1.2_dp),500),2e-14_dp)
    call chk(j_invariant(tau),(1.0_dp,0.0_dp),3e-14_dp)
    l0=modular_lambda((0.4_dp,1.1_dp))
    call chk(modular_lambda((1.4_dp,1.1_dp)),l0/(l0-1.0_dp),2e-12_dp)
    if(divisor_sigma(140_i8)/=336_i8)then;fails=fails+1;print *,'FAIL divisor';end if
    if(totient(25_i8)/=20_i8 .or. mobius(30_i8)/=-1 .or. liouville(12_i8)/=-1)then
        fails=fails+1;print *,'FAIL arithmetic'
    end if
    b=[(1.0_dp,0.0_dp),(0.2_dp,1.3_dp)]
    call chk(g2_from_periods(b),g2_divisor(b,100),2e-12_dp)
    call chk(g2_from_periods(b),g2_lambert(b,100),2e-12_dp)
    call chk(g3_from_periods(b),g3_divisor(b,100),2e-12_dp)
    call chk(g3_from_periods(b),g3_lambert(b,100),2e-12_dp)
    if(fails>0)error stop 1
    print *, 'test_modular_arithmetic: PASS'
contains
    subroutine chk(a,bv,tol)
        complex(dp),intent(in)::a,bv;real(dp),intent(in)::tol
        if(abs(a-bv)>tol)then;fails=fails+1;print *,'FAIL',a,bv,abs(a-bv),tol;end if
    end subroutine
end program
