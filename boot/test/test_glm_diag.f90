program test_glm_diag
    use boot_kinds, only : dp
    use boot_glm_diag
    implicit none
    real(dp)::dev(3),pear(3),h(3),res(3),rd(3),rp(3),cook(3),s
    dev=[1.0_dp,-2.0_dp,0.5_dp]
    pear=[0.8_dp,-1.5_dp,0.4_dp]
    h=[0.1_dp,0.2_dp,0.3_dp]
    s=2.0_dp
    call glm_diagnostics(dev,pear,h,2,s,res,rd,rp,cook)
    if(abs(rd(1)-0.5_dp/sqrt(0.9_dp))>1.0e-12_dp)error stop 1
    if(any(cook<0.0_dp))error stop 2
    if(abs(gaussian_scale(8.0_dp,2)-2.0_dp)>1.0e-12_dp)error stop 3
    print '(a)', 'test_glm_diag: PASS'
end program test_glm_diag
