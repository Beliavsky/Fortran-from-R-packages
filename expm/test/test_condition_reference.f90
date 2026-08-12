program test_condition_reference
    use expm_module
    implicit none
    real(dp) :: a(10,10),c1,cf,c1e,cfe
    integer :: it
    a=0.0_dp
    a(1,2:10)=[3,10,7,3,4,9,5,9,6]
    a(2,3:10)=[5,4,3,0,5,6,3,6]
    a(3,4:10)=[5,7,7,3,7,5,6]
    a(4,5:10)=[3,7,6,8,2,7]
    a(5,6:10)=[9,5,2,7,6]
    a(6,7:10)=[8,5,4,6]
    a(7,8:10)=[5,5,3]
    a(8,9:10)=[3,5]
    a(9,10)=3
    call expm_cond_exact(a,c1,cf)
    if(abs(c1-137.455837652872_dp)>2.0e-9_dp) error stop "exact 1-norm condition mismatch"
    if(abs(cf-566.582631819923_dp)>2.0e-9_dp) error stop "exact Frobenius condition mismatch"
    call expm_cond_1_est(a,c1e)
    if(abs(c1e-c1)>1.0e-8_dp) error stop "1-norm estimate mismatch"
    call expm_cond_f_est(a,cfe,it,abstol=0.01_dp,reltol=1.0e-12_dp)
    if(abs(cfe-cf)>1.0e-7_dp) error stop "Frobenius estimate mismatch"
    print *, "test_condition_reference: PASS", c1, cf, c1e, cfe, it
end program test_condition_reference
