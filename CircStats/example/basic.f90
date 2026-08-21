program basic
    use circstats
    implicit none
    real(dp) :: x(6)
    type(vm_fit_result) :: fit
    type(test_result) :: ray
    x=[0.1_dp,0.35_dp,0.5_dp,0.8_dp,0.9_dp,1.1_dp]
    fit=vm_ml(x)
    ray=rayleigh_test(x)
    print '(a,f10.6)', 'mean direction = ',fit%mu
    print '(a,f10.6)', 'kappa          = ',fit%kappa
    print '(a,f10.6)', 'Rayleigh p     = ',ray%p_value
end program
