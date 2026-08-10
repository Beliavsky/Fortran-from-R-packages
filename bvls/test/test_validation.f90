program test_validation
    use bvls, only : dp, bvls_fit, bvls_result, bvls_inconsistent_bounds, &
        bvls_no_free_variables
    implicit none
    real(dp) :: a(2,2), b(2), lower(2), upper(2)
    type(bvls_result) :: fit

    a = reshape([1.0_dp,0.0_dp,0.0_dp,1.0_dp], [2,2])
    b = [4.0_dp,-2.0_dp]
    lower = [1.5_dp,-3.0_dp]
    upper = lower
    call bvls_fit(a, b, lower, upper, fit)
    if (fit%status /= bvls_no_free_variables) error stop 'validation: fixed status'
    if (maxval(abs(fit%x-lower)) > tiny(1.0_dp)) error stop 'validation: fixed x'

    lower = [0.0_dp,2.0_dp]
    upper = [1.0_dp,1.0_dp]
    call bvls_fit(a, b, lower, upper, fit)
    if (fit%status /= bvls_inconsistent_bounds) error stop 'validation: bounds'
end program test_validation
