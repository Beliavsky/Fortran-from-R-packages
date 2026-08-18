program test_moments
    use bivgeo, only : dp, bivgeo_params, make_bivgeo_params, cfbivgeo, covbivgeo, corbivgeo, mean_bivgeo, variance_bivgeo
    implicit none

    type(bivgeo_params) :: theta
    real(dp) :: mu(2), var(2), covar

    theta = make_bivgeo_params(0.5_dp, 0.5_dp, 0.7_dp)
    call check(abs(cfbivgeo(theta) - 2.5174825174825175_dp) < 1.0e-13_dp, 1)
    call check(abs(covbivgeo(theta) - 0.1506186121570737_dp) < 1.0e-13_dp, 2)
    call check(abs(corbivgeo(theta) - 0.1818181818181818_dp) < 1.0e-13_dp, 3)

    mu = mean_bivgeo(theta)
    var = variance_bivgeo(theta)
    covar = covbivgeo(theta)
    call check(abs(cfbivgeo(theta) - mu(1) * mu(2) - covar) < 2.0e-13_dp, 4)
    call check(abs(covar / sqrt(var(1) * var(2)) - corbivgeo(theta)) < 2.0e-13_dp, 5)

    theta = make_bivgeo_params(0.9_dp, 0.9_dp, 0.9_dp)
    call check(abs(cfbivgeo(theta) - 35.15245678772579_dp) < 2.0e-12_dp, 6)
    call check(abs(covbivgeo(theta) - 7.451625762795031_dp) < 2.0e-12_dp, 7)
    call check(abs(corbivgeo(theta) - 0.3321033210332104_dp) < 2.0e-13_dp, 8)

    print '(a)', 'test_moments: PASS'

contains

    subroutine check(condition, code)
        logical, intent(in) :: condition
        integer, intent(in) :: code
        if (.not. condition) then
            print '(a,i0)', 'test_moments: FAIL ', code
            error stop 1
        end if
    end subroutine check

end program test_moments
