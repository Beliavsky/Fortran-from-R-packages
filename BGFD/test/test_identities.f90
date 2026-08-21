program test_identities
    use bgfd
    implicit none
    integer :: failures, family, comp, j
    real(dp) :: params(4), p, x, cdf, pdf, sf, hazard, h
    logical :: complementary

    failures = 0
    do family = bgfd_e, bgfd_bx
        params = 1.0_dp
        select case (family)
        case (bgfd_e)
            params(1:2) = [1.2_dp, 0.7_dp]
        case (bgfd_ee)
            params(1:3) = [1.2_dp, 1.4_dp, 0.7_dp]
        case (bgfd_w)
            params(1:3) = [0.9_dp, 1.6_dp, 0.7_dp]
        case (bgfd_ew)
            params = [0.9_dp, 1.6_dp, 1.3_dp, 0.7_dp]
        case (bgfd_f)
            params(1:3) = [1.1_dp, 2.2_dp, 0.7_dp]
        case (bgfd_l)
            params(1:3) = [1.5_dp, 2.3_dp, 0.7_dp]
        case (bgfd_b)
            params = [1.1_dp, 2.2_dp, 1.4_dp, 0.7_dp]
        case (bgfd_bx)
            params(1:2) = [1.3_dp, 0.7_dp]
        end select
        do comp = 0, 1
            complementary = comp == 1
            do j = 1, 9
                p = real(j,dp)/10.0_dp
                x = bgfd_quantile(family, complementary, p, params(:bgfd_npar(family)))
                cdf = bgfd_cdf(family, complementary, x, params(:bgfd_npar(family)))
                call check_close('cdf(quantile)', cdf, p, 2.0e-11_dp, failures)
                pdf = bgfd_pdf(family, complementary, x, params(:bgfd_npar(family)))
                sf = bgfd_survival(family, complementary, x, params(:bgfd_npar(family)))
                hazard = bgfd_hazard(family, complementary, x, params(:bgfd_npar(family)))
                h = pdf/sf
                if (pdf <= 0.0_dp .or. sf <= 0.0_dp) failures = failures + 1
                call check_close('hazard', hazard, h, 5.0e-12_dp, failures)
                call check_close('cdf+survival', cdf+sf, 1.0_dp, 5.0e-13_dp, failures)
            end do
        end do
    end do

    ! Explicitly verifies the corrected Weibull inverse exponent.
    x = q_bell_w(0.37_dp, 0.9_dp, 1.6_dp, 0.7_dp)
    call check_close('Bell-W quantile', x, 0.42517192754750094_dp, 5.0e-12_dp, failures)
    x = q_bell_ew(0.37_dp, 0.9_dp, 1.6_dp, 1.3_dp, 0.7_dp)
    call check_close('Bell-EW quantile', x, 0.55401334311064776_dp, 5.0e-12_dp, failures)

    if (failures /= 0) error stop 1
    print '(a)', 'test_identities: PASS'
contains
    subroutine check_close(name, got, expected, tol, failures)
        character(len=*), intent(in) :: name
        real(dp), intent(in) :: got, expected, tol
        integer, intent(inout) :: failures
        if (abs(got-expected) > tol*max(1.0_dp,abs(expected))) then
            print '(a,2es24.15)', trim(name)//' FAIL: ', got, expected
            failures = failures + 1
        end if
    end subroutine check_close
end program test_identities
