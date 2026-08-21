program test_derivative
    use bgfd
    implicit none
    integer :: failures, family, comp
    real(dp) :: params(4), x, h, fd, pdf
    logical :: complementary

    failures = 0
    x = 1.1_dp
    h = 1.0e-5_dp
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
            fd = (bgfd_cdf(family, complementary, x+h, params(:bgfd_npar(family))) - &
                  bgfd_cdf(family, complementary, x-h, params(:bgfd_npar(family))))/(2.0_dp*h)
            pdf = bgfd_pdf(family, complementary, x, params(:bgfd_npar(family)))
            if (abs(fd-pdf) > 3.0e-8_dp*max(1.0_dp,abs(pdf))) then
                print '(a,2i3,2es20.10)', 'derivative FAIL ', family, comp, fd, pdf
                failures = failures + 1
            end if
        end do
    end do

    if (failures /= 0) error stop 1
    print '(a)', 'test_derivative: PASS'
end program test_derivative
