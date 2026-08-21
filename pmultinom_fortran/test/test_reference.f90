program test_reference
    use pmultinom_module, only : dp, pmultinom
    implicit none

    integer :: failures, n
    real(dp) :: p3(3)
    real(dp) :: lower1(3), upper1(3), lower2(3), upper2(3)

    failures = 0
    p3 = [0.2_dp, 0.3_dp, 0.5_dp]
    lower1 = [-1.0_dp, -1.0_dp, -1.0_dp]
    upper1 = [2.0_dp, 3.0_dp, 8.0_dp]
    lower2 = [0.0_dp, -1.0_dp, 1.0_dp]
    upper2 = [4.0_dp, 5.0_dp, 6.0_dp]

    do n = 0, 8
        call compare_case('upper-only', n, p3, lower1, upper1)
        call compare_case('two-sided', n, p3, lower2, upper2)
    end do

    if (failures /= 0) then
        print '(a,i0)', 'test_reference: FAIL ', failures
        error stop 1
    end if
    print '(a)', 'test_reference: PASS'

contains

    subroutine compare_case(label, n, p, lower, upper)
        character(len=*), intent(in) :: label
        integer, intent(in) :: n
        real(dp), intent(in) :: p(3), lower(3), upper(3)
        real(dp) :: got, expected

        got = pmultinom(n, p, lower=lower, upper=upper)
        expected = brute_probability3(n, p, lower, upper)
        if (abs(got - expected) > 2.0e-11_dp) then
            failures = failures + 1
            print '(a,1x,i0,2(1x,es24.16))', trim(label), n, got, expected
        end if
    end subroutine compare_case

    real(dp) function brute_probability3(n, p, lower, upper) result(ans)
        integer, intent(in) :: n
        real(dp), intent(in) :: p(3), lower(3), upper(3)
        integer :: c1, c2, c3
        real(dp) :: lp

        ans = 0.0_dp
        do c1 = 0, n
            do c2 = 0, n - c1
                c3 = n - c1 - c2
                if (.not. (lower(1) < real(c1, dp) .and. real(c1, dp) <= upper(1))) cycle
                if (.not. (lower(2) < real(c2, dp) .and. real(c2, dp) <= upper(2))) cycle
                if (.not. (lower(3) < real(c3, dp) .and. real(c3, dp) <= upper(3))) cycle

                lp = log_gamma(real(n + 1, dp)) - log_gamma(real(c1 + 1, dp)) - &
                     log_gamma(real(c2 + 1, dp)) - log_gamma(real(c3 + 1, dp))
                if (c1 > 0) lp = lp + real(c1, dp) * log(p(1))
                if (c2 > 0) lp = lp + real(c2, dp) * log(p(2))
                if (c3 > 0) lp = lp + real(c3, dp) * log(p(3))
                ans = ans + exp(lp)
            end do
        end do
    end function brute_probability3

end program test_reference
