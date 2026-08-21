! SPDX-License-Identifier: GPL-2.0-or-later
program test_ks
    use adequacy_kinds, only: dp
    use adequacy_math, only: kolmogorov_pvalue
    implicit none

    call check(kolmogorov_pvalue(0.2_dp, 10), 0.74871904_dp, 2.0e-12_dp)
    call check(kolmogorov_pvalue(0.15_dp, 20), 0.704467154944287_dp, 2.0e-12_dp)
    call check(kolmogorov_pvalue(0.3_dp, 5), 0.664_dp, 2.0e-12_dp)
    print '(a)', 'test_ks: PASS'
contains
    subroutine check(a, b, tol)
        real(dp), intent(in) :: a, b, tol
        if (abs(a-b) > tol) then
            print *, a, b
            error stop 1
        end if
    end subroutine check
end program test_ks
