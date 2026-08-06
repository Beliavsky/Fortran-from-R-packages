! SPDX-License-Identifier: Artistic-2.0
! Derived from the R matlab package; see COPYRIGHTS and upstream/.
program test_numeric
    use matlab, only : dp, ceil, fix, matlab_mod, rem, nextpow2, pow2, &
                       pow2_scaled, linspace, logspace, std, std_cols, sum_cols
    use test_support
    implicit none

    real(dp), allocatable :: x(:), a(:, :), s(:)

    call assert_close(ceil(-1.9_dp), -1.0_dp, 0.0_dp, 'ceil negative')
    call assert_close(fix(-1.9_dp), -1.0_dp, 0.0_dp, 'fix negative')
    call assert_close(matlab_mod(-5.0_dp, 3.0_dp), 1.0_dp, 1.0e-14_dp, 'mod sign')
    call assert_close(rem(-5.0_dp, 3.0_dp), -2.0_dp, 1.0e-14_dp, 'rem sign')
    call assert_int_equal(nextpow2(9.0_dp), 4, 'nextpow2')
    call assert_close(pow2(3.0_dp), 8.0_dp, 0.0_dp, 'pow2')
    call assert_close(pow2_scaled(1.5_dp, 2.0_dp), 6.0_dp, 0.0_dp, 'scaled pow2')

    x = linspace(1.0_dp, 5.0_dp, 5)
    call assert_all_close(x, [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp], &
                          1.0e-14_dp, 'linspace')
    x = logspace(0.0_dp, 2.0_dp, 3)
    call assert_all_close(x, [1.0_dp, 10.0_dp, 100.0_dp], 1.0e-12_dp, 'logspace')
    call assert_close(std([1.0_dp, 5.0_dp, 9.0_dp]), 4.0_dp, 1.0e-14_dp, 'std')

    allocate(a(2, 3))
    a = reshape([1.0_dp, 7.0_dp, 5.0_dp, 15.0_dp, 9.0_dp, 22.0_dp], [2, 3])
    s = std_cols(a)
    call assert_all_close(s, [4.242640687119285_dp, 7.071067811865476_dp, &
                              9.192388155425117_dp], 1.0e-12_dp, 'std columns')
    s = sum_cols(a)
    call assert_all_close(s, [8.0_dp, 20.0_dp, 31.0_dp], 1.0e-14_dp, 'sum columns')

    write(*, '(a)') 'test_numeric: PASS'
end program test_numeric
