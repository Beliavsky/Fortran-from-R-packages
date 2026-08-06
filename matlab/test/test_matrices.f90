! SPDX-License-Identifier: Artistic-2.0
! Derived from the R matlab package; see COPYRIGHTS and upstream/.
program test_matrices
    use matlab, only : dp, hilb, vander, magic, pascal, rosser
    use test_support
    implicit none

    real(dp), allocatable :: a(:, :), rows(:), cols(:)
    integer, allocatable :: r(:, :)
    integer :: n, i

    a = hilb(3)
    call assert_close(a(3, 3), 0.2_dp, 1.0e-14_dp, 'Hilbert matrix')

    a = vander([1.0_dp, 2.0_dp, 3.0_dp])
    call assert_all_close(a(2, :), [4.0_dp, 2.0_dp, 1.0_dp], 0.0_dp, 'Vandermonde')

    do n = 3, 6
        a = magic(n)
        allocate(rows(n), cols(n))
        rows = sum(a, dim=2)
        cols = sum(a, dim=1)
        call assert_true(maxval(abs(rows - rows(1))) < 1.0e-12_dp, 'magic row sums')
        call assert_true(maxval(abs(cols - rows(1))) < 1.0e-12_dp, 'magic column sums')
        call assert_close(sum([(a(i, i), i=1,n)]), rows(1), 1.0e-12_dp, 'magic diagonal')
        deallocate(rows, cols)
    end do

    a = pascal(4)
    call assert_all_close(a(1, :), [1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp], &
                          0.0_dp, 'Pascal first row')
    call assert_all_close(a(:, 1), [1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp], &
                          0.0_dp, 'Pascal first column')

    r = rosser()
    call assert_int_equal(r(1, 1), 611, 'Rosser first value')
    call assert_int_equal(r(8, 8), 99, 'Rosser last value')

    write(*, '(a)') 'test_matrices: PASS'
end program test_matrices
