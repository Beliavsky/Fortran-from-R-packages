! SPDX-License-Identifier: Artistic-2.0
! Derived from the R matlab package; see COPYRIGHTS and upstream/.
program demo_matlab
    use matlab, only : dp, linspace, magic, padarray, factors
    use, intrinsic :: iso_fortran_env, only : int64
    implicit none

    real(dp), allocatable :: x(:), m(:, :), padded(:, :)
    integer(int64), allocatable :: f(:)
    integer :: i

    x = linspace(0.0_dp, 1.0_dp, 6)
    write(*, '(a,*(f6.2,1x))') 'linspace: ', x

    m = magic(5)
    write(*, '(a)') 'magic(5):'
    do i = 1, size(m, 1)
        write(*, '(*(f6.0,1x))') m(i, :)
    end do

    padded = padarray(m(1:2, 1:2), 1, 1, 'symmetric', 'both')
    write(*, '(a,2(i0,1x))') 'padded shape: ', shape(padded)

    f = factors(4294967295_int64)
    write(*, '(a,*(i0,1x))') 'factors(2^32-1): ', f
end program demo_matlab
