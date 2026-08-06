! SPDX-License-Identifier: Artistic-2.0
! Derived from the R matlab package; see COPYRIGHTS and upstream/.
program test_number_theory
    use, intrinsic :: iso_fortran_env, only : int64
    use matlab, only : primes, isprime, factors
    use test_support
    implicit none

    integer, allocatable :: p(:)
    integer(int64), allocatable :: f(:)

    p = primes(30)
    call assert_true(all(p == [2, 3, 5, 7, 11, 13, 17, 19, 23, 29]), 'primes')
    call assert_true(isprime(999983_int64), 'large prime')
    call assert_true(.not. isprime(999981_int64), 'large composite')
    f = factors(4294967295_int64)
    call assert_true(all(f == [3_int64, 5_int64, 17_int64, 257_int64, 65537_int64]), &
                     'factors 2^32-1')

    write(*, '(a)') 'test_number_theory: PASS'
end program test_number_theory
