! SPDX-License-Identifier: GPL-3.0-or-later
program test_exact
    use nbconv, only : dp, nb_sum_exact
    implicit none
    real(dp), parameter :: ps(2) = [0.05_dp, 0.44_dp]
    real(dp), parameter :: phis(2) = [5.0_dp, 8.0_dp]
    integer, parameter :: counts(9) = [0, 1, 5, 10, 25, 50, 100, 200, 500]
    real(dp), parameter :: ref(9) = [ &
        4.39006988288000260e-10_dp, 4.05203450189823796e-09_dp, &
        5.15730992668707602e-07_dp, 1.16058239435808571e-05_dp, &
        6.54737124484992155e-04_dp, 5.39557745235322989e-03_dp, &
        9.38992395085907171e-03_dp, 1.06277555872333904e-03_dp, &
        9.59405945164802966e-09_dp]
    real(dp), allocatable :: pmf(:)
    real(dp) :: kmass

    pmf = nb_sum_exact(ps, phis, counts, n_terms=1000, tolerance=1.0e-10_dp, k_mass=kmass)
    call assert_close_vec(pmf, ref, 2.0e-12_dp, 2.0e-14_dp)
    call assert_close(kmass, 1.0_dp, 1.0e-12_dp, 1.0e-12_dp)
    print *, "test_exact: PASS"

contains

    subroutine assert_close(actual, expected, rtol, atol)
        real(dp), intent(in) :: actual, expected, rtol, atol
        if (abs(actual - expected) > atol + rtol * abs(expected)) error stop 1
    end subroutine assert_close

    subroutine assert_close_vec(actual, expected, rtol, atol)
        real(dp), intent(in) :: actual(:), expected(:), rtol, atol
        integer :: i
        if (size(actual) /= size(expected)) error stop 1
        do i = 1, size(actual)
            call assert_close(actual(i), expected(i), rtol, atol)
        end do
    end subroutine assert_close_vec

end program test_exact
