! SPDX-License-Identifier: GPL-3.0-or-later
program test_cuba
    use cubature, only : dp, i8, cubature_result, cuhre, divonne, suave, vegas, &
        divonne_options, suave_options, vegas_options
    implicit none
    type(cubature_result) :: r
    type(divonne_options) :: dopt
    type(suave_options) :: sopt
    type(vegas_options) :: vopt
    integer :: fails
    real(dp), parameter :: expected = 1.0_dp / 6.0_dp
    fails = 0

    call cuhre(ffact, [0.0_dp, 0.0_dp, 0.0_dp], [1.0_dp, 1.0_dp, 1.0_dp], 1, r, 1.0e-6_dp, 1.0e-12_dp)
    if (abs(r%integral(1) - expected) > 1.0e-6_dp) fails = fails + 1

    dopt%max_eval = 250000_i8
    dopt%key1 = 96
    call divonne(ffact, [0.0_dp, 0.0_dp, 0.0_dp], [1.0_dp, 1.0_dp, 1.0_dp], 1, r, 2.0e-3_dp, 1.0e-8_dp, dopt)
    if (abs(r%integral(1) - expected) > 8.0e-4_dp) fails = fails + 1

    sopt%max_eval = 250000_i8
    sopt%nnew = 512
    call suave(ffact, [0.0_dp, 0.0_dp, 0.0_dp], [1.0_dp, 1.0_dp, 1.0_dp], 1, r, 2.0e-3_dp, 1.0e-8_dp, sopt)
    if (abs(r%integral(1) - expected) > 8.0e-4_dp) fails = fails + 1

    vopt%max_eval = 250000_i8
    vopt%nstart = 4000
    vopt%nincrease = 2000
    call vegas(ffact, [0.0_dp, 0.0_dp, 0.0_dp], [1.0_dp, 1.0_dp, 1.0_dp], 1, r, 4.0e-3_dp, 1.0e-8_dp, vopt)
    if (abs(r%integral(1) - expected) > 1.5e-3_dp) fails = fails + 1

    if (fails /= 0) error stop 1
    print *, 'test_cuba: PASS'
contains
    subroutine ffact(x, v)
        real(dp), intent(in) :: x(:)
        real(dp), intent(out) :: v(:)
        v(1) = 0.5_dp * x(1) * x(1)
    end subroutine ffact
end program test_cuba
