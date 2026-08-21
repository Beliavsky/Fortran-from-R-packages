! SPDX-License-Identifier: GPL-3.0-or-later
program test_cuba_phase
    use cubature, only : dp, cubature_result, cuhre, divonne, suave
    implicit none
    type(cubature_result) :: r
    real(dp), parameter :: expected(2) = [0.6646696797813771_dp, 0.3078074096213368_dp]
    integer :: fails
    fails = 0
    call cuhre(f, [0.0_dp, 0.0_dp, 0.0_dp], [1.0_dp, 1.0_dp, 1.0_dp], 2, r, 1.0e-5_dp, 0.0_dp)
    if (maxval(abs(r%integral - expected)) > 1.0e-5_dp) fails = fails + 1
    call divonne(f, [0.0_dp, 0.0_dp, 0.0_dp], [1.0_dp, 1.0_dp, 1.0_dp], 2, r, 1.0e-3_dp, 0.0_dp)
    if (maxval(abs(r%integral - expected)) > 7.0e-4_dp) fails = fails + 1
    call suave(f, [0.0_dp, 0.0_dp, 0.0_dp], [1.0_dp, 1.0_dp, 1.0_dp], 2, r, 1.0e-3_dp, 0.0_dp)
    if (maxval(abs(r%integral - expected)) > 7.0e-4_dp) fails = fails + 1
    if (fails /= 0) error stop 1
    print *, 'test_cuba_phase: PASS'
contains
    subroutine f(x, v)
        real(dp), intent(in) :: x(:)
        real(dp), intent(out) :: v(:)
        real(dp), parameter :: pi = acos(-1.0_dp)
        v(1) = sin(x(1)) * cos(x(2)) * exp(x(3))
        v(2) = 1.0_dp / (3.75_dp - cos(pi * x(1)) - cos(pi * x(2)) - cos(pi * x(3)))
    end subroutine f
end program test_cuba_phase
