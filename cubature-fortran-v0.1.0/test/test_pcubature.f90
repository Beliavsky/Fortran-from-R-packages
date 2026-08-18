! SPDX-License-Identifier: GPL-3.0-or-later
program test_pcubature
    use cubature, only : dp, cubature_result, pcubature, CUBATURE_SUCCESS
    implicit none
    type(cubature_result) :: r
    integer :: fails
    fails = 0
    call pcubature(fpoly, [0.0_dp, 0.0_dp, 0.0_dp], [1.0_dp, 1.0_dp, 1.0_dp], 1, r, 1.0e-10_dp, 1.0e-12_dp)
    if (abs(r%integral(1) - 1.0_dp) > 1.0e-12_dp) fails = fails + 1
    if (r%return_code /= CUBATURE_SUCCESS) fails = fails + 1
    call pcubature(fcos, [0.0_dp, 0.0_dp], [1.0_dp, 1.0_dp], 1, r, 1.0e-8_dp, 1.0e-12_dp)
    if (abs(r%integral(1) - sin(1.0_dp) ** 2) > 1.0e-8_dp) fails = fails + 1
    if (fails /= 0) error stop 1
    print *, 'test_pcubature: PASS'
contains
    subroutine fpoly(x, v)
        real(dp), intent(in) :: x(:)
        real(dp), intent(out) :: v(:)
        v(1) = product(2.0_dp * x)
    end subroutine fpoly
    subroutine fcos(x, v)
        real(dp), intent(in) :: x(:)
        real(dp), intent(out) :: v(:)
        v(1) = product(cos(x))
    end subroutine fcos
end program test_pcubature
