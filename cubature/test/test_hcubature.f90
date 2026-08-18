! SPDX-License-Identifier: GPL-3.0-or-later
program test_hcubature
    use cubature, only : dp, cubature_result, hcubature, CUBATURE_SUCCESS
    implicit none
    type(cubature_result) :: r
    real(dp), parameter :: expected = 0.7080734182735712_dp
    integer :: fails
    fails = 0
    call hcubature(fcos, [0.0_dp, 0.0_dp], [1.0_dp, 1.0_dp], 1, r, 1.0e-7_dp, 1.0e-12_dp)
    if (abs(r%integral(1) - expected) > 2.0e-7_dp) fails = fails + 1
    if (r%return_code /= CUBATURE_SUCCESS) fails = fails + 1
    call hcubature(fpoly, [0.0_dp, 0.0_dp, 0.0_dp], [1.0_dp, 1.0_dp, 1.0_dp], 1, r, 1.0e-10_dp, 1.0e-12_dp)
    if (abs(r%integral(1) - 1.0_dp) > 1.0e-12_dp) fails = fails + 1
    call hcubature(fwang, [-2.0_dp], [2.0_dp], 1, r, 1.0e-9_dp, 1.0e-12_dp)
    if (abs(r%integral(1) - 1.63564436296_dp) > 2.0e-10_dp) fails = fails + 1
    if (fails /= 0) error stop 1
    print *, 'test_hcubature: PASS'
contains
    subroutine fcos(x, v)
        real(dp), intent(in) :: x(:)
        real(dp), intent(out) :: v(:)
        v(1) = product(cos(x))
    end subroutine fcos
    subroutine fpoly(x, v)
        real(dp), intent(in) :: x(:)
        real(dp), intent(out) :: v(:)
        v(1) = product(2.0_dp * x)
    end subroutine fpoly
    subroutine fwang(x, v)
        real(dp), intent(in) :: x(:)
        real(dp), intent(out) :: v(:)
        v(1) = sin(4.0_dp * x(1)) * x(1) * (x(1) * (x(1) * (x(1) * x(1) - 4.0_dp) + 1.0_dp) - 1.0_dp)
    end subroutine fwang
end program test_hcubature
