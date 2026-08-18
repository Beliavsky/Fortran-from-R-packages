! SPDX-License-Identifier: GPL-3.0-or-later
program test_infinite
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_positive_inf, ieee_negative_inf
    use cubature, only : dp, cubature_result, cubintegrate
    implicit none
    type(cubature_result) :: r
    real(dp) :: lo(1), hi(1)
    lo(1) = ieee_value(0.0_dp, ieee_negative_inf)
    hi(1) = ieee_value(0.0_dp, ieee_positive_inf)
    call cubintegrate(fgauss, lo, hi, 1, 'hcubature', r, 1.0e-7_dp, 1.0e-10_dp)
    if (abs(r%integral(1) - 1.0_dp) > 2.0e-7_dp) error stop 1
    print *, 'test_infinite: PASS'
contains
    subroutine fgauss(x, v)
        real(dp), intent(in) :: x(:)
        real(dp), intent(out) :: v(:)
        v(1) = exp(-x(1) * x(1)) / sqrt(acos(-1.0_dp))
    end subroutine fgauss
end program test_infinite
