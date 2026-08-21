! SPDX-License-Identifier: GPL-3.0-or-later
program test_vector
    use cubature, only : dp, cubature_result, hcubature_v, cuhre_v
    implicit none
    type(cubature_result) :: r
    integer :: fails
    real(dp) :: e1, e2
    fails = 0
    e1 = (1.0_dp - cos(1.0_dp)) * sin(1.0_dp) * (exp(1.0_dp) - 1.0_dp)
    e2 = 0.3078074096213368_dp
    call hcubature_v(fvec, [0.0_dp, 0.0_dp, 0.0_dp], [1.0_dp, 1.0_dp, 1.0_dp], 2, r, 1.0e-6_dp, 1.0e-12_dp)
    if (abs(r%integral(1) - e1) > 2.0e-6_dp) fails = fails + 1
    if (abs(r%integral(2) - e2) > 5.0e-6_dp) fails = fails + 1
    call cuhre_v(fvec, [0.0_dp, 0.0_dp, 0.0_dp], [1.0_dp, 1.0_dp, 1.0_dp], 2, r, 1.0e-5_dp, 1.0e-12_dp)
    if (abs(r%integral(1) - e1) > 2.0e-5_dp) fails = fails + 1
    if (fails /= 0) error stop 1
    print *, 'test_vector: PASS'
contains
    subroutine fvec(x, v)
        real(dp), intent(in) :: x(:, :)
        real(dp), intent(out) :: v(:, :)
        integer :: k
        do k = 1, size(x, 2)
            v(1, k) = sin(x(1, k)) * cos(x(2, k)) * exp(x(3, k))
            v(2, k) = 1.0_dp / (3.75_dp - cos(acos(-1.0_dp) * x(1, k)) - &
                cos(acos(-1.0_dp) * x(2, k)) - cos(acos(-1.0_dp) * x(3, k)))
        end do
    end subroutine fvec
end program test_vector
