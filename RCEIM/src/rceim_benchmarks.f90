! SPDX-License-Identifier: GPL-2.0-or-later
module rceim_benchmarks
    use rceim_kinds, only : dp
    implicit none
    private
    public :: test_fun_optimization, test_fun_optimization_2d
contains
    real(dp) function test_fun_optimization(x) result(f)
        real(dp), intent(in) :: x(:)
        if (size(x) < 1) error stop "test_fun_optimization: x must have at least one element"
        f = exp(-((x(1)-2.0_dp)**2)) + 0.9_dp*exp(-((x(1)+2.0_dp)**2)) &
            + 0.5_dp*sin(8.0_dp*x(1)) + 0.25_dp*cos(2.0_dp*x(1))
    end function test_fun_optimization

    real(dp) function test_fun_optimization_2d(x) result(f)
        real(dp), intent(in) :: x(:)
        if (size(x) < 2) error stop "test_fun_optimization_2d: x must have at least two elements"
        f = ((x(1)-4.0_dp)**2 + (x(2)+2.0_dp)**2)/50.0_dp &
            - ((x(1)+2.0_dp)**2 + (x(2)+4.0_dp)**2)/90.0_dp &
            - exp(-((x(1)-2.0_dp)**2)) - 0.9_dp*exp(-((x(2)+2.0_dp)**2)) &
            - 0.5_dp*sin(8.0_dp*x(1)) - 0.25_dp*cos(2.0_dp*x(2)) &
            + 0.25_dp*sin(x(1)*x(2)/2.0_dp) + 0.5_dp*cos(x(1)*x(2)/2.5_dp)
    end function test_fun_optimization_2d
end module rceim_benchmarks
