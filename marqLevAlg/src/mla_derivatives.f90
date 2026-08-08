! SPDX-License-Identifier: GPL-2.0-or-later
module mla_derivatives
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use mla_kinds, only : dp
    use mla_interfaces, only : objective_fn, gradient_fn
    implicit none
    private

    public :: numerical_derivatives, numerical_information_from_gradient
    public :: numerical_derivatives_scaled, numerical_information_from_gradient_scaled

contains

    subroutine numerical_derivatives(x, fn, f0, info, grad)
        real(dp), intent(in) :: x(:)
        procedure(objective_fn) :: fn
        real(dp), intent(out) :: f0
        real(dp), intent(out) :: info(:, :), grad(:)
        real(dp), allocatable :: xp(:), xm(:), fp(:), fm(:), h(:)
        real(dp) :: fij
        integer :: i, j, n

        n = size(x)
        allocate(xp(n), xm(n), fp(n), fm(n), h(n))
        f0 = fn(x)
        do i = 1, n
            h(i) = max(1.0e-7_dp, 1.0e-4_dp * abs(x(i)))
            xp = x
            xm = x
            xp(i) = xp(i) + h(i)
            xm(i) = xm(i) - h(i)
            fp(i) = fn(xp)
            fm(i) = fn(xm)
            grad(i) = (fp(i) - fm(i)) / (2.0_dp * h(i))
        end do

        info = 0.0_dp
        do i = 1, n
            do j = 1, i
                xp = x
                xp(i) = xp(i) + h(i)
                xp(j) = xp(j) + h(j)
                fij = fn(xp)
                info(j, i) = -(fij - fp(i) - fp(j) + f0) / (h(i) * h(j))
                info(i, j) = info(j, i)
            end do
        end do
    end subroutine numerical_derivatives

    subroutine numerical_information_from_gradient(x, gr, info)
        real(dp), intent(in) :: x(:)
        procedure(gradient_fn) :: gr
        real(dp), intent(out) :: info(:, :)
        real(dp), allocatable :: xp(:), xm(:), gp(:), gm(:), h(:)
        integer :: i, j, n

        n = size(x)
        allocate(xp(n), xm(n), gp(n), gm(n), h(n))
        do i = 1, n
            h(i) = max(1.0e-7_dp, 1.0e-4_dp * abs(x(i)))
        end do
        info = 0.0_dp
        do i = 1, n
            xp = x
            xm = x
            xp(i) = xp(i) + h(i)
            xm(i) = xm(i) - h(i)
            call gr(xp, gp)
            call gr(xm, gm)
            do j = i, n
                info(i, j) = (gm(j) - gp(j)) / (2.0_dp * h(i))
                info(j, i) = info(i, j)
            end do
        end do
    end subroutine numerical_information_from_gradient

    subroutine numerical_derivatives_scaled(x, fn, scale, f0, info, grad)
        real(dp), intent(in) :: x(:), scale
        procedure(objective_fn) :: fn
        real(dp), intent(out) :: f0
        real(dp), intent(out) :: info(:, :), grad(:)
        real(dp), allocatable :: xp(:), xm(:), fp(:), fm(:), h(:)
        real(dp) :: fij
        integer :: i, j, n

        n = size(x)
        allocate(xp(n), xm(n), fp(n), fm(n), h(n))
        f0 = scale * fn(x)
        do i = 1, n
            h(i) = max(1.0e-7_dp, 1.0e-4_dp * abs(x(i)))
            xp = x
            xm = x
            xp(i) = xp(i) + h(i)
            xm(i) = xm(i) - h(i)
            fp(i) = scale * fn(xp)
            fm(i) = scale * fn(xm)
            grad(i) = (fp(i) - fm(i)) / (2.0_dp * h(i))
        end do

        info = 0.0_dp
        do i = 1, n
            do j = 1, i
                xp = x
                xp(i) = xp(i) + h(i)
                xp(j) = xp(j) + h(j)
                fij = scale * fn(xp)
                info(j, i) = -(fij - fp(i) - fp(j) + f0) / (h(i) * h(j))
                info(i, j) = info(j, i)
            end do
        end do
    end subroutine numerical_derivatives_scaled

    subroutine numerical_information_from_gradient_scaled(x, gr, scale, info)
        real(dp), intent(in) :: x(:), scale
        procedure(gradient_fn) :: gr
        real(dp), intent(out) :: info(:, :)
        real(dp), allocatable :: xp(:), xm(:), gp(:), gm(:), h(:)
        integer :: i, j, n

        n = size(x)
        allocate(xp(n), xm(n), gp(n), gm(n), h(n))
        do i = 1, n
            h(i) = max(1.0e-7_dp, 1.0e-4_dp * abs(x(i)))
        end do
        info = 0.0_dp
        do i = 1, n
            xp = x
            xm = x
            xp(i) = xp(i) + h(i)
            xm(i) = xm(i) - h(i)
            call gr(xp, gp)
            call gr(xm, gm)
            do j = i, n
                info(i, j) = scale * (gm(j) - gp(j)) / (2.0_dp * h(i))
                info(j, i) = info(i, j)
            end do
        end do
    end subroutine numerical_information_from_gradient_scaled

end module mla_derivatives
