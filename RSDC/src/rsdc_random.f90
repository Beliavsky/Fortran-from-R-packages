! SPDX-License-Identifier: GPL-3.0-only
module rsdc_random
    use rsdc_kinds, only: dp, rsdc_pi
    use rsdc_linalg, only: cholesky_lower
    implicit none
    private
    public :: seed_rsdc, random_normal, random_multivariate_normal, categorical_draw

contains

    subroutine seed_rsdc(seed)
        integer, intent(in) :: seed
        integer :: n, i
        integer, allocatable :: put(:)
        call random_seed(size=n)
        allocate(put(n))
        do i = 1, n
            put(i) = modulo(seed + 104729 * i, huge(1) - 1)
            if (put(i) <= 0) put(i) = i
        end do
        call random_seed(put=put)
    end subroutine seed_rsdc

    real(dp) function random_normal() result(z)
        real(dp) :: u1, u2
        call random_number(u1)
        call random_number(u2)
        u1 = max(u1, tiny(1.0_dp))
        z = sqrt(-2.0_dp * log(u1)) * cos(2.0_dp * rsdc_pi * u2)
    end function random_normal

    subroutine random_multivariate_normal(mean, covariance, x, ok)
        real(dp), intent(in) :: mean(:), covariance(:, :)
        real(dp), intent(out) :: x(:)
        logical, intent(out) :: ok
        real(dp), allocatable :: l(:, :), z(:)
        integer :: i, n
        n = size(mean)
        allocate(l(n, n), z(n))
        call cholesky_lower(covariance, l, ok)
        if (.not. ok) return
        do i = 1, n
            z(i) = random_normal()
        end do
        x = mean + matmul(l, z)
    end subroutine random_multivariate_normal

    integer function categorical_draw(probabilities) result(index)
        real(dp), intent(in) :: probabilities(:)
        real(dp) :: u, cumulative, total
        integer :: i
        total = sum(max(probabilities, 0.0_dp))
        if (total <= 0.0_dp) then
            index = 1
            return
        end if
        call random_number(u)
        u = u * total
        cumulative = 0.0_dp
        do i = 1, size(probabilities)
            cumulative = cumulative + max(probabilities(i), 0.0_dp)
            if (u <= cumulative) then
                index = i
                return
            end if
        end do
        index = size(probabilities)
    end function categorical_draw
end module rsdc_random
