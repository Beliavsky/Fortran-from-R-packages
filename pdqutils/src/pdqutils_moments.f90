! SPDX-License-Identifier: LGPL-3.0-or-later
module pdqutils_moments
    use pdqutils_kinds, only : dp
    use pdqutils_special, only : binomial_dp
    implicit none
    private
    public :: moment2cumulant, cumulant2moment

contains

    pure function moment2cumulant(moms) result(kappa)
        real(dp), intent(in) :: moms(:)
        real(dp), allocatable :: kappa(:)
        integer :: n, m
        allocate(kappa(size(moms)))
        kappa = moms
        do n = 2, size(moms)
            do m = 1, n-1
                kappa(n) = kappa(n) - binomial_dp(n-1,m-1)*kappa(m)*moms(n-m)
            end do
        end do
    end function moment2cumulant

    pure function cumulant2moment(kappa) result(moms)
        real(dp), intent(in) :: kappa(:)
        real(dp), allocatable :: moms(:)
        integer :: n, m
        allocate(moms(size(kappa)))
        moms = kappa
        do n = 2, size(kappa)
            do m = 1, n-1
                moms(n) = moms(n) + binomial_dp(n-1,m-1)*kappa(m)*moms(n-m)
            end do
        end do
    end function cumulant2moment

end module pdqutils_moments
