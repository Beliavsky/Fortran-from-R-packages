! SPDX-License-Identifier: GPL-3.0-only
! Derived from LSMonteCarlo 1.0 by Mikhail A. Beketov.
! Copyright (C) 2013 Mikhail A. Beketov.
module lsmc_surface
    use lsmc_kinds, only : dp
    use lsmc_pricing, only : amer_put_lsm, asian_amer_put_lsm, quanto_amer_put_lsm
    use lsmc_types, only : option_result, price_surface
    implicit none
    private

    interface american_put_lsmc_price_surface
        module procedure amer_put_lsm_price_surface
    end interface american_put_lsmc_price_surface

    interface asian_american_put_lsmc_price_surface
        module procedure asian_amer_put_lsm_price_surface
    end interface asian_american_put_lsmc_price_surface

    interface quanto_american_put_lsmc_price_surface
        module procedure quanto_amer_put_lsm_price_surface
    end interface quanto_american_put_lsmc_price_surface

    public :: amer_put_lsm_price_surface
    public :: asian_amer_put_lsm_price_surface
    public :: quanto_amer_put_lsm_price_surface
    public :: american_put_lsmc_price_surface
    public :: asian_american_put_lsmc_price_surface
    public :: quanto_american_put_lsmc_price_surface

contains

    function amer_put_lsm_price_surface(volatilities, strikes, spot, n, m, rate, dividend, maturity, seed) result(surface)
        real(dp), intent(in) :: volatilities(:)
        real(dp), intent(in) :: strikes(:)
        real(dp), intent(in), optional :: spot
        integer, intent(in), optional :: n
        integer, intent(in), optional :: m
        real(dp), intent(in), optional :: rate
        real(dp), intent(in), optional :: dividend
        real(dp), intent(in), optional :: maturity
        integer, intent(in), optional :: seed
        type(price_surface) :: surface
        type(option_result) :: fit
        real(dp) :: d
        real(dp) :: q
        real(dp) :: s
        real(dp) :: t
        integer :: base_seed
        integer :: i
        integer :: j
        integer :: np
        integer :: nt

        call resolve_surface_inputs(spot, n, m, rate, dividend, maturity, seed, s, np, nt, d, q, t, base_seed)
        call allocate_surface(surface, volatilities, strikes)
        do i = 1, size(volatilities)
            do j = 1, size(strikes)
                fit = amer_put_lsm(s, volatilities(i), np, nt, strikes(j), d, q, t, cell_seed(base_seed, i, j))
                surface%values(i, j) = fit%price
            end do
        end do
    end function amer_put_lsm_price_surface

    function asian_amer_put_lsm_price_surface(volatilities, strikes, spot, n, m, rate, dividend, maturity, seed) &
            result(surface)
        real(dp), intent(in) :: volatilities(:)
        real(dp), intent(in) :: strikes(:)
        real(dp), intent(in), optional :: spot
        integer, intent(in), optional :: n
        integer, intent(in), optional :: m
        real(dp), intent(in), optional :: rate
        real(dp), intent(in), optional :: dividend
        real(dp), intent(in), optional :: maturity
        integer, intent(in), optional :: seed
        type(price_surface) :: surface
        type(option_result) :: fit
        real(dp) :: d
        real(dp) :: q
        real(dp) :: s
        real(dp) :: t
        integer :: base_seed
        integer :: i
        integer :: j
        integer :: np
        integer :: nt

        call resolve_surface_inputs(spot, n, m, rate, dividend, maturity, seed, s, np, nt, d, q, t, base_seed)
        call allocate_surface(surface, volatilities, strikes)
        do i = 1, size(volatilities)
            do j = 1, size(strikes)
                fit = asian_amer_put_lsm(s, volatilities(i), np, nt, strikes(j), d, q, t, cell_seed(base_seed, i, j))
                surface%values(i, j) = fit%price
            end do
        end do
    end function asian_amer_put_lsm_price_surface

    function quanto_amer_put_lsm_price_surface(volatilities, strikes, spot, n, m, rate, dividend, maturity, spot2, &
            sigma2, rate2, dividend2, rho, seed) result(surface)
        real(dp), intent(in) :: volatilities(:)
        real(dp), intent(in) :: strikes(:)
        real(dp), intent(in), optional :: spot
        integer, intent(in), optional :: n
        integer, intent(in), optional :: m
        real(dp), intent(in), optional :: rate
        real(dp), intent(in), optional :: dividend
        real(dp), intent(in), optional :: maturity
        real(dp), intent(in), optional :: spot2
        real(dp), intent(in), optional :: sigma2
        real(dp), intent(in), optional :: rate2
        real(dp), intent(in), optional :: dividend2
        real(dp), intent(in), optional :: rho
        integer, intent(in), optional :: seed
        type(price_surface) :: surface
        type(option_result) :: fit
        real(dp) :: correlation
        real(dp) :: d
        real(dp) :: d2
        real(dp) :: q
        real(dp) :: q2
        real(dp) :: s
        real(dp) :: s2
        real(dp) :: t
        real(dp) :: v2
        integer :: base_seed
        integer :: i
        integer :: j
        integer :: np
        integer :: nt

        call resolve_surface_inputs(spot, n, m, rate, dividend, maturity, seed, s, np, nt, d, q, t, base_seed)
        s2 = 1.0_dp
        v2 = 0.2_dp
        d2 = 0.0_dp
        q2 = 0.0_dp
        correlation = 0.0_dp
        if (present(spot2)) s2 = spot2
        if (present(sigma2)) v2 = sigma2
        if (present(rate2)) d2 = rate2
        if (present(dividend2)) q2 = dividend2
        if (present(rho)) correlation = rho
        call allocate_surface(surface, volatilities, strikes)
        do i = 1, size(volatilities)
            do j = 1, size(strikes)
                fit = quanto_amer_put_lsm(s, volatilities(i), np, nt, strikes(j), d, q, t, s2, v2, d2, q2, &
                    correlation, cell_seed(base_seed, i, j))
                surface%values(i, j) = fit%price
            end do
        end do
    end function quanto_amer_put_lsm_price_surface

    subroutine allocate_surface(surface, volatilities, strikes)
        type(price_surface), intent(out) :: surface
        real(dp), intent(in) :: volatilities(:)
        real(dp), intent(in) :: strikes(:)

        if (size(volatilities) == 0 .or. size(strikes) == 0) then
            error stop "price surface: volatility and strike vectors must be nonempty"
        end if
        if (any(volatilities < 0.0_dp)) error stop "price surface: volatilities must be nonnegative"
        if (any(strikes <= 0.0_dp)) error stop "price surface: strikes must be positive"
        allocate(surface%volatilities(size(volatilities)), surface%strikes(size(strikes)))
        allocate(surface%values(size(volatilities), size(strikes)))
        surface%volatilities = volatilities
        surface%strikes = strikes
    end subroutine allocate_surface

    subroutine resolve_surface_inputs(spot, n, m, rate, dividend, maturity, seed, s, np, nt, d, q, t, base_seed)
        real(dp), intent(in), optional :: spot
        integer, intent(in), optional :: n
        integer, intent(in), optional :: m
        real(dp), intent(in), optional :: rate
        real(dp), intent(in), optional :: dividend
        real(dp), intent(in), optional :: maturity
        integer, intent(in), optional :: seed
        real(dp), intent(out) :: s
        integer, intent(out) :: np
        integer, intent(out) :: nt
        real(dp), intent(out) :: d
        real(dp), intent(out) :: q
        real(dp), intent(out) :: t
        integer, intent(out) :: base_seed

        s = 1.0_dp
        np = 1000
        nt = 365
        d = 0.06_dp
        q = 0.0_dp
        t = 1.0_dp
        base_seed = 271828
        if (present(spot)) s = spot
        if (present(n)) np = n
        if (present(m)) nt = m
        if (present(rate)) d = rate
        if (present(dividend)) q = dividend
        if (present(maturity)) t = maturity
        if (present(seed)) base_seed = seed
    end subroutine resolve_surface_inputs

    pure integer function cell_seed(base_seed, i, j) result(value)
        integer, intent(in) :: base_seed
        integer, intent(in) :: i
        integer, intent(in) :: j
        integer(kind=8) :: raw

        raw = abs(int(base_seed, kind=8)) + 1009_8 * int(i, kind=8) + 9176_8 * int(j, kind=8)
        value = int(modulo(raw, 2147483000_8)) + 1
    end function cell_seed

end module lsmc_surface
