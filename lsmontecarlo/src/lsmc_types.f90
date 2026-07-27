! SPDX-License-Identifier: GPL-3.0-only
! Derived from LSMonteCarlo 1.0 by Mikhail A. Beketov.
! Copyright (C) 2013 Mikhail A. Beketov.
module lsmc_types
    use lsmc_kinds, only : dp
    implicit none
    private

    type, public :: option_result
        real(dp) :: price = 0.0_dp
        real(dp) :: standard_error = 0.0_dp
        real(dp) :: spot = 1.0_dp
        real(dp) :: strike = 1.1_dp
        real(dp) :: sigma = 0.2_dp
        real(dp) :: rate = 0.06_dp
        real(dp) :: dividend = 0.0_dp
        real(dp) :: maturity = 1.0_dp
        real(dp) :: spot2 = 0.0_dp
        real(dp) :: sigma2 = 0.0_dp
        real(dp) :: rate2 = 0.0_dp
        real(dp) :: dividend2 = 0.0_dp
        real(dp) :: rho = 0.0_dp
        integer :: n_paths = 0
        integer :: effective_paths = 0
        integer :: n_steps = 0
        character(len=40) :: option_type = ""
        character(len=40) :: method = ""
    end type option_result

    type, public :: price_surface
        real(dp), allocatable :: volatilities(:)
        real(dp), allocatable :: strikes(:)
        real(dp), allocatable :: values(:, :)
    end type price_surface

    interface price
        module procedure option_result_price
    end interface price

    public :: price
    public :: surface_mean
    public :: surface_minimum
    public :: surface_maximum

contains

    pure function option_result_price(result) result(value)
        type(option_result), intent(in) :: result
        real(dp) :: value

        value = result%price
    end function option_result_price

    pure function surface_mean(surface) result(value)
        type(price_surface), intent(in) :: surface
        real(dp) :: value

        if (.not. allocated(surface%values) .or. size(surface%values) == 0) then
            value = 0.0_dp
        else
            value = sum(surface%values) / real(size(surface%values), dp)
        end if
    end function surface_mean

    pure function surface_minimum(surface) result(value)
        type(price_surface), intent(in) :: surface
        real(dp) :: value

        if (.not. allocated(surface%values) .or. size(surface%values) == 0) then
            value = 0.0_dp
        else
            value = minval(surface%values)
        end if
    end function surface_minimum

    pure function surface_maximum(surface) result(value)
        type(price_surface), intent(in) :: surface
        real(dp) :: value

        if (.not. allocated(surface%values) .or. size(surface%values) == 0) then
            value = 0.0_dp
        else
            value = maxval(surface%values)
        end if
    end function surface_maximum

end module lsmc_types
