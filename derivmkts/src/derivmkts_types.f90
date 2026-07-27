! SPDX-License-Identifier: MIT
! Copyright (c) 2014-2022 Robert L. McDonald
module derivmkts_types
    use derivmkts_kinds, only: dp
    implicit none
    private
    public :: option_pair, perpetual_result, compound_result, greek_result
    public :: binomial_result, asian_mc_result, simulation_result, quincunx_result

    type :: option_pair
        real(dp) :: call = 0.0_dp
        real(dp) :: put = 0.0_dp
    end type option_pair

    type :: perpetual_result
        real(dp) :: price = 0.0_dp
        real(dp) :: barrier = 0.0_dp
    end type perpetual_result

    type :: compound_result
        real(dp) :: price = 0.0_dp
        real(dp) :: critical_spot = 0.0_dp
    end type compound_result

    type :: greek_result
        real(dp) :: premium = 0.0_dp
        real(dp) :: delta = 0.0_dp
        real(dp) :: gamma = 0.0_dp
        real(dp) :: vega = 0.0_dp
        real(dp) :: rho = 0.0_dp
        real(dp) :: theta = 0.0_dp
        real(dp) :: psi = 0.0_dp
        real(dp) :: elasticity = 0.0_dp
    end type greek_result

    type :: binomial_result
        real(dp) :: price = 0.0_dp
        real(dp) :: delta = 0.0_dp
        real(dp) :: gamma = 0.0_dp
        real(dp) :: theta = 0.0_dp
        real(dp) :: p = 0.0_dp
        real(dp) :: up = 0.0_dp
        real(dp) :: down = 0.0_dp
        real(dp) :: h = 0.0_dp
        logical :: valid = .false.
        real(dp), allocatable :: stock_tree(:, :)
        real(dp), allocatable :: option_tree(:, :)
        real(dp), allocatable :: probability_tree(:, :)
        real(dp), allocatable :: delta_tree(:, :)
        real(dp), allocatable :: bond_tree(:, :)
        logical, allocatable :: exercise_tree(:, :)
    end type binomial_result

    type :: asian_mc_result
        real(dp) :: avg_price_call = 0.0_dp
        real(dp) :: avg_price_put = 0.0_dp
        real(dp) :: avg_strike_call = 0.0_dp
        real(dp) :: avg_strike_put = 0.0_dp
        real(dp) :: vanilla_call = 0.0_dp
        real(dp) :: vanilla_put = 0.0_dp
        real(dp) :: sd_avg_price_call = 0.0_dp
        real(dp) :: sd_avg_price_put = 0.0_dp
        real(dp) :: sd_avg_strike_call = 0.0_dp
        real(dp) :: sd_avg_strike_put = 0.0_dp
        real(dp) :: sd_vanilla_call = 0.0_dp
        real(dp) :: sd_vanilla_put = 0.0_dp
        real(dp) :: exact_avg_price_call = 0.0_dp
        real(dp) :: exact_avg_price_put = 0.0_dp
        real(dp) :: exact_avg_strike_call = 0.0_dp
        real(dp) :: exact_avg_strike_put = 0.0_dp
        real(dp) :: beta = 0.0_dp
    end type asian_mc_result

    type :: simulation_result
        real(dp), allocatable :: price(:, :, :)
        integer, allocatable :: jumps(:, :, :)
        logical :: valid = .false.
    end type simulation_result

    type :: quincunx_result
        integer, allocatable :: counts(:)
        real(dp), allocatable :: expected(:)
    end type quincunx_result

end module derivmkts_types
