! SPDX-License-Identifier: GPL-3.0-only
module nmof_types
   use nmof_kinds, only: dp
   implicit none
   private
   integer, parameter, public :: nmof_ok = 0
   integer, parameter, public :: nmof_invalid_input = 1
   integer, parameter, public :: nmof_max_iter = 2
   integer, parameter, public :: nmof_linear_solve_failed = 3
   integer, parameter, public :: nmof_infeasible = 4
   integer, parameter, public :: nmof_not_bracketed = 5
   integer, parameter, public :: nmof_numerical_failure = 6

   type, public :: optimization_result
      real(dp), allocatable :: xbest(:)
      real(dp) :: ofvalue = huge(1.0_dp)
      real(dp), allocatable :: population_values(:)
      real(dp), allocatable :: history(:, :)
      integer :: iterations = 0
      integer :: status = nmof_ok
   end type optimization_result

   type, public :: binary_optimization_result
      logical, allocatable :: xbest(:)
      real(dp) :: ofvalue = huge(1.0_dp)
      real(dp), allocatable :: population_values(:)
      real(dp), allocatable :: history(:, :)
      integer :: iterations = 0
      integer :: status = nmof_ok
   end type binary_optimization_result

   type, public :: option_result
      real(dp) :: value = 0.0_dp
      real(dp) :: delta = 0.0_dp
      real(dp) :: gamma = 0.0_dp
      real(dp) :: theta = 0.0_dp
      real(dp) :: vega = 0.0_dp
      real(dp) :: rho = 0.0_dp
      real(dp) :: rho_div = 0.0_dp
      real(dp) :: dvega_dspot = 0.0_dp
      real(dp) :: dvega_dvol = 0.0_dp
      integer :: status = nmof_ok
   end type option_result

   type, public :: cppi_result
      real(dp), allocatable :: value(:), cushion(:), bond(:), floor(:)
      real(dp), allocatable :: exposure(:), units(:), spot(:)
      integer :: status = nmof_ok
   end type cppi_result

   type, public :: drawdown_summary
      real(dp) :: maximum = 0.0_dp
      real(dp) :: high = 0.0_dp
      real(dp) :: low = 0.0_dp
      integer :: high_position = 0
      integer :: low_position = 0
   end type drawdown_summary

   type, public :: quadrature_rule
      real(dp), allocatable :: nodes(:), weights(:)
      integer :: status = nmof_ok
   end type quadrature_rule

   type, public :: frontier_result
      real(dp), allocatable :: returns(:), volatility(:), portfolios(:, :)
      integer :: status = nmof_ok
   end type frontier_result

   type, public :: pbo_result
      real(dp) :: pbo = 0.0_dp
      real(dp), allocatable :: lambda(:), in_sample(:), out_of_sample(:)
      integer :: status = nmof_ok
   end type pbo_result

   type, public :: qtable_result
      real(dp), allocatable :: whiskers(:, :)
      real(dp), allocatable :: median(:), minimum(:), maximum(:)
      integer :: status = nmof_ok
   end type qtable_result

   type, public :: bond_return_result
      real(dp), allocatable :: returns(:), duration(:), convexity(:)
      integer :: status = nmof_ok
   end type bond_return_result
end module nmof_types
