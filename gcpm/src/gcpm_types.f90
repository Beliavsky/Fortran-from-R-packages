! SPDX-License-Identifier: GPL-2.0-only
!
! Modern Fortran translation of computational methods from GCPM 1.2.2.
! Original software copyright (C) 2015 Kevin Jakob and Dr. Matthias Fischer.
! Fortran translation copyright (C) 2026.
module gcpm_types
   use gcpm_kinds, only: dp, name_len, sector_name_len
   implicit none
   private

   integer, parameter, public :: default_bernoulli = 1
   integer, parameter, public :: default_poisson = 2
   integer, parameter, public :: link_crp = 1
   integer, parameter, public :: link_cm = 2
   integer, parameter, public :: model_analytical_crp = 1
   integer, parameter, public :: model_simulation = 2

   type, public :: credit_portfolio
      integer :: n_counterparties = 0
      integer :: n_sectors = 0
      integer, allocatable :: number(:)
      character(len=name_len), allocatable :: name(:)
      character(len=sector_name_len), allocatable :: sector_name(:)
      integer, allocatable :: default_kind(:)
      real(dp), allocatable :: ead(:)
      real(dp), allocatable :: lgd(:)
      real(dp), allocatable :: pd(:)
      real(dp), allocatable :: weight(:,:)
      real(dp), allocatable :: idiosyncratic(:)
      real(dp), allocatable :: potential_loss(:)
      integer, allocatable :: loss_multiple(:)
      real(dp), allocatable :: discretized_loss(:)
      real(dp), allocatable :: discretized_pd(:)
      real(dp) :: analytical_expected_loss = 0.0_dp
      real(dp) :: analytical_sd_diversifiable = 0.0_dp
      real(dp) :: analytical_sd_systematic = 0.0_dp
      real(dp) :: analytical_sd = 0.0_dp
   end type credit_portfolio

   type, public :: gcpm_model
      integer :: model_kind = model_analytical_crp
      integer :: link_kind = link_crp
      real(dp) :: loss_unit = 1.0_dp
      real(dp) :: alpha_max = 0.9999_dp
      real(dp) :: loss_threshold = huge(1.0_dp)
      integer :: max_stored_scenarios = 1000
      integer :: seed = 12345
      integer :: n_simulations = 0
      real(dp), allocatable :: sector_variance(:)
   end type gcpm_model

   type, public :: loss_distribution
      integer :: model_kind = model_analytical_crp
      real(dp), allocatable :: loss(:)
      real(dp), allocatable :: pdf(:)
      real(dp), allocatable :: cdf(:)
      real(dp), allocatable :: recursion_a(:)
      real(dp), allocatable :: recursion_b(:,:)
      real(dp) :: expected_loss = 0.0_dp
      real(dp) :: standard_deviation = 0.0_dp
      real(dp) :: reached_alpha = 0.0_dp
      real(dp), allocatable :: simulated_losses(:)
      real(dp), allocatable :: likelihood(:)
      integer, allocatable :: stored_scenario(:)
      real(dp), allocatable :: stored_counterparty_losses(:,:)
      logical :: contributions_available = .false.
   end type loss_distribution

   type, public :: risk_measures
      real(dp), allocatable :: alpha(:)
      real(dp), allocatable :: var(:)
      real(dp), allocatable :: economic_capital(:)
      real(dp), allocatable :: expected_shortfall(:)
   end type risk_measures

end module gcpm_types
