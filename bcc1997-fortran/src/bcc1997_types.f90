! SPDX-License-Identifier: GPL-2.0-or-later
! Based on BCC1997 0.1.1, Copyright (C) 2017 Haoran Zhang.
module bcc1997_types
   use bcc1997_kinds, only : dp
   implicit none
   private

   type, public :: bcc_parameters
      real(dp) :: kappa_v = 0.0_dp
      real(dp) :: kappa_r = 0.0_dp
      real(dp) :: theta_v = 0.0_dp
      real(dp) :: theta_r = 0.0_dp
      real(dp) :: sigma_v = 1.0e-7_dp
      real(dp) :: sigma_r = 1.0e-7_dp
      real(dp) :: mu_j = 0.0_dp
      real(dp) :: sigma_j = 1.0e-7_dp
      real(dp) :: rho = 0.0_dp
      real(dp) :: lambda = 0.0_dp
      real(dp) :: spot = 100.0_dp
      real(dp) :: strike = 100.0_dp
      real(dp) :: variance0 = 0.04_dp
      real(dp) :: rate0 = 0.01_dp
      real(dp) :: maturity = 1.0_dp
   end type bcc_parameters

   type, public :: integration_settings
      real(dp) :: abs_tolerance = 1.0e-9_dp
      real(dp) :: rel_tolerance = 1.0e-8_dp
      real(dp) :: panel_width = 8.0_dp
      real(dp) :: maximum_upper_bound = 512.0_dp
      integer :: maximum_depth = 24
      integer :: minimum_panels = 4
      integer :: tail_panels = 4
   end type integration_settings

   type, public :: bcc_result
      real(dp) :: call = 0.0_dp
      real(dp) :: put = 0.0_dp
      real(dp) :: probability1 = 0.0_dp
      real(dp) :: probability2 = 0.0_dp
      real(dp) :: integral1 = 0.0_dp
      real(dp) :: integral2 = 0.0_dp
      real(dp) :: error1 = 0.0_dp
      real(dp) :: error2 = 0.0_dp
      real(dp) :: upper_bound1 = 0.0_dp
      real(dp) :: upper_bound2 = 0.0_dp
      integer :: evaluations1 = 0
      integer :: evaluations2 = 0
      logical :: converged = .false.
      integer :: status = 0
      character(len=160) :: message = ''
   end type bcc_result
end module bcc1997_types
