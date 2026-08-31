! SPDX-License-Identifier: GPL-2.0-or-later
! Upstream mitml 0.4-5 (2023-03-08), authored by Simon Grund,
! Alexander Robitzsch, and Oliver Luedtke; upstream license GPL (>= 2).
! Modern free-form Fortran translation for Fortran-from-R-packages.
! Derived from computational routines in R package mitml 0.4-5.
module mitml_types
   use r_kinds, only : dp
   implicit none
   private

   integer, parameter, public :: MITML_OK = 0
   integer, parameter, public :: MITML_ERR_DIMENSION = 1
   integer, parameter, public :: MITML_ERR_ARGUMENT = 2
   integer, parameter, public :: MITML_ERR_LINALG = 3
   integer, parameter, public :: MITML_ERR_NUMERIC = 4

   type, public :: pooled_estimates
      real(dp), allocatable :: estimate(:)
      real(dp), allocatable :: std_error(:)
      real(dp), allocatable :: t_value(:)
      real(dp), allocatable :: df(:)
      real(dp), allocatable :: p_value(:)
      real(dp), allocatable :: riv(:)
      real(dp), allocatable :: fmi(:)
      real(dp), allocatable :: ubar(:, :)
      real(dp), allocatable :: between(:, :)
      real(dp), allocatable :: total(:, :)
      integer :: m = 0
      integer :: status = MITML_OK
      character(len=:), allocatable :: message
   end type pooled_estimates

   type, public :: mi_test_result
      real(dp) :: f_value = 0.0_dp
      real(dp) :: df1 = 0.0_dp
      real(dp) :: df2 = 0.0_dp
      real(dp) :: p_value = 1.0_dp
      real(dp) :: riv = 0.0_dp
      integer :: status = MITML_OK
      character(len=:), allocatable :: message
   end type mi_test_result

   type, public :: multilevel_r2_result
      real(dp) :: rb1 = 0.0_dp
      real(dp) :: rb2 = 0.0_dp
      real(dp) :: sb = 0.0_dp
      real(dp) :: mvp = 0.0_dp
      logical :: has_reduction_measures = .false.
      integer :: status = MITML_OK
      character(len=:), allocatable :: message
   end type multilevel_r2_result

contains

end module mitml_types
