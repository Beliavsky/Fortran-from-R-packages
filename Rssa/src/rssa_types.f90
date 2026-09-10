! SPDX-License-Identifier: GPL-2.0-or-later
! Typed containers for the modern Fortran translation of Rssa 1.1.
module rssa_types
   use rssa_kinds, only : dp
   implicit none
   private
   integer, parameter, public :: rssa_success = 0
   integer, parameter, public :: rssa_invalid_input = -3001
   integer, parameter, public :: rssa_numerical_failure = -3002
   integer, parameter, public :: rssa_not_supported = -3003
   type, public :: hmat_type
      real(dp), allocatable :: values(:, :)
   end type hmat_type
   type, public :: hbhmat_type
      real(dp), allocatable :: values(:, :)
      integer, allocatable :: field_shape(:)
      integer, allocatable :: window_shape(:)
   end type hbhmat_type
   type, public :: tmat_type
      real(dp), allocatable :: values(:, :)
      real(dp), allocatable :: lags(:)
   end type tmat_type
   type, public :: ssa_result
      real(dp), allocatable :: series(:)
      real(dp), allocatable :: trajectory(:, :)
      real(dp), allocatable :: sigma(:)
      real(dp), allocatable :: u(:, :)
      real(dp), allocatable :: v(:, :)
      integer :: window = 0
      integer :: n_special_right = 0
      integer :: n_special_left = 0
      integer :: info = rssa_success
      logical :: toeplitz = .false.
      logical :: circular = .false.
   end type ssa_result
   type, public :: mssa_result
      real(dp), allocatable :: series(:, :)
      integer, allocatable :: lengths(:)
      real(dp), allocatable :: trajectory(:, :)
      real(dp), allocatable :: sigma(:)
      real(dp), allocatable :: u(:, :)
      real(dp), allocatable :: v(:, :)
      integer :: window = 0
      integer :: info = rssa_success
   end type mssa_result
   type, public :: ssa2d_result
      real(dp), allocatable :: field(:, :)
      real(dp), allocatable :: trajectory(:, :)
      real(dp), allocatable :: sigma(:)
      real(dp), allocatable :: u(:, :)
      real(dp), allocatable :: v(:, :)
      integer :: window(2) = 0
      integer :: info = rssa_success
   end type ssa2d_result
   type, public :: cssa_result
      complex(dp), allocatable :: series(:)
      complex(dp), allocatable :: trajectory(:, :)
      real(dp), allocatable :: sigma(:)
      complex(dp), allocatable :: u(:, :)
      complex(dp), allocatable :: v(:, :)
      integer :: window = 0
      integer :: info = rssa_success
   end type cssa_result
   type, public :: period_estimate
      complex(dp), allocatable :: roots(:)
      real(dp), allocatable :: periods(:)
      real(dp), allocatable :: frequencies(:)
      real(dp), allocatable :: rates(:)
      real(dp), allocatable :: moduli(:)
      integer :: info = rssa_success
   end type period_estimate
   type, public :: gap_summary
      integer :: n_missing = 0
      integer :: n_left = 0
      integer :: n_right = 0
      integer :: n_internal = 0
      integer :: n_dense = 0
      integer :: n_sparse = 0
      integer, allocatable :: missing_indices(:)
   end type gap_summary
   type, public :: grouping_result
      integer, allocatable :: group(:)
      real(dp), allocatable :: score(:)
      integer :: n_groups = 0
      integer :: info = rssa_success
   end type grouping_result
end module rssa_types
