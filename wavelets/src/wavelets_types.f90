! SPDX-License-Identifier: GPL-2.0-or-later
! Derived from the R package wavelets 0.3-0.2 by Eric Aldrich.
! Modern Fortran translation and modifications: 2026-09-10.
module wavelets_types
   use wavelets_kinds, only : dp
   implicit none
   private

   type, public :: wt_filter_type
      integer :: l = 0
      integer :: level = 1
      real(dp), allocatable :: h(:)
      real(dp), allocatable :: g(:)
      character(len=:), allocatable :: wt_class
      character(len=:), allocatable :: wt_name
      character(len=:), allocatable :: transform
   end type wt_filter_type

   type, public :: coefficient_level_type
      real(dp), allocatable :: values(:, :)
   end type coefficient_level_type

   type, public :: wavelet_transform_type
      type(coefficient_level_type), allocatable :: w(:)
      type(coefficient_level_type), allocatable :: v(:)
      type(wt_filter_type) :: filter
      integer :: level = 0
      integer, allocatable :: n_boundary(:)
      character(len=:), allocatable :: boundary
      real(dp), allocatable :: series(:, :)
      logical :: aligned = .false.
      logical :: coe = .false.
   end type wavelet_transform_type

   type, public :: mra_type
      type(coefficient_level_type), allocatable :: d(:)
      type(coefficient_level_type), allocatable :: s(:)
      type(wt_filter_type) :: filter
      integer :: level = 0
      character(len=:), allocatable :: boundary
      real(dp), allocatable :: series(:, :)
      character(len=:), allocatable :: method
   end type mra_type

end module wavelets_types
