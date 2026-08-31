! SPDX-License-Identifier: GPL-2.0-only
! Computational translation of CRAN mitools 2.4 by Thomas Lumley.
! Fortran translation and modifications: 2026-08-30.
module mitools_types
   use r_kinds, only : dp
   implicit none
   private

   integer, parameter, public :: mitools_success = 0
   integer, parameter, public :: mitools_invalid_shape = 1
   integer, parameter, public :: mitools_insufficient_imputations = 2
   integer, parameter, public :: mitools_invalid_probability = 3
   integer, parameter, public :: mitools_invalid_index = 4

   type, public :: mi_result
      real(dp), allocatable :: coefficients(:)
      real(dp), allocatable :: variance(:, :)
      integer :: nimp = 0
      real(dp), allocatable :: df(:)
      real(dp), allocatable :: missinfo(:)
   end type mi_result

   type, public :: imputation_list
      real(dp), allocatable :: values(:, :, :)
   end type imputation_list

end module mitools_types
