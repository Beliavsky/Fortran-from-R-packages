! SPDX-License-Identifier: GPL-3.0-only
module garchito_types
   use garchito_kinds, only : dp
   implicit none
   private

   integer, parameter, public :: garchito_success = 0
   integer, parameter, public :: garchito_max_iterations = 1
   integer, parameter, public :: garchito_invalid_input = 2
   integer, parameter, public :: garchito_numerical_failure = 3

   type, public :: garchito_control
      integer :: max_iterations = 5000
      integer :: max_evaluations = 100000
      real(dp) :: tolerance = 1.0e-9_dp
      real(dp) :: simplex_scale = 0.10_dp
      integer :: trace = 0
   end type garchito_control

   type, public :: garchito_result
      real(dp), allocatable :: coefficients(:)
      character(len=16), allocatable :: coefficient_names(:)
      real(dp), allocatable :: sigma(:)
      real(dp) :: pred = 0.0_dp
      real(dp) :: objective = huge(1.0_dp)
      integer :: convergence = garchito_invalid_input
      integer :: iterations = 0
      integer :: evaluations = 0
      character(len=160) :: message = 'model not fitted'
   end type garchito_result
end module garchito_types
