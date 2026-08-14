! SPDX-License-Identifier: GPL-2.0-or-later
module lsei_types
   use lsei_kinds, only : dp
   implicit none
   private
   integer, parameter, public :: LSEI_SUCCESS = 1
   integer, parameter, public :: LSEI_BAD_DIMENSIONS = 2
   integer, parameter, public :: LSEI_ITERATION_LIMIT = 3
   integer, parameter, public :: LSEI_INFEASIBLE = 4
   integer, parameter, public :: LSEI_NUMERICAL = 5

   type, public :: ls_result
      real(dp), allocatable :: x(:)
      real(dp), allocatable :: residuals(:)
      real(dp), allocatable :: dual(:)
      integer, allocatable :: index(:)
      real(dp) :: rnorm = 0.0_dp
      real(dp) :: objective = 0.0_dp
      integer :: mode = LSEI_BAD_DIMENSIONS
      integer :: iterations = 0
      integer :: rank = 0
      integer :: k = 0
   contains
      procedure :: succeeded
   end type ls_result

   type, public :: hfti_result
      real(dp), allocatable :: x(:,:)
      real(dp), allocatable :: transformed_b(:,:)
      integer, allocatable :: pivot(:)
      real(dp), allocatable :: rnorm(:)
      integer :: krank = 0
      integer :: mode = LSEI_BAD_DIMENSIONS
   end type hfti_result
contains
   pure logical function succeeded(self)
      class(ls_result), intent(in) :: self
      succeeded = self%mode == LSEI_SUCCESS
   end function succeeded
end module lsei_types
