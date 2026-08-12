! SPDX-License-Identifier: MIT
module optmatch_types
   use optmatch_kinds, only : dp
   implicit none
   private

   type, public :: distance_spec
      real(dp), allocatable :: value(:, :)
      logical, allocatable :: allowed(:, :)
   contains
      procedure :: n_treatment
      procedure :: n_control
   end type distance_spec

   type, public :: sparse_distance
      integer :: n_treatment = 0
      integer :: n_control = 0
      integer, allocatable :: treatment(:)
      integer, allocatable :: control(:)
      real(dp), allocatable :: value(:)
   end type sparse_distance

   type, public :: match_result
      logical :: feasible = .false.
      real(dp) :: objective = 0.0_dp
      integer :: n_selected = 0
      integer, allocatable :: treatment_group(:)
      integer, allocatable :: control_group(:)
      logical, allocatable :: selected(:, :)
   end type match_result

   type, public :: stratum_summary
      integer :: n_strata = 0
      integer, allocatable :: treatments(:)
      integer, allocatable :: controls(:)
      integer, allocatable :: frequency(:)
      real(dp) :: effective_sample_size = 0.0_dp
   end type stratum_summary

contains

integer function n_treatment(self) result(n)
   class(distance_spec), intent(in) :: self
   if (allocated(self%value)) then
      n = size(self%value, 1)
   else
      n = 0
   end if
end function n_treatment

integer function n_control(self) result(n)
   class(distance_spec), intent(in) :: self
   if (allocated(self%value)) then
      n = size(self%value, 2)
   else
      n = 0
   end if
end function n_control

end module optmatch_types
