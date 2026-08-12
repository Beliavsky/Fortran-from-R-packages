! SPDX-License-Identifier: GPL-2.0-or-later
module nilde_types
   use nilde_kinds, only : i8, dp
   implicit none
   private

   type, public :: integer_solutions_t
      integer(i8), allocatable :: x(:,:)  ! variables x solutions
      integer :: nsol = 0
   contains
      procedure :: clear => integer_solutions_clear
   end type integer_solutions_t

   type, public :: knapsack_result_t
      integer(i8), allocatable :: x(:,:)
      integer :: nsol = 0
      real(dp) :: objective = -huge(1.0_dp)
      logical :: legacy_unbounded_all = .false.
   end type knapsack_result_t

   type, public :: bin_packing_result_t
      integer, allocatable :: assignment(:,:) ! item x optimal solution, canonical bin numbers
      integer(i8), allocatable :: bin_ineff(:,:) ! bin x optimal solution
      integer(i8), allocatable :: total_ineff(:)
      integer :: min_bins = 0
      integer :: nsol = 0
   end type bin_packing_result_t

   type, public :: tsp_result_t
      integer, allocatable :: tours(:,:) ! city position x optimal tour
      integer :: ntours = 0
      integer(i8) :: tour_length = -1_i8
      integer(i8) :: initial_lower_bound = -1_i8
      integer(i8) :: final_lower_bound = -1_i8
      integer(i8) :: upper_bound = -1_i8
      integer :: iterations = 0
      integer(i8) :: edge_subsets_tested = 0_i8
      integer(i8) :: degree_feasible = 0_i8
   end type tsp_result_t

contains

   subroutine integer_solutions_clear(self)
      class(integer_solutions_t), intent(inout) :: self
      if (allocated(self%x)) deallocate(self%x)
      self%nsol = 0
   end subroutine integer_solutions_clear

end module nilde_types
