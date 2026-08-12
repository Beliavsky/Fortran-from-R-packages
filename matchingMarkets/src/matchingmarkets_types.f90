module matchingmarkets_types
   use matchingmarkets_kinds, only : dp
   implicit none
   private
   public :: assignment_result_t, eadam_result_t, stable_set_t, plp_result_t, khb_result_t

   type :: assignment_result_t
      integer, allocatable :: assignment(:)
      integer, allocatable :: pairs(:,:)
      integer, allocatable :: singles(:)
      integer, allocatable :: free_capacity(:)
      integer :: iterations = 0
   end type assignment_result_t

   type :: eadam_result_t
      integer, allocatable :: assignment(:)
      integer, allocatable :: pairs(:,:)
      integer, allocatable :: singles(:)
      integer, allocatable :: student_prefs(:,:)
      integer, allocatable :: interrupting_pairs(:,:)
      integer :: iterations = 0
      integer :: ea_iterations = 0
   end type eadam_result_t

   type :: stable_set_t
      integer, allocatable :: assignments(:,:)
      integer :: count = 0
   end type stable_set_t

   type :: plp_result_t
      integer, allocatable :: assignment(:,:)
      integer, allocatable :: pairs(:,:)
      real(dp) :: objective = 0.0_dp
      integer :: status = 0
   end type plp_result_t

   type :: khb_result_t
      real(dp), allocatable :: p_value(:)
      real(dp), allocatable :: reduced_coef(:)
      real(dp), allocatable :: full_coef(:)
      real(dp), allocatable :: rescaled_coef(:)
      logical :: converged = .false.
   end type khb_result_t
end module matchingmarkets_types
