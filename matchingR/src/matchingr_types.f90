module matchingr_types
   use matchingr_kinds, only : dp
   implicit none
   private

   type, public :: marriage_result_t
      integer, allocatable :: proposals(:)       ! proposer -> reviewer; 0 means unmatched
      integer, allocatable :: engagements(:)     ! reviewer -> proposer; 0 means unmatched
      integer, allocatable :: single_proposers(:)
      integer, allocatable :: single_reviewers(:)
   end type marriage_result_t

   type, public :: college_result_t
      integer, allocatable :: matched_students(:)   ! student -> college; 0 means unmatched
      integer, allocatable :: matched_colleges(:,:) ! college x max(slots); 0 means vacant
      integer, allocatable :: slots(:)
      integer, allocatable :: unmatched_students(:)
      integer, allocatable :: unmatched_colleges(:) ! one entry per vacant slot
   end type college_result_t

   type, public :: roommate_result_t
      integer, allocatable :: matching(:) ! 0 means unmatched/no stable matching for that agent
      logical :: stable_exists = .false.
   end type roommate_result_t

end module matchingr_types
