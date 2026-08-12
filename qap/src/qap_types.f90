module qap_types
   use qap_kinds, only : dp, i64
   implicit none
   private

   type, public :: qap_control_t
      integer :: rep = 1
      integer :: miter = -1
      real(dp) :: fiter = 1.1_dp
      real(dp) :: ft = 0.5_dp
      integer :: maxsteps = 50
      integer(i64) :: seed = 1234567_i64
      logical :: verbose = .false.
   end type qap_control_t

   type, public :: qap_result_t
      integer, allocatable :: permutation(:)
      real(dp) :: objective = huge(1.0_dp)
      integer :: best_rep = 0
      integer(i64) :: attempted_swaps = 0_i64
      integer(i64) :: accepted_swaps = 0_i64
      integer(i64) :: duplicate_trials = 0_i64
      integer(i64) :: cooling_steps = 0_i64
   end type qap_result_t

   type, public :: qap_problem_t
      real(dp), allocatable :: A(:,:)
      real(dp), allocatable :: B(:,:)
      integer, allocatable :: solution(:)
      real(dp) :: opt = huge(1.0_dp)
      logical :: has_solution = .false.
   end type qap_problem_t

end module qap_types
