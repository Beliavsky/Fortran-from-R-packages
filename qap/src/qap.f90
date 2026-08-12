module qap
   use qap_kinds, only : dp, i64
   use qap_types, only : qap_control_t, qap_result_t, qap_problem_t
   use qap_core, only : qap_obj, qap_swap_delta, qap_is_permutation
   use qap_io, only : read_qaplib
   use qap_sa_mod, only : qap_solve
   implicit none
   private
   public :: dp, i64
   public :: qap_control_t, qap_result_t, qap_problem_t
   public :: qap_obj, qap_swap_delta, qap_is_permutation
   public :: read_qaplib, qap_solve
end module qap
