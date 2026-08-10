! Public facade for the modern Fortran translation of Rcsdp/CSDP.
! See LICENSE (CPL-1.0).
module rcsdp
   use rcsdp_kinds, only : dp
   use rcsdp_types
   use rcsdp_triplet
   use rcsdp_problem_mod
   use rcsdp_sparse_ops
   use rcsdp_fill_ops
   use rcsdp_solver
   use rcsdp_io
   implicit none
   public
end module rcsdp
