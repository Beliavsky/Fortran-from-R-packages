! SPDX-License-Identifier: Artistic-2.0
module mts
   !! Convenience umbrella module for the MTS Fortran computational API.
   use mts_kinds
   use mts_types
   use mts_linalg
   use mts_stats
   use mts_rng
   use mts_optimize
   use mts_var
   use mts_varma
   use mts_regression
   use mts_diagnostics
   use mts_volatility
   use mts_factor
   use mts_ecm_missing
   use mts_structural
   implicit none
   public
end module mts
