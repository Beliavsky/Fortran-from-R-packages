! SPDX-License-Identifier: GPL-2.0-only
! Public facade for the computational translation of CRAN mitools 2.4.
! Fortran translation and modifications: 2026-08-30.
module mitools
   use r_kinds, only : dp
   use mitools_types, only : imputation_list, mi_result
   use mitools_types, only : mitools_insufficient_imputations, mitools_invalid_index
   use mitools_types, only : mitools_invalid_probability, mitools_invalid_shape, mitools_success
   use mitools_combine, only : mi_combine, mi_confidence_intervals, mi_standard_errors, mi_summary
   use mitools_imputation, only : imputation_cbind, imputation_dimensions, imputation_get
   use mitools_imputation, only : imputation_list_from_array, imputation_rbind
   use mitools_pv, only : pv_materialize, pv_select
   implicit none
   public
end module mitools
