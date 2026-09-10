! SPDX-License-Identifier: MIT
! SPDX-FileComment: Public facade for the Fortran forcats translation.
module forcats
   !! Public umbrella module for the restricted modern-Fortran forcats API.
   use forcats_levels, only: lvls_expand, lvls_reorder, lvls_revalue, lvls_union
   use forcats_stats, only: fct_count, fct_lump, fct_lump_lowfreq, fct_lump_min, &
      fct_lump_n, fct_lump_prop, fct_reorder, fct_reorder2, first2, last2
   use forcats_transform, only: fct_c, fct_collapse, fct_cross, fct_drop, fct_expand, fct_explicit_na, &
      fct_infreq, fct_inorder, fct_inseq, fct_match, fct_na_level_to_value, &
      fct_na_value_to_level, fct_other, fct_recode, fct_relevel, fct_rev, &
      fct_shift, fct_unify, fct_unique
   use forcats_types, only: as_factor, dp, factor_count_type, factor_type, fct, new_factor
   implicit none
   public
end module forcats
