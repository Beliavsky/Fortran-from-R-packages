! SPDX-License-Identifier: MIT
module vctrs
   !! Public umbrella module for the restricted modern-Fortran vctrs API.
   use vctrs_types, only: dp, new_vctr, vctr_type, vctrs_character, &
      vctrs_data_frame, vctrs_integer, vctrs_logical, vctrs_real
   use vctrs_core, only: vec_any_missing, vec_cast, vec_c, vec_compare, &
      vec_copy_element, vec_detect_missing, vec_duplicate_any, vec_duplicate_detect, vec_equal, &
      vec_group_id, vec_identify_runs, vec_if_else, vec_in, vec_insert, vec_is, vec_match, &
      vec_init, vec_order, vec_ptype2, vec_ptype_common, vec_recycle, vec_recycle_common, &
      vec_rep, vec_rep_each, vec_size, vec_size_common, vec_slice, vec_sort, &
      vec_unique, vec_unique_count, vec_unique_loc
   use vctrs_frames, only: new_data_frame, validate_data_frame, vec_cbind, vec_rbind
   implicit none
   public
end module vctrs
