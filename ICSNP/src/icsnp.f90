! SPDX-License-Identifier: GPL-2.0-or-later
module icsnp
   use icsnp_kinds, only : dp
   use icsnp_status, only : icsnp_ok, icsnp_invalid_input, icsnp_singular, &
      icsnp_iteration_limit, icsnp_numerical_error, icsnp_status_message
   use icsnp_types, only : test_result, location_scatter_result, spatial_sign_result
   use icsnp_pairs, only : pair_diff, pair_sum, pair_prod, spatial_ranks, signed_ranks
   use icsnp_estimators, only : spatial_median, spatial_sign, tyler_shape, &
      duembgen_shape, duembgen_shape_wt, symm_huber, symm_huber_wt, &
      HR_Mest, HP1_shape, hl_loc, vdw_loc
   use icsnp_tests, only : HotellingsT2, rank_ctest, rank_ctest_groups, &
      rank_ictest, ind_ctest, ind_ictest, HP_loc_test
   implicit none
   public
end module icsnp
