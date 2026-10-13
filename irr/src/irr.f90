! SPDX-License-Identifier: GPL-2.0-or-later
module irr
   use irr_common, only: dp, irr_result, irr_ok, irr_invalid, irr_degenerate, irr_resource
   use irr_agreement, only: agree, maxwell, meancor, meanrho, kendall, finn, robinson, iota, kripp_alpha
   use irr_kappa, only: kappa2, kappa2_table, kappam_fleiss, kappam_light
   use irr_marginal, only: bhapkar, bhapkar_table, stuart_maxwell_mh, stuart_maxwell_ratings, rater_bias, rater_bias_ratings
   use irr_icc, only: icc
   use irr_designs, only: n_cohen_kappa, n2_cohen_kappa, rel_inter_intra, rel_inter_intra_result
   implicit none
   public
end module
