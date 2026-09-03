! Copyright (C) 2000 Robert Gray
! Modern Fortran translation maintained for Fortran-from-R-packages.
! SPDX-License-Identifier: GPL-2.0-or-later
module cmprsk
   use r_kinds, only : dp
   use cmprsk_status
   use cmprsk_cuminc, only : cuminc_curve, gray_test_result, cumulative_incidence, gray_test, &
                             timepoint_indices, curve_timepoints
   use cmprsk_api, only : cuminc_entry, cuminc_result, fit_cuminc, cuminc_timepoints
   use cmprsk_crr, only : crr_result, fit_crr, predict_crr
   use cmprsk_summary, only : crr_summary_result, summarize_crr
   implicit none
   public
end module cmprsk
