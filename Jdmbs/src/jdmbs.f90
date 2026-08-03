! SPDX-License-Identifier: GPL-2.0-or-later
module jdmbs
   use jdmbs_kinds, only : dp, int64
   use jdmbs_status, only : jdmbs_success, jdmbs_invalid_argument, &
      jdmbs_nonfinite_input, jdmbs_numerical_warning
   use jdmbs_model, only : jdmbs_control, jdmbs_result, normal_bs, jdm_bs, jdm_new_bs
   implicit none
   public
end module jdmbs
