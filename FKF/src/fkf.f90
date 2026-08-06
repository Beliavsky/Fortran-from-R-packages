! SPDX-License-Identifier: GPL-2.0-or-later
module fkf_module
  use fkf_kinds, only : dp
  use fkf_types, only : fkf_model, fkf_result, fks_result, &
    fkf_success, fkf_invalid_input, fkf_non_pos_def
  use fkf_filter, only : fkf, kalman_filter, validate_model
  use fkf_smoother, only : fks, kalman_smooth
  implicit none
  public
end module fkf_module
