! SPDX-License-Identifier: MIT
! Copyright (c) 2023 Bernardo Reckziegel
module epo
  use epo_core, only : anchored_epo, epo_from_covariance, epo_optimize, &
    simple_epo
  use epo_kinds, only : dp
  use epo_statistics, only : covariance_to_correlation, sample_covariance
  use epo_types, only : epo_invalid_input, epo_normalization_failure, &
    epo_result, epo_singular_matrix, epo_success
  implicit none
  public

end module epo
