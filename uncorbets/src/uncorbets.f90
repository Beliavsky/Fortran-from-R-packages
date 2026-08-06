! SPDX-License-Identifier: MIT
module uncorbets
  use uncorbets_kinds, only : dp
  use uncorbets_types, only : status_type, torsion_result, effective_bets_result, &
      max_effective_bets_result, uncorbets_ok, uncorbets_invalid_input, &
      uncorbets_not_pos_semidefinite, uncorbets_singular_matrix, &
      uncorbets_no_convergence
  use uncorbets_core, only : sqrtm, torsion, torsion_pca, torsion_minimum, &
      effective_bets, effective_bets_gradient, max_effective_bets
  implicit none
  public
end module uncorbets
