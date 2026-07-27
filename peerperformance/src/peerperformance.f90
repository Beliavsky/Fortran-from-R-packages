! SPDX-License-Identifier: GPL-2.0-or-later
! Derived from PeerPerformance 2.4.0, copyright 2012-2023 David Ardia and Kris Boudt.
module peerperformance
  use peerperformance_kinds, only: dp
  use peerperformance_types, only: peer_control, test_result, screening_result, &
                                   rolling_result, valid_control
  use peerperformance_stats, only: sharpe, modified_sharpe, alpha_coefficients, &
                                   alpha_testing, sharpe_testing_asymptotic, &
                                   modified_sharpe_testing_asymptotic, &
                                   sharpe_difference, modified_sharpe_difference, &
                                   sharpe_standard_error, modified_sharpe_standard_error
  use peerperformance_bootstrap, only: bootstrap_indices, sharpe_testing_bootstrap, &
                                       modified_sharpe_testing_bootstrap, &
                                       sharpe_block_size, modified_sharpe_block_size
  use peerperformance_pi, only: adjust_pi, compute_pizero, optimal_lambda, compute_peer_ratios
  use peerperformance_screening, only: alpha_screening, sharpe_screening, &
                                       modified_sharpe_screening, target_peer_performance, &
                                       roll_screening, exposure_heterogeneity
  implicit none
  public
end module peerperformance
