! SPDX-License-Identifier: GPL-2.0-only
! Copyright (C) 2005 Antonio Fabio Di Narzo
! Copyright (C) 2026 Modern Fortran translation contributors
! This file is part of a translation of tseriesChaos and is distributed
! under the GNU General Public License version 2 only.
module tserieschaos
  use chaos_kinds, only : dp
  use chaos_embedding, only : delay_embed, delay_embed_lags, delay_embed_matrix
  use chaos_systems, only : rhs_proc, observation_proc, lorenz_rhs, rossler_rhs, duffing_rhs, &
    integrate_rk4, simulate_observed, simulate_lorenz, simulate_rossler, simulate_duffing
  use chaos_metrics, only : correlation_integral, correlation_dimension_curve, &
    average_mutual_information, recurrence_distance_matrix, space_time_separation
  use chaos_neighbors, only : false_nearest_fraction, false_nearest_curve, find_k_nearests, &
    follow_neighbor_points, lyapunov_stretching, lyapunov_linear_fit
  implicit none
  public
end module tserieschaos
