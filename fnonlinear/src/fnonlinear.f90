! SPDX-License-Identifier: GPL-2.0-or-later
! Copyright (C) 1988-1990 Blake LeBaron
! Copyright (C) 2000 Adrian Trapletti
! Copyright (C) 2005 Antonio Fabio Di Narzo
! Copyright (C) 2005-2026 Rmetrics contributors
! Copyright (C) 2026 Modern Fortran translation contributors
! This file is part of a modern Fortran translation of fNonlinear and is
! distributed under the GNU General Public License version 2 or later.
module fnonlinear
  use chaos_kinds, only : dp
  use chaos_embedding, only : delay_embed, delay_embed_lags, delay_embed_matrix
  use fnonlinear_rng, only : rng_state, rng_seed, rng_uniform, rng_normal, fill_uniform, fill_normal
  use fnonlinear_maps, only : ode_rhs, rk4_integrate_times, tent_sim, henon_sim, ikeda_sim, &
    logistic_sim, lorenz_rhs, rossler_rhs, lorenz_sim, rossler_sim, lorentz_sim, roessler_sim
  use fnonlinear_statistics, only : mutual_information_curve, false_nearest_neighbors, &
    recurrence_matrix, recurrence_distance_matrix, space_time_separation, &
    lyapunov_stretching, lyapunov_linear_fit, correlation_integral, &
    correlation_dimension_curve, find_k_nearests
  use fnonlinear_tests, only : bds_test_result, neural_test_result, runs_test_result, &
    generic_test_result, bds_test, white_neural_test, terasvirta_neural_test, &
    runs_test, ts_test
  implicit none
  public
end module fnonlinear
