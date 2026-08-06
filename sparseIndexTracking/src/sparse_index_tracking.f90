! sparseIndexTracking modern Fortran translation
! Copyright (C) 2026 OpenAI
! SPDX-License-Identifier: GPL-3.0-only

module sparse_index_tracking
   use sparse_index_tracking_kinds, only : dp
   use sparse_index_tracking_projection, only : project_capped_simplex, bisection
   use sparse_index_tracking_core, only : sparse_index_fit, fit_sparse_index_tracking, &
      sp_index_track, spIndexTrack, tracking_objective, parse_measure, &
      sit_success, sit_invalid_argument, sit_dimension_error, sit_infeasible_bounds, &
      sit_degenerate_data, sit_iteration_limit, sit_numerical_error, &
      measure_ete, measure_dr, measure_hete, measure_hdr
   implicit none
   public
end module sparse_index_tracking
