! SPDX-License-Identifier: GPL-3.0-only
module bb
  use bb_kinds, only: dp
  use bb_interfaces, only: bb_scalar_fn, bb_gradient_fn, bb_vector_fn, bb_projection_fn
  use bb_types, only: spg_control, spg_result, sane_control, sane_result, &
    bboptim_control, bbsolve_control, multistart_result, bb_success
  use bb_projection, only: project_box, project_linear
  use bb_spg, only: spg, spg_box, spg_projected, spg_linear
  use bb_nonlinear, only: sane, dfsane
  use bb_drivers, only: bboptim, bboptim_box, bboptim_projected, bboptim_linear, &
    bbsolve, multistart_solve, multistart_optimize, multistart_optimize_box, &
    multistart_optimize_projected, multistart_optimize_linear
  implicit none
  public
end module bb
