! SPDX-License-Identifier: GPL-2.0-or-later
module rootsolve
  use rootsolve_kinds, only : dp
  use rootsolve_types
  use rootsolve_derivatives, only : perturb_value, gradient, hessian, jacobian_full, jacobian_band
  use rootsolve_roots, only : multiroot, uniroot_all, brent_root
  use rootsolve_steady, only : stode, stodes, steady
  use rootsolve_runsteady, only : runsteady
  use rootsolve_pde, only : steady_1d, steady_2d, steady_3d, steady_band, multiroot_1d
  use rootsolve_sparse, only : build_grid_pattern, discover_pattern
  implicit none
  public
end module rootsolve
