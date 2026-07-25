! SPDX-License-Identifier: GPL-2.0-or-later
! Copyright (C) 2013-2026 the original stochvol and factorstochvol authors
! and the modern Fortran translation contributors.
! This file is free software: you may redistribute it and/or modify it
! under the terms of the GNU General Public License as published by the
! Free Software Foundation, either version 2, or (at your option) any later version.
! This file is distributed without any warranty; see LICENSE for details.
program custom_restriction
  use sv_kinds, only : dp
  use sv_rng, only : seed_rng
  use fsv_types, only : fsv_sim_result, fsv_options, fsv_draws
  use fsv_core, only : simulate_fsv, fit_fsv
  implicit none
  real(dp) :: loadings(3,1), idio_params(3,3), factor_params(1,2)
  logical :: restriction(3,1)
  type(fsv_sim_result) :: sim
  type(fsv_options) :: options
  type(fsv_draws) :: draws

  call seed_rng(4123)
  loadings(:,1) = [1.0_dp, 0.0_dp, 0.6_dp]
  idio_params(:,1) = [-1.4_dp, -1.2_dp, -1.0_dp]
  idio_params(:,2) = [0.92_dp, 0.94_dp, 0.96_dp]
  idio_params(:,3) = [0.20_dp, 0.18_dp, 0.16_dp]
  factor_params(1,:) = [0.97_dp, 0.12_dp]
  restriction(:,1) = [.true., .false., .true.]

  call simulate_fsv(100, loadings, idio_params, factor_params, sim)
  options%draws = 20
  options%burnin = 10
  options%normal_gamma = .true.
  call fit_fsv(sim%y, 1, options, draws, restriction)

  print '(a,3f10.4)', 'Posterior mean loading column: ', &
    sum(draws%loadings(:,1,:), dim=2) / real(draws%ndraws, dp)
end program custom_restriction
