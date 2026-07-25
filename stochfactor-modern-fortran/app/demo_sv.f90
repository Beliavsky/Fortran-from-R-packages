! SPDX-License-Identifier: GPL-2.0-or-later
! Copyright (C) 2013-2026 the original stochvol and factorstochvol authors
! and the modern Fortran translation contributors.
! This file is free software: you may redistribute it and/or modify it
! under the terms of the GNU General Public License as published by the
! Free Software Foundation, either version 2, or (at your option) any later version.
! This file is distributed without any warranty; see LICENSE for details.
program demo_sv
  use sv_kinds,only:dp
  use sv_rng,only:seed_rng
  use sv_types
  use sv_core
  implicit none
  type(sv_params)::p
  type(sv_sim_result)::sim
  type(sv_prior)::prior
  type(sv_mcmc_options)::opt
  type(sv_draws)::draws
  call seed_rng(2026);p%mu=-1.0_dp;p%phi=.96_dp;p%sigma=.2_dp;p%nu=8.0_dp
  call simulate_sv(250,p,sim)
  opt%draws=100;opt%burnin=50;opt%sample_nu=.true.;opt%store_latent=.true.;opt%store_tau=.true.
  call svsample(sim%y,p,prior,opt,draws)
  print '(a,3f10.4)','Posterior mean mu, phi, sigma: ',sum(draws%para(1,:))/draws%ndraws,sum(draws%para(2,:))/draws%ndraws,sum(draws%para(3,:))/draws%ndraws
  print '(a,f10.4)','Posterior mean nu: ',sum(draws%para(4,:))/draws%ndraws
end program demo_sv
