! SPDX-License-Identifier: GPL-2.0-or-later
! Copyright (C) 2013-2026 the original stochvol and factorstochvol authors
! and the modern Fortran translation contributors.
! This file is free software: you may redistribute it and/or modify it
! under the terms of the GNU General Public License as published by the
! Free Software Foundation, either version 2, or (at your option) any later version.
! This file is distributed without any warranty; see LICENSE for details.
module sv_types
  use sv_kinds, only : dp
  implicit none
  private
  public :: sv_params, sv_prior, sv_mcmc_options, sv_draws, sv_sim_result, sv_prediction
  type :: sv_params
    real(dp) :: mu=-10.0_dp, phi=0.98_dp, sigma=0.2_dp
    real(dp) :: nu=huge(1.0_dp), rho=0.0_dp
    real(dp), allocatable :: beta(:)
  end type sv_params
  type :: sv_prior
    real(dp) :: mu_mean=0.0_dp, mu_sd=100.0_dp
    real(dp) :: phi_a=20.0_dp, phi_b=1.5_dp
    real(dp) :: sigma_shape=0.5_dp, sigma_rate=0.5_dp
    real(dp) :: nu_rate=0.1_dp
    real(dp) :: rho_sd=1.0_dp
    real(dp) :: beta_sd=10.0_dp
  end type sv_prior
  type :: sv_mcmc_options
    integer :: draws=1000,burnin=500,thin=1
    integer :: latent_sweeps=1
    logical :: sample_mu=.true.,sample_phi=.true.,sample_sigma=.true.
    logical :: sample_nu=.false.,sample_rho=.false.,sample_beta=.true.
    logical :: use_mixture=.true.,store_latent=.true.,store_tau=.false.
    real(dp) :: offset=1.0e-12_dp
    real(dp) :: step_mu=0.12_dp,step_phi=0.08_dp,step_logsigma=0.06_dp
    real(dp) :: step_rho=0.06_dp,step_lognu=0.08_dp,step_beta=0.04_dp
    real(dp) :: step_latent=0.35_dp
  end type sv_mcmc_options
  type :: sv_draws
    real(dp), allocatable :: para(:,:), beta(:,:), latent(:,:), latent0(:), tau(:,:)
    real(dp), allocatable :: data(:), design(:,:)
    real(dp), allocatable :: accept(:)
    integer :: nobs=0, ndraws=0
  end type sv_draws
  type :: sv_sim_result
    real(dp), allocatable :: y(:), latent(:), vol(:), tau(:)
    real(dp) :: latent0=0.0_dp,vol0=0.0_dp
    type(sv_params) :: params
  end type sv_sim_result
  type :: sv_prediction
    real(dp), allocatable :: y(:,:,:), latent(:,:,:), vola(:,:,:)
  end type sv_prediction
end module sv_types
