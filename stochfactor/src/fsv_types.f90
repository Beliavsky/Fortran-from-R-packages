! SPDX-License-Identifier: GPL-2.0-or-later
! Copyright (C) 2013-2026 the original stochvol and factorstochvol authors
! and the modern Fortran translation contributors.
! This file is free software: you may redistribute it and/or modify it
! under the terms of the GNU General Public License as published by the
! Free Software Foundation, either version 2, or (at your option) any later version.
! This file is distributed without any warranty; see LICENSE for details.
module fsv_types
  use sv_kinds, only : dp
  implicit none
  private
  public :: fsv_sim_result,fsv_options,fsv_draws,fsv_prediction
  type :: fsv_sim_result
    real(dp),allocatable::y(:,:),factors(:,:),loadings(:,:)
    real(dp),allocatable::h_idio(:,:),h_factor(:,:),h0_idio(:),h0_factor(:)
    real(dp),allocatable::idio_params(:,:),factor_params(:,:)
  end type fsv_sim_result
  type :: fsv_options
    integer::draws=500,burnin=200,thin=1,sv_sweeps=1
    logical::lower_triangular=.true.,normal_gamma=.false.,store_factors=.true.,store_latent=.true.
    real(dp)::loading_prior_sd=1.0_dp,global_shrinkage=1.0_dp
    real(dp)::ng_a=1.0_dp,ng_c=1.0_dp,ng_d=1.0_dp
  end type fsv_options
  type :: fsv_draws
    real(dp),allocatable::loadings(:,:,:),factors(:,:,:),latent(:,:,:),para(:,:,:)
    real(dp),allocatable::local_scale(:,:,:),global_shrinkage(:,:)
    integer::nobs=0,nseries=0,nfactors=0,ndraws=0
  end type fsv_draws
  type :: fsv_prediction
    real(dp),allocatable::h(:,:,:),cov(:,:,:,:),cor(:,:,:,:),precision(:,:,:,:)
    real(dp),allocatable::logdet_precision(:,:)
  end type fsv_prediction
end module fsv_types
