! SPDX-License-Identifier: GPL-3.0-only
! Derived from BLModel 1.0.2, Copyright (C) 2017 Andrzej Palczewski and Jan Palczewski.
module blmodel
  use blmodel_kinds, only : dp, pi
  use blmodel_types, only : moment_result, equilibrium_result, posterior_result
  use blmodel_utils, only : diag_of, make_diag
  use blmodel_distributions, only : view_density_interface, observ_normal, observ_powerexp, observ_student_t, observ_ts
  use blmodel_equilibrium, only : discrete_variance, equilibrium_mean, equilibrium_mean_elliptic
  use blmodel_posterior, only : post_distribution, bl_post_distribution
  implicit none
  private

  public :: dp, pi
  public :: moment_result, equilibrium_result, posterior_result
  public :: diag_of, make_diag
  public :: view_density_interface
  public :: observ_normal, observ_powerexp, observ_student_t, observ_ts
  public :: discrete_variance, equilibrium_mean, equilibrium_mean_elliptic
  public :: post_distribution, bl_post_distribution

end module blmodel
