! SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only
module multiatsm_types
  use multiatsm_kinds, only : dp
  implicit none
  private

  integer, parameter, public :: VARX_UNCONSTRAINED = 0
  integer, parameter, public :: VARX_SPANNED_RESTRICTED = 1
  integer, parameter, public :: VARX_FACTOR_RESTRICTED = 2

  type, public :: var_model
    real(dp), allocatable :: intercept(:)
    real(dp), allocatable :: phi(:, :)
    real(dp), allocatable :: sigma(:, :)
    real(dp), allocatable :: residuals(:, :)
  end type var_model

  type, public :: varx_country_model
    real(dp), allocatable :: phi0(:)
    real(dp), allocatable :: phi1(:, :)
    real(dp), allocatable :: phi1_star(:, :)
    real(dp), allocatable :: phi_global(:, :)
    real(dp), allocatable :: phi0_star(:, :)
    real(dp), allocatable :: sigma(:, :)
    real(dp), allocatable :: residuals(:, :)
  end type varx_country_model

  type, public :: gvar_model
    type(varx_country_model), allocatable :: country(:)
    type(var_model) :: global_model
    real(dp), allocatable :: gy0(:, :)
    real(dp), allocatable :: f0(:)
    real(dp), allocatable :: f1(:, :)
    real(dp), allocatable :: sigma_y(:, :)
  end type gvar_model

  type, public :: jll_model
    real(dp), allocatable :: pi_matrix(:, :)
    real(dp), allocatable :: orthogonal_factors(:, :)
    real(dp), allocatable :: k0_e(:)
    real(dp), allocatable :: k1_e(:, :)
    real(dp), allocatable :: k0(:)
    real(dp), allocatable :: k1(:, :)
    real(dp), allocatable :: sigma_ortho(:, :)
    real(dp), allocatable :: sigma_nonortho(:, :)
    real(dp), allocatable :: chol_ortho(:, :)
    real(dp), allocatable :: chol_nonortho(:, :)
    logical, allocatable :: feedback_free(:, :)
    logical, allocatable :: chol_free(:, :)
  end type jll_model

  type, public :: affine_loadings
    real(dp), allocatable :: a(:)
    real(dp), allocatable :: b(:, :)
    real(dp), allocatable :: b_adjustment(:)
  end type affine_loadings

  type, public :: atsm_likelihood_result
    real(dp) :: negative_log_likelihood = huge(1.0_dp)
    real(dp) :: mean_negative_log_likelihood = huge(1.0_dp)
    real(dp), allocatable :: a_loadings(:)
    real(dp), allocatable :: b_loadings(:, :)
    real(dp), allocatable :: yield_error_variance(:)
    real(dp), allocatable :: factor_residuals(:, :)
    real(dp), allocatable :: yield_residuals(:, :)
  end type atsm_likelihood_result

  type, public :: response_result
    real(dp), allocatable :: factors(:, :, :)
    real(dp), allocatable :: yields(:, :, :)
  end type response_result

  type, public :: variance_decomposition_result
    real(dp), allocatable :: factors(:, :, :)
    real(dp), allocatable :: yields(:, :, :)
  end type variance_decomposition_result

  type, public :: forecast_result
    real(dp), allocatable :: factors(:, :)
    real(dp), allocatable :: yields(:, :)
  end type forecast_result

  type, public :: bootstrap_result
    real(dp), allocatable :: phi_draws(:, :, :)
    real(dp), allocatable :: sigma_draws(:, :, :)
    real(dp), allocatable :: lower(:, :)
    real(dp), allocatable :: median(:, :)
    real(dp), allocatable :: upper(:, :)
  end type bootstrap_result

end module multiatsm_types
