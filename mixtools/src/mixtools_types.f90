! SPDX-License-Identifier: GPL-2.0-or-later
module mixtools_types
  use mixtools_kinds, only : dp
  use mixtools_status, only : MIXTOOLS_SUCCESS
  implicit none
  private
  public :: em_control, mixture_result, mv_mixture_result, gamma_mixture_result
  public :: multinomial_mixture_result, regression_mixture_result, semiparametric_result
  public :: reliability_mixture_result, bootstrap_result, model_selection_result
  public :: mcmc_result, constraint_result

  type :: em_control
    real(dp) :: tolerance = 1.0e-8_dp
    integer :: max_iterations = 1000
    real(dp) :: minimum_scale = 1.0e-8_dp
    real(dp) :: ridge = 1.0e-8_dp
    integer :: seed = 12345
    logical :: verbose = .false.
  end type em_control

  type :: mixture_result
    real(dp), allocatable :: lambda(:)
    real(dp), allocatable :: mu(:)
    real(dp), allocatable :: sigma(:)
    real(dp), allocatable :: posterior(:,:)
    real(dp), allocatable :: loglik_history(:)
    real(dp) :: loglik = -huge(1.0_dp)
    integer :: iterations = 0
    integer :: status = MIXTOOLS_SUCCESS
    logical :: converged = .false.
  end type mixture_result

  type :: mv_mixture_result
    real(dp), allocatable :: lambda(:)
    real(dp), allocatable :: mu(:,:)
    real(dp), allocatable :: sigma(:,:,:)
    real(dp), allocatable :: posterior(:,:)
    real(dp), allocatable :: loglik_history(:)
    real(dp) :: loglik = -huge(1.0_dp)
    integer :: iterations = 0
    integer :: status = MIXTOOLS_SUCCESS
    logical :: converged = .false.
  end type mv_mixture_result

  type :: gamma_mixture_result
    real(dp), allocatable :: lambda(:)
    real(dp), allocatable :: shape(:)
    real(dp), allocatable :: scale(:)
    real(dp), allocatable :: posterior(:,:)
    real(dp), allocatable :: loglik_history(:)
    real(dp) :: loglik = -huge(1.0_dp)
    integer :: iterations = 0
    integer :: status = MIXTOOLS_SUCCESS
    logical :: converged = .false.
  end type gamma_mixture_result

  type :: multinomial_mixture_result
    real(dp), allocatable :: lambda(:)
    real(dp), allocatable :: theta(:,:)
    real(dp), allocatable :: posterior(:,:)
    real(dp), allocatable :: loglik_history(:)
    real(dp) :: loglik = -huge(1.0_dp)
    integer :: iterations = 0
    integer :: status = MIXTOOLS_SUCCESS
    logical :: converged = .false.
  end type multinomial_mixture_result

  type :: regression_mixture_result
    real(dp), allocatable :: lambda(:)
    real(dp), allocatable :: beta(:,:)
    real(dp), allocatable :: sigma(:)
    real(dp), allocatable :: posterior(:,:)
    real(dp), allocatable :: auxiliary(:,:)
    real(dp), allocatable :: loglik_history(:)
    real(dp) :: loglik = -huge(1.0_dp)
    integer :: iterations = 0
    integer :: status = MIXTOOLS_SUCCESS
    logical :: converged = .false.
  end type regression_mixture_result

  type :: semiparametric_result
    real(dp), allocatable :: lambda(:)
    real(dp), allocatable :: posterior(:,:)
    real(dp), allocatable :: grid(:)
    real(dp), allocatable :: density(:,:)
    real(dp), allocatable :: bandwidth(:)
    real(dp), allocatable :: location(:)
    real(dp), allocatable :: loglik_history(:)
    real(dp) :: loglik = -huge(1.0_dp)
    integer :: iterations = 0
    integer :: status = MIXTOOLS_SUCCESS
    logical :: converged = .false.
  end type semiparametric_result

  type :: reliability_mixture_result
    real(dp), allocatable :: lambda(:)
    real(dp), allocatable :: shape(:)
    real(dp), allocatable :: scale(:)
    real(dp), allocatable :: posterior(:,:)
    real(dp), allocatable :: loglik_history(:)
    real(dp) :: loglik = -huge(1.0_dp)
    integer :: iterations = 0
    integer :: status = MIXTOOLS_SUCCESS
    logical :: converged = .false.
  end type reliability_mixture_result


  type :: mcmc_result
    real(dp), allocatable :: lambda_draws(:,:)
    real(dp), allocatable :: beta_draws(:,:,:)
    real(dp), allocatable :: sigma_draws(:,:)
    integer, allocatable :: allocation_draws(:,:)
    integer :: status = MIXTOOLS_SUCCESS
  end type mcmc_result

  type :: constraint_result
    integer, allocatable :: category(:)
    real(dp), allocatable :: alpha(:)
    integer :: status = MIXTOOLS_SUCCESS
  end type constraint_result

  type :: bootstrap_result
    real(dp), allocatable :: statistic(:)
    real(dp) :: observed = 0.0_dp
    real(dp) :: p_value = 1.0_dp
    integer :: successful = 0
    integer :: status = MIXTOOLS_SUCCESS
  end type bootstrap_result

  type :: model_selection_result
    integer, allocatable :: components(:)
    real(dp), allocatable :: aic(:), bic(:), caic(:), icl(:), loglik(:)
    integer :: best_aic = 0
    integer :: best_bic = 0
    integer :: status = MIXTOOLS_SUCCESS
  end type model_selection_result
end module mixtools_types
