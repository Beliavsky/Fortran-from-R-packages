! SPDX-License-Identifier: GPL-3.0-only
module rsdc_types
    use rsdc_kinds, only: dp
    implicit none
    private

    integer, parameter, public :: rsdc_const = 1
    integer, parameter, public :: rsdc_nox   = 2
    integer, parameter, public :: rsdc_tvtp  = 3

    type, public :: rsdc_filter_result
        logical :: ok = .false.
        real(dp) :: log_likelihood = -huge(1.0_dp)
        real(dp), allocatable :: filtered(:, :)
        real(dp), allocatable :: predicted(:, :)
        real(dp), allocatable :: smoothed(:, :)
        real(dp), allocatable :: loglik_t(:)
    end type rsdc_filter_result

    type, public :: rsdc_control
        integer :: seed = 123
        integer :: population_size = 0
        integer :: max_global_iterations = 250
        integer :: max_local_iterations = 250
        integer :: n_starts = 1
        real(dp) :: tolerance = 1.0e-7_dp
        real(dp) :: initial_step = 0.15_dp
        logical :: trace = .false.
        logical :: compute_vcov = .false.
        real(dp), allocatable :: start(:)
    end type rsdc_control

    type, public :: rsdc_model
        integer :: method = rsdc_const
        integer :: n_regimes = 1
        integer :: n_series = 0
        integer :: n_covariates = 0
        integer :: convergence = 1
        integer :: iterations = 0
        real(dp) :: log_likelihood = -huge(1.0_dp)
        real(dp) :: log_likelihood_oos = -huge(1.0_dp)
        real(dp), allocatable :: beta(:, :)
        real(dp), allocatable :: correlations(:, :)
        real(dp), allocatable :: transition_matrix(:, :)
        real(dp), allocatable :: covariance(:, :, :)
        real(dp), allocatable :: parameters(:)
        real(dp), allocatable :: vcov(:, :)
        real(dp), allocatable :: standard_errors(:)
        real(dp), allocatable :: average_x(:)
    end type rsdc_model

    type, public :: rsdc_forecast_result
        real(dp), allocatable :: regime_probabilities(:, :)
        real(dp), allocatable :: predicted_correlations(:, :)
        real(dp), allocatable :: covariance(:, :, :)
        real(dp) :: bic = huge(1.0_dp)
    end type rsdc_forecast_result

    type, public :: rsdc_simulation_result
        integer, allocatable :: states(:)
        real(dp), allocatable :: observations(:, :)
        real(dp), allocatable :: transition_matrices(:, :, :)
    end type rsdc_simulation_result

    type, public :: rsdc_portfolio_result
        real(dp), allocatable :: weights(:, :)
        real(dp), allocatable :: returns(:)
        real(dp), allocatable :: diversification_ratios(:)
        real(dp) :: realized_volatility = 0.0_dp
        real(dp) :: mean_return = 0.0_dp
        real(dp) :: mean_diversification = 0.0_dp
        integer :: n_fallback = 0
    end type rsdc_portfolio_result

    type, public :: rsdc_bootstrap_result
        real(dp), allocatable :: parameter_draws(:, :)
        real(dp), allocatable :: covariance(:, :)
        real(dp), allocatable :: standard_errors(:)
        real(dp), allocatable :: confidence_interval(:, :)
        real(dp) :: level = 0.95_dp
        integer :: n_success = 0
    end type rsdc_bootstrap_result


    type, public :: rsdc_starts_result
        real(dp), allocatable :: starts(:, :)
        real(dp), allocatable :: initial_log_likelihood(:)
        integer, allocatable :: regime_split(:)
        real(dp), allocatable :: shrinkage(:)
    end type rsdc_starts_result

    type, public :: rsdc_bands_result
        real(dp), allocatable :: fit(:, :)
        real(dp), allocatable :: lower(:, :)
        real(dp), allocatable :: upper(:, :)
        integer :: n_used = 0
        real(dp) :: level = 0.95_dp
    end type rsdc_bands_result

    type, public :: rsdc_diagnostics_result
        real(dp), allocatable :: stay_probability(:)
        real(dp), allocatable :: expected_duration(:)
        real(dp), allocatable :: ergodic_probability(:)
    end type rsdc_diagnostics_result
end module rsdc_types
