! SPDX-License-Identifier: GPL-2.0-or-later
module tsa_types
  use tsa_kinds, only : dp
  implicit none
  private

  public :: tsa_test_result, spectrum_result, ar_fit_result, arimax_result
  public :: tar_result, runs_result, outlier_result, bootstrap_result
  public :: transfer_spec, spectral_estimate, tar_multi_result

  type :: tsa_test_result
    real(dp) :: statistic = 0.0_dp
    real(dp) :: p_value = 1.0_dp
    integer :: df = 0
    integer :: lag = 0
    integer :: order = 0
    integer :: status = 0
  end type tsa_test_result

  type :: spectrum_result
    real(dp), allocatable :: frequency(:)
    real(dp), allocatable :: spectrum(:)
    integer :: status = 0
  end type spectrum_result

  type :: ar_fit_result
    integer :: order = 0
    logical :: intercept = .true.
    real(dp), allocatable :: coefficients(:)
    real(dp), allocatable :: residuals(:)
    real(dp), allocatable :: fitted(:)
    real(dp) :: variance = huge(1.0_dp)
    real(dp) :: aic = huge(1.0_dp)
    integer :: status = 0
  end type ar_fit_result

  type :: transfer_spec
    integer :: ar_order = 0
    integer :: ma_order = 0
  end type transfer_spec

  type :: spectral_estimate
    real(dp), allocatable :: frequency(:)
    real(dp), allocatable :: spectrum(:)
    real(dp), allocatable :: spectrum_matrix(:,:)
    real(dp), allocatable :: coherence(:,:), phase(:,:)
    real(dp), allocatable :: lower(:), upper(:)
    real(dp), allocatable :: degrees_freedom_series(:), taper_series(:)
    real(dp) :: degrees_freedom = 0.0_dp
    real(dp) :: bandwidth = 0.0_dp
    real(dp) :: taper = 0.0_dp
    integer :: n_used = 0, orig_n = 0, pad = 0, order = 0, n_series = 1
    logical :: detrend = .false., demean = .false.
    character(len=24) :: method = ''
    integer :: status = 0
  end type spectral_estimate

  type :: arimax_result
    integer :: p = 0, d = 0, q = 0
    integer :: seasonal_p = 0, seasonal_d = 0, seasonal_q = 0, period = 1
    logical :: include_mean = .true.
    real(dp), allocatable :: ar(:), ma(:), sar(:), sma(:), regression(:)
    real(dp), allocatable :: transfer(:), series(:)
    logical, allocatable :: estimated(:)
    integer :: n_xreg = 0, n_io = 0, n_transfer = 0
    integer :: n_cond = 0
    character(len=6) :: method = 'CSS'
    real(dp), allocatable :: coefficients(:), covariance(:,:)
    real(dp), allocatable :: residuals(:), fitted(:)
    real(dp) :: sigma2 = huge(1.0_dp)
    real(dp) :: loglik = -huge(1.0_dp)
    real(dp) :: aic = huge(1.0_dp)
    integer :: iterations = 0
    integer :: status = 0
  end type arimax_result

  type :: tar_result
    integer :: p1 = 0, p2 = 0, d = 1
    logical :: constant1 = .true., constant2 = .true.
    real(dp) :: threshold = 0.0_dp
    integer :: threshold_index = 0
    real(dp), allocatable :: phi1(:), phi2(:)
    real(dp), allocatable :: residuals(:), standardized_residuals(:)
    real(dp), allocatable :: fitted(:)
    real(dp) :: rms1 = huge(1.0_dp), rms2 = huge(1.0_dp)
    integer :: n1 = 0, n2 = 0
    real(dp) :: aic = huge(1.0_dp)
    real(dp) :: transform_mean = 0.0_dp, transform_sd = 1.0_dp
    logical :: centered = .false., standardized = .false.
    character(len=5) :: transform = 'no'
    integer :: status = 0
  end type tar_result


  type :: tar_multi_result
    integer :: p1 = 0, p2 = 0, d = 1, nseries = 0
    logical :: constant1 = .true., constant2 = .true.
    real(dp) :: threshold = 0.0_dp
    integer :: threshold_index = 0
    real(dp), allocatable :: coefficients1(:), coefficients2(:)
    real(dp), allocatable :: rms1(:), rms2(:)
    integer, allocatable :: n1(:), n2(:)
    real(dp), allocatable :: residuals(:,:), standardized_residuals(:,:), fitted(:,:)
    real(dp), allocatable :: transform_mean(:), transform_sd(:)
    real(dp) :: aic = huge(1.0_dp)
    real(dp), allocatable :: aic_series(:)
    logical :: centered = .false., standardized = .false.
    character(len=5) :: transform = 'no'
    integer :: status = 0
  end type tar_multi_result

  type :: runs_result
    integer :: observed_runs = 0
    real(dp) :: expected_runs = 0.0_dp
    integer :: n1 = 0, n2 = 0
    real(dp) :: threshold = 0.0_dp
    real(dp) :: p_value = -1.0_dp
  end type runs_result

  type :: outlier_result
    integer, allocatable :: index(:)
    real(dp), allocatable :: statistic(:)
    real(dp) :: cutoff = 0.0_dp
  end type outlier_result

  type :: bootstrap_result
    real(dp), allocatable :: coefficients(:,:)
    real(dp), allocatable :: sigma2(:)
    integer :: successful = 0
  end type bootstrap_result
end module tsa_types
