! SPDX-License-Identifier: MIT
! Copyright (c) 2026 Charles Coverdale
module yc_types
  use yc_kinds, only : dp
  implicit none
  private

  type, public :: status_t
    logical :: ok = .true.
    character(len=256) :: message = ''
  end type status_t

  type, public :: curve_t
    logical :: ok = .true.
    character(len=256) :: message = ''
    character(len=24) :: rate_type = 'zero'
    character(len=24) :: method = 'observed'
    character(len=16) :: spline_method = ''
    real(dp) :: beta0 = 0.0_dp
    real(dp) :: beta1 = 0.0_dp
    real(dp) :: beta2 = 0.0_dp
    real(dp) :: beta3 = 0.0_dp
    real(dp) :: tau = 0.0_dp
    real(dp) :: tau1 = 0.0_dp
    real(dp) :: tau2 = 0.0_dp
    real(dp) :: objective = 0.0_dp
    integer :: iterations = 0
    real(dp), allocatable :: maturities(:)
    real(dp), allocatable :: rates(:)
    real(dp), allocatable :: fitted(:)
    real(dp), allocatable :: residuals(:)
    real(dp), allocatable :: spline_b(:)
    real(dp), allocatable :: spline_c(:)
    real(dp), allocatable :: spline_d(:)
  end type curve_t

  type, public :: series_t
    logical :: ok = .true.
    character(len=256) :: message = ''
    real(dp), allocatable :: x(:)
    real(dp), allocatable :: y(:)
  end type series_t

  type, public :: duration_result_t
    logical :: ok = .true.
    character(len=256) :: message = ''
    real(dp), allocatable :: maturity(:)
    real(dp), allocatable :: macaulay(:)
    real(dp), allocatable :: modified(:)
    real(dp), allocatable :: convexity(:)
  end type duration_result_t

  type, public :: bond_duration_result_t
    logical :: ok = .true.
    character(len=256) :: message = ''
    real(dp) :: macaulay_duration = 0.0_dp
    real(dp) :: modified_duration = 0.0_dp
    real(dp) :: convexity = 0.0_dp
    real(dp) :: price = 0.0_dp
  end type bond_duration_result_t

  type, public :: zspread_result_t
    logical :: ok = .true.
    character(len=256) :: message = ''
    real(dp) :: zspread = 0.0_dp
    real(dp) :: price = 0.0_dp
    real(dp) :: model_price = 0.0_dp
    integer :: iterations = 0
  end type zspread_result_t

  type, public :: carry_result_t
    logical :: ok = .true.
    character(len=256) :: message = ''
    real(dp), allocatable :: maturity(:)
    real(dp), allocatable :: carry(:)
    real(dp), allocatable :: rolldown(:)
    real(dp), allocatable :: total(:)
  end type carry_result_t

  type, public :: slope_result_t
    logical :: ok = .true.
    character(len=256) :: message = ''
    real(dp) :: spread_2s10s = 0.0_dp
    real(dp) :: spread_2s30s = 0.0_dp
    real(dp) :: spread_5s30s = 0.0_dp
    real(dp) :: spread_3m10y = 0.0_dp
    real(dp) :: butterfly_2s5s10s = 0.0_dp
  end type slope_result_t

  type, public :: factor_result_t
    logical :: ok = .true.
    character(len=256) :: message = ''
    real(dp) :: level = 0.0_dp
    real(dp) :: slope = 0.0_dp
    real(dp) :: curvature = 0.0_dp
  end type factor_result_t

  type, public :: pca_result_t
    logical :: ok = .true.
    character(len=256) :: message = ''
    integer :: n_components = 0
    real(dp), allocatable :: loadings(:,:)
    real(dp), allocatable :: scores(:,:)
    real(dp), allocatable :: variance_explained(:)
    real(dp), allocatable :: cumulative_variance(:)
    real(dp), allocatable :: sdev(:)
    real(dp), allocatable :: center(:)
    real(dp), allocatable :: scale(:)
  end type pca_result_t

end module yc_types
