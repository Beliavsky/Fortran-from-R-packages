! SPDX-License-Identifier: GPL-3.0-only
module indgen_types
  use indgen_kinds, only : dp
  implicit none
  private

  integer, parameter, public :: indgen_success = 0
  integer, parameter, public :: indgen_invalid_argument = 1
  integer, parameter, public :: indgen_numerical_error = 2

  type, public :: lag_test_result
    real(dp), allocatable :: stat(:)
    integer, allocatable :: lags(:,:)
    real(dp) :: aggregate = 0.0_dp
    real(dp) :: p_aggregate = 1.0_dp
    integer :: n = 0
    integer :: status = indgen_success
  end type lag_test_result

  type, public :: cvm_test_result
    real(dp), allocatable :: cvm(:)
    real(dp), allocatable :: p_cvm(:)
    integer, allocatable :: lags(:,:)
    real(dp) :: wstat = 0.0_dp
    real(dp) :: fstat = 0.0_dp
    real(dp) :: p_wstat = 1.0_dp
    real(dp) :: p_fstat = 1.0_dp
    integer :: n = 0
    integer :: status = indgen_success
  end type cvm_test_result

  type, public :: four_lag_test_result
    type(lag_test_result) :: xy
    type(lag_test_result) :: xz
    type(lag_test_result) :: yz
    type(lag_test_result) :: xyz
    real(dp) :: aggregate = 0.0_dp
    real(dp) :: p_aggregate = 1.0_dp
    integer :: status = indgen_success
  end type four_lag_test_result

  type, public :: cvm_three_result
    type(cvm_test_result) :: xy
    type(cvm_test_result) :: xz
    type(cvm_test_result) :: yz
    type(cvm_test_result) :: xyz
    real(dp) :: wstat = 0.0_dp
    real(dp) :: fstat = 0.0_dp
    real(dp) :: p_wstat = 1.0_dp
    real(dp) :: p_fstat = 1.0_dp
    integer :: status = indgen_success
  end type cvm_three_result

  type, public :: dependence_two_result
    type(lag_test_result) :: spearman
    type(lag_test_result) :: vdw
    type(lag_test_result) :: savage
    integer :: status = indgen_success
  end type dependence_two_result

  type, public :: dependence_three_result
    type(four_lag_test_result) :: spearman
    type(four_lag_test_result) :: vdw
    type(four_lag_test_result) :: savage
    integer :: status = indgen_success
  end type dependence_three_result

end module indgen_types
