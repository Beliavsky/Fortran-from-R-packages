! SPDX-License-Identifier: GPL-3.0-only
module mixedind_types
  use mixedind_kinds, only : dp
  implicit none
  private

  integer, parameter, public :: mixedind_success = 0
  integer, parameter, public :: mixedind_invalid_argument = 1
  integer, parameter, public :: mixedind_allocation_error = 2
  integer, parameter, public :: mixedind_numerical_error = 3

  type, public :: prepared_data_result
    real(dp), allocatable :: values(:)
    real(dp), allocatable :: cdf(:)
    real(dp), allocatable :: pdf(:)
    integer :: status = mixedind_success
  end type prepared_data_result

  type, public :: pair_dependence_result
    real(dp) :: tau = 0.0_dp
    real(dp) :: rho = 0.0_dp
    real(dp) :: scale = 0.0_dp
    integer :: status = mixedind_success
  end type pair_dependence_result

  type, public :: dependence_result
    real(dp), allocatable :: tau(:,:)
    real(dp), allocatable :: rho(:,:)
    real(dp), allocatable :: p_tau(:,:)
    real(dp), allocatable :: p_rho(:,:)
    real(dp) :: lb_tau = 0.0_dp
    real(dp) :: lb_rho = 0.0_dp
    real(dp) :: p_lb_tau = 100.0_dp
    real(dp) :: p_lb_rho = 100.0_dp
    integer :: status = mixedind_success
  end type dependence_result

  type, public :: serial_dependence_result
    real(dp), allocatable :: tau(:)
    real(dp), allocatable :: rho(:)
    real(dp), allocatable :: p_tau(:)
    real(dp), allocatable :: p_rho(:)
    real(dp) :: scale = 0.0_dp
    real(dp) :: lb_tau = 0.0_dp
    real(dp) :: lb_rho = 0.0_dp
    real(dp) :: p_lb_tau = 100.0_dp
    real(dp) :: p_lb_rho = 100.0_dp
    integer :: status = mixedind_success
  end type serial_dependence_result

  type, public :: sn_result
    real(dp), allocatable :: stats(:)
    integer, allocatable :: cardinality(:)
    integer, allocatable :: subsets(:,:)
    real(dp), allocatable :: multiplier(:,:,:)
    real(dp), allocatable :: sn_multiplier(:,:)
    real(dp) :: sn = 0.0_dp
    integer :: status = mixedind_success
  end type sn_result

  type, public :: bootstrap_result
    real(dp), allocatable :: cvm(:)
    real(dp) :: sn = 0.0_dp
    integer :: status = mixedind_success
  end type bootstrap_result

  type, public :: copula_test_result
    real(dp), allocatable :: cvm(:)
    real(dp), allocatable :: p_cvm(:)
    integer, allocatable :: cardinality(:)
    integer, allocatable :: subsets(:,:)
    real(dp) :: sn = 0.0_dp
    real(dp) :: tn = 0.0_dp
    real(dp) :: tn2 = 0.0_dp
    real(dp) :: p_sn = 100.0_dp
    real(dp) :: p_tn = 100.0_dp
    real(dp) :: p_tn2 = 100.0_dp
    integer :: status = mixedind_success
  end type copula_test_result

  type, public :: moebius_result
    real(dp), allocatable :: spearman(:)
    real(dp), allocatable :: vdw(:)
    real(dp), allocatable :: savage(:)
    integer, allocatable :: cardinality(:)
    integer, allocatable :: subsets(:,:)
    real(dp) :: ln_spearman = 0.0_dp
    real(dp) :: ln_vdw = 0.0_dp
    real(dp) :: ln_savage = 0.0_dp
    real(dp) :: ln2_spearman = 0.0_dp
    real(dp) :: ln2_vdw = 0.0_dp
    real(dp) :: ln2_savage = 0.0_dp
    real(dp) :: p_ln_spearman = 100.0_dp
    real(dp) :: p_ln_vdw = 100.0_dp
    real(dp) :: p_ln_savage = 100.0_dp
    real(dp) :: p_ln2_spearman = 100.0_dp
    real(dp) :: p_ln2_vdw = 100.0_dp
    real(dp) :: p_ln2_savage = 100.0_dp
    integer :: status = mixedind_success
  end type moebius_result

end module mixedind_types
