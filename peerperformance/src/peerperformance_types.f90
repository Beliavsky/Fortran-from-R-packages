! SPDX-License-Identifier: GPL-2.0-or-later
! Derived from PeerPerformance 2.4.0, copyright 2012-2023 David Ardia and Kris Boudt.
module peerperformance_types
  use peerperformance_kinds, only: dp
  implicit none
  private

  public :: valid_control

  type, public :: peer_control
    integer :: test_type = 1
    integer :: ttype = 2
    logical :: hac = .false.
    integer :: n_boot = 499
    integer :: block_length = 1
    integer :: p_boot = 1
    integer :: min_obs = 10
    integer :: min_obs_pi = 1
    logical :: has_lambda = .false.
    real(dp) :: lambda = 0.5_dp
    logical :: screen_beta = .false.
    logical :: fast_adjust = .true.
    real(dp) :: gamma_pos = 0.4_dp
    real(dp) :: gamma_neg = 0.6_dp
    integer :: seed = 12345
  end type peer_control

  type, public :: test_result
    integer :: n = 0
    integer :: status = 0
    character(len=160) :: message = ''
    real(dp), allocatable :: estimate(:,:)
    real(dp), allocatable :: difference(:)
    real(dp), allocatable :: standard_error(:)
    real(dp), allocatable :: tstat(:)
    real(dp), allocatable :: pvalue(:)
  end type test_result

  type, public :: screening_result
    integer :: status = 0
    character(len=160) :: message = ''
    integer :: ncoef = 0
    integer :: n_focal = 0
    integer :: n_peer_group = 0
    logical :: cross_group = .false.
    integer, allocatable :: nobs(:)
    integer, allocatable :: npeer(:,:)
    real(dp), allocatable :: estimate(:,:)
    real(dp), allocatable :: difference(:,:,:)
    real(dp), allocatable :: standard_error(:,:,:)
    real(dp), allocatable :: tstat(:,:,:)
    real(dp), allocatable :: pvalue(:,:,:)
    real(dp), allocatable :: lambda(:,:)
    real(dp), allocatable :: pizero(:,:)
    real(dp), allocatable :: pipos(:,:)
    real(dp), allocatable :: pineg(:,:)
  end type screening_result

  type, public :: rolling_result
    integer :: status = 0
    character(len=160) :: message = ''
    integer :: nwindow = 0
    integer :: ncoef = 0
    integer, allocatable :: window(:)
    integer, allocatable :: end_index(:)
    real(dp), allocatable :: pizero(:,:)
    real(dp), allocatable :: pipos(:,:)
    real(dp), allocatable :: pineg(:,:)
    real(dp), allocatable :: heterogeneity(:,:)
  end type rolling_result

contains

  pure logical function valid_control(control) result(ok)
    type(peer_control), intent(in) :: control
    ok = (control%test_type == 1 .or. control%test_type == 2) .and. &
         (control%ttype == 1 .or. control%ttype == 2) .and. &
         (control%p_boot == 1 .or. control%p_boot == 2) .and. &
         control%n_boot >= 1 .and. control%block_length >= 0 .and. &
         control%min_obs >= 1 .and. control%min_obs_pi >= 0 .and. &
         control%gamma_pos > 0.0_dp .and. control%gamma_pos < 1.0_dp .and. &
         control%gamma_neg > 0.0_dp .and. control%gamma_neg < 1.0_dp .and. &
         (.not. control%has_lambda .or. &
          (control%lambda >= 0.0_dp .and. control%lambda < 1.0_dp))
  end function valid_control

end module peerperformance_types
