! SPDX-License-Identifier: GPL-2.0-or-later
module fkf_types
  use fkf_kinds, only : dp
  implicit none
  private

  integer, parameter, public :: fkf_success = 0
  integer, parameter, public :: fkf_invalid_input = 1
  integer, parameter, public :: fkf_non_pos_def = 2

  type, public :: fkf_model
    real(dp), allocatable :: a0(:)
    real(dp), allocatable :: p0(:, :)
    real(dp), allocatable :: dt(:, :)
    real(dp), allocatable :: ct(:, :)
    real(dp), allocatable :: tt(:, :, :)
    real(dp), allocatable :: zt(:, :, :)
    real(dp), allocatable :: hht(:, :, :)
    real(dp), allocatable :: ggt(:, :, :)
  end type fkf_model

  type, public :: fkf_result
    real(dp), allocatable :: att(:, :)
    real(dp), allocatable :: at(:, :)
    real(dp), allocatable :: ptt(:, :, :)
    real(dp), allocatable :: pt(:, :, :)
    real(dp), allocatable :: vt(:, :)
    real(dp), allocatable :: ft(:, :, :)
    real(dp), allocatable :: ftinv(:, :, :)
    real(dp), allocatable :: kt(:, :, :)
    real(dp) :: log_likelihood = 0.0_dp
    integer :: status = fkf_success
    integer :: failure_time = 0
    character(len=:), allocatable :: message
  end type fkf_result

  type, public :: fks_result
    real(dp), allocatable :: ahatt(:, :)
    real(dp), allocatable :: vt(:, :, :)
    integer :: status = fkf_success
    character(len=:), allocatable :: message
  end type fks_result

end module fkf_types
