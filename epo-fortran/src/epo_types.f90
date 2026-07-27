! SPDX-License-Identifier: MIT
! Copyright (c) 2023 Bernardo Reckziegel
module epo_types
  use epo_kinds, only : dp
  implicit none
  private

  integer, parameter, public :: epo_success = 0
  integer, parameter, public :: epo_invalid_input = 1
  integer, parameter, public :: epo_singular_matrix = 2
  integer, parameter, public :: epo_normalization_failure = 3

  type, public :: epo_result
    real(dp), allocatable :: weights(:)
    real(dp), allocatable :: covariance(:,:)
    real(dp), allocatable :: correlation(:,:)
    real(dp), allocatable :: shrunk_correlation(:,:)
    real(dp), allocatable :: shrunk_covariance(:,:)
    real(dp) :: gamma = 0.0_dp
    real(dp) :: weight_sum_before_normalization = 0.0_dp
    integer :: status = epo_success
    logical :: ok = .true.
    logical :: normalized = .false.
    logical :: endogenous = .false.
    character(len=:), allocatable :: method
    character(len=:), allocatable :: message
  end type epo_result

  public :: set_epo_error

contains

  subroutine set_epo_error(result, status, message)
    type(epo_result), intent(inout) :: result
    integer, intent(in) :: status
    character(len=*), intent(in) :: message

    result%ok = .false.
    result%status = status
    result%message = trim(message)
  end subroutine set_epo_error

end module epo_types
