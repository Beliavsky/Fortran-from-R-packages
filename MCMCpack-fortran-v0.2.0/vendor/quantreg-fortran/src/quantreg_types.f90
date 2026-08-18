! SPDX-License-Identifier: GPL-2.0-or-later
module quantreg_types
  use quantreg_kinds, only : dp
  implicit none
  private

  type, public :: rq_result
    real(dp), allocatable :: coefficients(:)
    real(dp), allocatable :: residuals(:)
    real(dp), allocatable :: dual(:)
    real(dp) :: tau = 0.5_dp
    integer :: iterations = 0
    integer :: corrector_steps = 0
    integer :: info = 0
  end type rq_result

  type, public :: rq_multi_result
    real(dp), allocatable :: coefficients(:,:)
    real(dp), allocatable :: tau(:)
    integer :: iterations = 0
    integer :: corrector_steps = 0
    integer :: info = 0
  end type rq_multi_result

  type, public :: nlrq_result
    real(dp), allocatable :: coefficients(:)
    real(dp), allocatable :: residuals(:)
    real(dp) :: objective = huge(1.0_dp)
    integer :: iterations = 0
    integer :: info = 0
  end type nlrq_result

  type, public :: lprq_result
    real(dp), allocatable :: x(:)
    real(dp), allocatable :: fitted(:)
    real(dp), allocatable :: derivative(:)
    integer :: info = 0
  end type lprq_result

end module quantreg_types
