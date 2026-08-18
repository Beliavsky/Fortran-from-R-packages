! SPDX-License-Identifier: GPL-2.0-or-later
module tolerance_types
  use tolerance_kinds, only : dp
  implicit none
  private

  type, public :: tolerance_interval
    real(dp) :: alpha = 0.05_dp
    real(dp) :: p = 0.99_dp
    real(dp) :: estimate = 0.0_dp
    real(dp) :: lower = 0.0_dp
    real(dp) :: upper = 0.0_dp
  end type tolerance_interval

  type, public :: discrete_tolerance_interval
    real(dp) :: alpha = 0.05_dp
    real(dp) :: p = 0.99_dp
    real(dp) :: estimate = 0.0_dp
    integer :: lower = 0
    integer :: upper = 0
    logical :: upper_infinite = .false.
  end type discrete_tolerance_interval

  type, public :: regression_band
    real(dp), allocatable :: fit(:)
    real(dp), allocatable :: lower(:)
    real(dp), allocatable :: upper(:)
  end type regression_band

  type, public :: mv_tolerance_region
    real(dp), allocatable :: center(:)
    real(dp), allocatable :: covariance(:,:)
    real(dp) :: radius2 = 0.0_dp
    real(dp) :: alpha = 0.05_dp
    real(dp) :: p = 0.99_dp
  end type mv_tolerance_region

end module tolerance_types
