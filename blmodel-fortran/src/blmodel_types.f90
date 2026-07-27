! SPDX-License-Identifier: GPL-3.0-only
! Derived from BLModel 1.0.2, Copyright (C) 2017 Andrzej Palczewski and Jan Palczewski.
module blmodel_types
  use blmodel_kinds, only : dp
  implicit none
  private

  type, public :: moment_result
    real(dp), allocatable :: mean(:)
    real(dp), allocatable :: covariance(:,:)
    logical :: ok = .false.
    character(len=:), allocatable :: message
  end type moment_result

  type, public :: equilibrium_result
    real(dp), allocatable :: market_returns(:)
    real(dp), allocatable :: portfolio(:)
    logical :: ok = .false.
    character(len=:), allocatable :: message
  end type equilibrium_result

  type, public :: posterior_result
    real(dp), allocatable :: returns(:,:)
    real(dp), allocatable :: probabilities(:)
    real(dp), allocatable :: equilibrium_returns(:)
    real(dp), allocatable :: view_covariance(:,:)
    logical :: ok = .false.
    character(len=:), allocatable :: message
  end type posterior_result

end module blmodel_types
