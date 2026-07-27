! SPDX-License-Identifier: MIT
! Based on bidask 2.1.5, Copyright (c) 2024 Emanuele Guidotti.
module bidask_types
  use bidask_kinds, only: dp
  implicit none
  private

  type, public :: ohlc_data
    real(dp), allocatable :: open(:)
    real(dp), allocatable :: high(:)
    real(dp), allocatable :: low(:)
    real(dp), allocatable :: close(:)
  contains
    procedure :: size => ohlc_size
    procedure :: valid_dimensions => ohlc_valid_dimensions
  end type ohlc_data

  type, public :: spread_result
    character(len=16), allocatable :: method(:)
    real(dp), allocatable :: value(:)
    logical :: ok = .true.
    character(len=160) :: message = ''
  end type spread_result

  type, public :: spread_series_result
    character(len=16), allocatable :: method(:)
    real(dp), allocatable :: value(:, :)
    logical :: ok = .true.
    character(len=160) :: message = ''
  end type spread_series_result

contains

  integer function ohlc_size(self) result(n)
    class(ohlc_data), intent(in) :: self
    n = 0
    if (allocated(self%open)) n = size(self%open)
  end function ohlc_size

  logical function ohlc_valid_dimensions(self) result(ok)
    class(ohlc_data), intent(in) :: self
    ok = allocated(self%open) .and. allocated(self%high) .and. &
      allocated(self%low) .and. allocated(self%close)
    if (.not. ok) return
    ok = size(self%high) == size(self%open) .and. &
      size(self%low) == size(self%open) .and. &
      size(self%close) == size(self%open)
  end function ohlc_valid_dimensions

end module bidask_types
