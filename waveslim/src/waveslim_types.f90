! SPDX-License-Identifier: BSD-3-Clause
! Modern Fortran computational translation of waveslim.
module waveslim_types
  use waveslim_kinds, only : dp
  use waveslim_status, only : status_type
  implicit none
  private

  type, public :: real_vector
    real(dp), allocatable :: values(:)
  end type real_vector

  type, public :: complex_vector
    complex(dp), allocatable :: values(:)
  end type complex_vector

  type, public :: wavelet_filter_type
    character(len=16) :: name = ''
    real(dp), allocatable :: hpf(:)
    real(dp), allocatable :: lpf(:)
    real(dp), allocatable :: dhpf(:)
    real(dp), allocatable :: dlpf(:)
  contains
    procedure :: length => filter_length
  end type wavelet_filter_type

  type, public :: hilbert_filter_type
    character(len=16) :: name = ''
    real(dp), allocatable :: h0(:), g0(:), h1(:), g1(:)
  contains
    procedure :: length => hilbert_length
  end type hilbert_filter_type

  type, public :: wavelet_transform
    type(real_vector), allocatable :: detail(:)
    real(dp), allocatable :: smooth(:)
    character(len=16) :: wavelet = ''
    character(len=16) :: boundary = 'periodic'
    character(len=16) :: method = ''
    integer :: original_length = 0
    type(status_type) :: status
  contains
    procedure :: levels => transform_levels
  end type wavelet_transform

  type, public :: complex_wavelet_transform
    type(complex_vector), allocatable :: detail(:)
    complex(dp), allocatable :: smooth(:)
    character(len=16) :: wavelet = ''
    character(len=16) :: boundary = 'periodic'
    character(len=16) :: method = ''
    integer :: original_length = 0
    type(status_type) :: status
  contains
    procedure :: levels => complex_transform_levels
  end type complex_wavelet_transform

  type, public :: packet_level
    type(real_vector), allocatable :: node(:)
  end type packet_level

  type, public :: packet_transform
    type(packet_level), allocatable :: level(:)
    character(len=16) :: wavelet = ''
    character(len=16) :: boundary = 'periodic'
    character(len=16) :: method = ''
    integer :: original_length = 0
    type(status_type) :: status
  contains
    procedure :: levels => packet_levels
  end type packet_transform

  type, public :: mra_result
    type(real_vector), allocatable :: detail(:)
    real(dp), allocatable :: smooth(:)
    character(len=16) :: wavelet = ''
    character(len=16) :: boundary = 'periodic'
    character(len=16) :: method = ''
    type(status_type) :: status
  end type mra_result

  type, public :: wavelet_level_2d
    real(dp), allocatable :: lh(:,:), hl(:,:), hh(:,:)
  end type wavelet_level_2d

  type, public :: wavelet_transform_2d
    type(wavelet_level_2d), allocatable :: level(:)
    real(dp), allocatable :: smooth(:,:)
    character(len=16) :: wavelet = ''
    character(len=16) :: boundary = 'periodic'
    character(len=16) :: method = ''
    integer :: original_shape(2) = 0
    type(status_type) :: status
  end type wavelet_transform_2d

  type, public :: wavelet_level_3d
    real(dp), allocatable :: band(:,:,:,:)
  end type wavelet_level_3d

  type, public :: wavelet_transform_3d
    type(wavelet_level_3d), allocatable :: level(:)
    real(dp), allocatable :: smooth(:,:,:)
    character(len=16) :: wavelet = ''
    character(len=16) :: boundary = 'periodic'
    character(len=16) :: method = ''
    integer :: original_shape(3) = 0
    type(status_type) :: status
  end type wavelet_transform_3d

  type, public :: test_result
    real(dp) :: statistic = 0.0_dp
    real(dp) :: p_value = 1.0_dp
    integer :: degrees_freedom = 0
    type(status_type) :: status
  end type test_result

  type, public :: estimate_result
    real(dp), allocatable :: estimate(:)
    real(dp), allocatable :: covariance(:,:)
    real(dp) :: objective = 0.0_dp
    integer :: iterations = 0
    logical :: converged = .false.
    type(status_type) :: status
  end type estimate_result

contains
  integer function filter_length(self)
    class(wavelet_filter_type), intent(in) :: self
    if (allocated(self%lpf)) then
      filter_length = size(self%lpf)
    else
      filter_length = 0
    end if
  end function filter_length

  integer function hilbert_length(self)
    class(hilbert_filter_type), intent(in) :: self
    if (allocated(self%h0)) then
      hilbert_length = size(self%h0)
    else
      hilbert_length = 0
    end if
  end function hilbert_length

  integer function transform_levels(self)
    class(wavelet_transform), intent(in) :: self
    if (allocated(self%detail)) then
      transform_levels = size(self%detail)
    else
      transform_levels = 0
    end if
  end function transform_levels

  integer function complex_transform_levels(self)
    class(complex_wavelet_transform), intent(in) :: self
    if (allocated(self%detail)) then
      complex_transform_levels = size(self%detail)
    else
      complex_transform_levels = 0
    end if
  end function complex_transform_levels

  integer function packet_levels(self)
    class(packet_transform), intent(in) :: self
    if (allocated(self%level)) then
      packet_levels = size(self%level) - 1
    else
      packet_levels = 0
    end if
  end function packet_levels
end module waveslim_types
