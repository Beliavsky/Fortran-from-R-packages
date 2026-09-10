! SPDX-License-Identifier: GPL-2.0-only
!
! Modern Fortran translation derived from CRAN signal 1.8-1.
module signal_types
    use signal_kinds, only : dp
    implicit none
    private

    type, public :: arma_filter
        real(dp), allocatable :: b(:)
        real(dp), allocatable :: a(:)
    end type arma_filter

    type, public :: zpg_filter
        complex(dp), allocatable :: zero(:)
        complex(dp), allocatable :: pole(:)
        real(dp) :: gain = 1.0_dp
    end type zpg_filter

    type, public :: filter_order
        integer :: n = 0
        real(dp), allocatable :: wc(:)
        character(len=8) :: filter_type = 'low'
        real(dp) :: rp = 0.0_dp
        real(dp) :: rs = 0.0_dp
        real(dp) :: beta = 0.0_dp
    end type filter_order

    type, public :: fft_filter
        real(dp), allocatable :: b(:)
        integer :: n = 0
    end type fft_filter

    type, public :: median_filter
        integer :: n = 3
    end type median_filter

    type, public :: frequency_response
        complex(dp), allocatable :: h(:)
        real(dp), allocatable :: f(:)
    end type frequency_response

    type, public :: analog_response
        complex(dp), allocatable :: h(:)
        real(dp), allocatable :: w(:)
    end type analog_response

    type, public :: group_delay_response
        real(dp), allocatable :: gd(:)
        real(dp), allocatable :: w(:)
    end type group_delay_response

    type, public :: impulse_response_data
        real(dp), allocatable :: x(:)
        real(dp), allocatable :: t(:)
    end type impulse_response_data

    type, public :: sgolay_filter
        real(dp), allocatable :: coefficients(:, :)
        integer :: p = 0
        integer :: n = 0
        integer :: m = 0
        real(dp) :: ts = 1.0_dp
    end type sgolay_filter

    type, public :: spectrogram_data
        complex(dp), allocatable :: s(:, :)
        real(dp), allocatable :: f(:)
        real(dp), allocatable :: t(:)
    end type spectrogram_data

end module signal_types
