! SPDX-License-Identifier: BSD-3-Clause
module metrics_kinds
    implicit none
    private

    integer, parameter, public :: dp = kind(1.0d0)

    integer, parameter, public :: metrics_success = 0
    integer, parameter, public :: metrics_invalid_size = 1
    integer, parameter, public :: metrics_invalid_argument = 2
    integer, parameter, public :: metrics_undefined = 3
end module metrics_kinds
