! SPDX-License-Identifier: GPL-3.0-only AND LicenseRef-ISDA-CDS-Standard-Model
module creditr_kinds
  implicit none
  private

  integer, parameter, public :: dp = kind(1.0d0)
  integer, parameter, public :: creditr_ok = 0
  integer, parameter, public :: creditr_invalid_input = 1
  integer, parameter, public :: creditr_no_bracket = 2
  integer, parameter, public :: creditr_singular = 3
  integer, parameter, public :: creditr_io_error = 4
  integer, parameter, public :: creditr_max_iter = 5
  public :: quiet_nan

contains

  pure real(kind=dp) function quiet_nan() result(x)
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
    x = ieee_value(0.0_dp, ieee_quiet_nan)
  end function quiet_nan

end module creditr_kinds
