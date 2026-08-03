! SPDX-License-Identifier: GPL-2.0-or-later
module lifeinsurer_kinds
  use, intrinsic :: iso_fortran_env, only : int64
  implicit none
  private
  integer, parameter, public :: dp = kind(1.0d0)
  integer, parameter, public :: i8 = int64
  integer, parameter, public :: lir_success = 0
  integer, parameter, public :: lir_invalid_argument = 1
  integer, parameter, public :: lir_dimension_error = 2
  integer, parameter, public :: lir_zero_denominator = 3
  integer, parameter, public :: lir_nonfinite = 4
  integer, parameter, public :: tariff_annuity = 1
  integer, parameter, public :: tariff_term_fix = 2
  integer, parameter, public :: tariff_dread_disease = 3
  integer, parameter, public :: tariff_endowment = 4
  integer, parameter, public :: tariff_pure_endowment = 5
  integer, parameter, public :: tariff_whole_life = 6
  integer, parameter, public :: tariff_endowment_dread = 7
  integer, parameter, public :: payment_advance = 1
  integer, parameter, public :: payment_arrears = 2
end module lifeinsurer_kinds
