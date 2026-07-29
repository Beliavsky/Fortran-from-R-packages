! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2018 Marius Hofert, Erik Hintz and Christiane Lemieux
module nvmix_kinds
  use iso_fortran_env, only : real64, int64
  implicit none
  private
  integer, parameter, public :: dp = real64
  integer, parameter, public :: i8 = int64
  real(dp), parameter, public :: pi = acos(-1.0_dp)
  real(dp), parameter, public :: log_two_pi = log(2.0_dp*pi)
  real(dp), parameter, public :: sqrt_two = sqrt(2.0_dp)
end module nvmix_kinds
