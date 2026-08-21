! SPDX-License-Identifier: GPL-2.0-or-later
module tsa_kinds
  implicit none
  private
  public :: dp
  integer, parameter :: dp = kind(1.0d0)
end module tsa_kinds
